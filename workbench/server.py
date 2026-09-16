"""Read-only, loopback molecular structure and evidence workbench.

Python's standard library is the only runtime dependency. No solver is imported.
The HTTP API accepts registered artifact IDs, never arbitrary filesystem paths.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
import os
from pathlib import Path, PurePosixPath
import shlex
import stat
import threading
from urllib.parse import parse_qs, quote, urlsplit
from uuid import uuid4


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STATIC_ROOT = Path(__file__).resolve().parent / "static"
MAX_FILE_BYTES = 8 * 1024 * 1024
MAX_SNAPSHOT_BYTES = 64 * 1024 * 1024
MAX_ATOMS = 1000
MAX_POSES = 100


class EvidenceError(ValueError):
    """Malformed, unavailable, or out-of-bounds evidence."""


@dataclass(frozen=True)
class SafeRoot:
    """An explicitly configured read boundary; child symlinks are disallowed."""

    path: Path
    label: str = "project"
    _identity: tuple[int, int] | None = None

    def __post_init__(self):
        if self._identity is not None:
            # directory() has already walked this path through anchored file
            # descriptors. Do not resolve a replacement symlink here.
            object.__setattr__(self, "path", Path(self.path))
            return
        root = Path(self.path).resolve(strict=True)
        if not root.is_dir():
            raise EvidenceError("Configured evidence root must be a directory")
        object.__setattr__(self, "path", root)
        descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            info = os.fstat(descriptor)
            object.__setattr__(self, "_identity", (info.st_dev, info.st_ino))
        finally:
            os.close(descriptor)

    def _open_root(self):
        descriptor = os.open(self.path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        info = os.fstat(descriptor)
        if (info.st_dev, info.st_ino) != self._identity:
            os.close(descriptor)
            raise EvidenceError("Configured evidence directory was replaced; restart to select a new boundary")
        return descriptor

    @staticmethod
    def parts(relative):
        if not isinstance(relative, (str, Path, PurePosixPath)):
            raise EvidenceError("Artifact path must be a relative string")
        text = str(relative)
        path = PurePosixPath(text)
        if not text or "\\" in text or "\x00" in text or path.is_absolute() or ".." in path.parts:
            raise EvidenceError("Absolute paths, traversal, and ambiguous path separators are not permitted")
        if not path.parts or path == PurePosixPath("."):
            raise EvidenceError("An artifact filename is required")
        return path.parts

    def directory(self, relative):
        parts = self.parts(relative)
        descriptor = self._open_root()
        try:
            for part in parts:
                child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=descriptor)
                os.close(descriptor)
                descriptor = child
            info = os.fstat(descriptor)
            return SafeRoot(self.path.joinpath(*parts), f"{self.label}/{PurePosixPath(*parts)}", (info.st_dev, info.st_ino))
        except OSError as exc:
            if isinstance(exc, FileNotFoundError):
                raise
            raise EvidenceError("Campaign directory is unreadable or uses a forbidden symlink") from exc
        finally:
            os.close(descriptor)

    def read_bytes(self, relative):
        # Walk directory descriptors with O_NOFOLLOW so an imported path cannot
        # escape through symlinks, including a link swapped after validation.
        parts = self.parts(relative)
        directory_fd = self._open_root()
        try:
            for part in parts[:-1]:
                next_fd = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory_fd)
                os.close(directory_fd)
                directory_fd = next_fd
            descriptor = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory_fd)
            with os.fdopen(descriptor, "rb") as stream:
                info = os.fstat(stream.fileno())
                if not stat.S_ISREG(info.st_mode) or info.st_size > MAX_FILE_BYTES:
                    raise EvidenceError("Artifact must be a regular file no larger than 8 MiB")
                content = stream.read(MAX_FILE_BYTES + 1)
                if len(content) > MAX_FILE_BYTES:
                    raise EvidenceError("Artifact grew beyond the 8 MiB read limit")
                return content
        except OSError as exc:
            if isinstance(exc, FileNotFoundError):
                raise
            raise EvidenceError("Artifact is unreadable or uses a forbidden symlink") from exc
        finally:
            os.close(directory_fd)


def decode_json(content):
    def invalid_constant(value):
        raise EvidenceError(f"Nonfinite JSON number {value} is not valid evidence")

    try:
        record = json.loads(content, parse_constant=invalid_constant)
        if not isinstance(record, dict):
            raise EvidenceError("Evidence JSON must contain an object")
        json.dumps(record, allow_nan=False)
        return record
    except (UnicodeError, json.JSONDecodeError, ValueError, RecursionError) as exc:
        raise EvidenceError(f"Unreadable evidence JSON: {exc}") from exc


def parse_xyz(content):
    """Read the first XYZ frame, retaining exactly its supplied coordinate values."""
    try:
        lines = content.decode("utf-8").splitlines()
        count = int(lines[0].strip())
        if not 1 <= count <= MAX_ATOMS or len(lines) < count + 2:
            raise EvidenceError(f"XYZ must contain 1–{MAX_ATOMS} complete atoms")
        properties = next((token.split("=", 1)[1] for token in shlex.split(lines[1]) if token.startswith("Properties=")), None)
        species_column, position_column = 0, 1
        if properties:
            fields = properties.split(":")
            if len(fields) % 3:
                raise EvidenceError("Invalid extended-XYZ Properties declaration")
            species_column = position_column = None
            offset = 0
            for index in range(0, len(fields), 3):
                name, kind, width_text = fields[index:index + 3]
                width = int(width_text)
                if width < 1:
                    raise EvidenceError("Invalid extended-XYZ property width")
                if name == "species" and kind == "S" and width == 1:
                    species_column = offset
                if name == "pos" and kind == "R" and width == 3:
                    position_column = offset
                offset += width
            if species_column is None or position_column is None:
                raise EvidenceError("Extended XYZ must declare species:S:1 and pos:R:3")
        symbols, positions = [], []
        for line in lines[2:count + 2]:
            fields = shlex.split(line)
            symbol = fields[species_column]
            if len(symbol) > 2 or not symbol.isalpha() or not symbol[0].isupper() or (len(symbol) == 2 and not symbol[1].islower()):
                raise EvidenceError("Invalid element symbol in XYZ")
            position = [float(fields[position_column + axis]) for axis in range(3)]
            if not all(math.isfinite(value) for value in position):
                raise EvidenceError("XYZ coordinates must be finite")
            symbols.append(symbol)
            positions.append(position)
        return symbols, positions
    except (IndexError, UnicodeError, ValueError) as exc:
        if isinstance(exc, EvidenceError):
            raise
        raise EvidenceError(f"Invalid XYZ geometry: {exc}") from exc


def validated_indices(value, count, field):
    if not isinstance(value, list) or any(type(index) is not int or not 0 <= index < count for index in value):
        raise EvidenceError(f"{field} must contain in-range integer atom indices")
    if len(set(value)) != len(value):
        raise EvidenceError(f"{field} contains duplicate atom indices")
    return value


def visual_bonds(symbols, positions, metadata, state):
    nominal = metadata.get("product_bonds" if state == "final" else "reactant_bonds")
    if nominal:
        if not isinstance(nominal, list):
            raise EvidenceError("Nominal bonds must be a list")
        pairs = []
        for bond in nominal:
            if not isinstance(bond, list) or len(bond) != 3:
                raise EvidenceError("Nominal bonds must be [atom, atom, order] triples")
            validated_indices(bond[:2], len(symbols), "nominal bond")
            if type(bond[2]) is not int or bond[2] not in (1, 2, 3):
                raise EvidenceError("Invalid nominal bond order")
            pairs.append(bond[:2])
        return pairs, "Nominal design connectivity; visual guides, not computed electronic bond orders."
    # Covalent-radius proximity is only a drawing aid, never a chemical result.
    radii = {"H": .31, "He": .28, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57,
             "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ge": 1.20, "Br": 1.20}
    pairs = []
    for left in range(len(symbols)):
        for right in range(left + 1, len(symbols)):
            if symbols[left] not in radii or symbols[right] not in radii:
                continue
            distance = math.dist(positions[left], positions[right])
            if .1 < distance <= 1.2 * (radii[symbols[left]] + radii[symbols[right]]):
                pairs.append([left, right])
    return pairs, "Distance-inferred visual guides (1.2 × covalent radii); not electronic bond orders. Unsupported radii are omitted."


def completed_result_errors(result, controls, atom_count, fixed):
    """Check basic numerical evidence, without certifying the physical model."""
    errors = []

    def finite(value):
        return type(value) in (int, float) and math.isfinite(value)

    def structure(summary):
        if not isinstance(summary, dict):
            errors.append("Completed result is missing a structure summary.")
            return None
        energy, residual = summary.get("energy_ev"), summary.get("free_force_max_ev_per_angstrom")
        forces = summary.get("forces_ev_per_angstrom")
        diagnostics = summary.get("quantum_diagnostics", {})
        valid_forces = (isinstance(forces, list) and len(forces) == atom_count and
            all(isinstance(force, list) and len(force) == 3 and all(finite(value) for value in force) for force in forces))
        if (not finite(energy) or not finite(residual) or residual < 0 or not valid_forces or
                not isinstance(diagnostics, dict) or diagnostics.get("scf_converged") is not True or
                diagnostics.get("gradient_completed") is not True):
            errors.append("Completed result lacks finite energy/forces and converged electronic evidence.")
            return None
        measured = max((math.sqrt(sum(component**2 for component in force)) for index, force in enumerate(forces) if index not in fixed), default=0)
        if not math.isclose(residual, measured, rel_tol=1e-7, abs_tol=1e-9):
            errors.append("Reported free-force residual disagrees with recorded forces.")
        return residual

    stage = controls.get("stage")
    if stage in ("singlepoint", "relax"):
        residual = structure(result.get("structure"))
        if stage == "relax" and (result.get("geometry_converged") is not True or
                residual is None or residual > controls.get("fmax", 0)):
            errors.append("Completed relaxation lacks matching geometry convergence.")
    elif stage == "path":
        endpoints, images = result.get("endpoints"), result.get("images")
        if (result.get("endpoints_converged") is not True or result.get("neb_converged") is not True or
                not isinstance(endpoints, list) or len(endpoints) != 2 or
                not isinstance(images, list) or len(images) != controls.get("images")):
            errors.append("Completed path lacks the required endpoint/image convergence evidence.")
        else:
            for endpoint in endpoints:
                if not isinstance(endpoint, dict):
                    errors.append("Completed path endpoint must be a structure object.")
                    continue
                residual = structure(endpoint)
                topology = endpoint.get("topology_screen", {})
                if (endpoint.get("geometry_converged") is not True or residual is None or
                        residual > controls.get("fmax", 0) or endpoint.get("endpoint_identity_ok") is False or
                        not isinstance(topology, dict) or topology.get("preserved") is False):
                    errors.append("Completed path has invalid endpoint evidence.")
            for image in images:
                structure(image)
            relative = result.get("relative_energies_ev")
            if not errors:
                expected = [image["energy_ev"] - images[0]["energy_ev"] for image in images]
                if (not isinstance(relative, list) or len(relative) != len(expected) or
                        any(not finite(value) or not math.isclose(value, expected[index], abs_tol=1e-8, rel_tol=0)
                            for index, value in enumerate(relative))):
                    errors.append("Recorded path profile disagrees with image energies.")
    else:
        errors.append("Unsupported calculation stage in the campaign plan.")
    return errors


def quantum_settings_errors(settings):
    """Check the supported serialized controls, without importing a solver."""
    if not isinstance(settings, dict):
        return ["Quantum settings must be an object."]
    supported = {"charge", "spin", "xc", "basis", "dispersion", "grid_level", "conv_tol",
                 "max_cycle", "threads", "memory_mb", "density_fit", "scf_initial_guess"}
    errors = []
    for key, value in settings.items():
        valid = key in supported
        if key in {"charge", "spin", "grid_level", "max_cycle", "threads", "memory_mb"}:
            valid = type(value) is int
            if valid and key == "spin":
                valid = value >= 0
            if valid and key == "grid_level":
                valid = 0 <= value <= 9
            if valid and key in {"max_cycle", "threads", "memory_mb"}:
                valid = value > 0
        elif key in {"xc", "basis"}:
            valid = isinstance(value, str) and bool(value.strip())
        elif key == "dispersion":
            valid = value is None or (isinstance(value, str) and value in {"d3bj", "d3zero"})
        elif key == "conv_tol":
            valid = type(value) in (int, float) and math.isfinite(value) and 0 < value < 1
        elif key == "density_fit":
            valid = isinstance(value, bool)
        elif key == "scf_initial_guess":
            valid = isinstance(value, str) and value in {"minao", "atom", "1e", "huckel"}
        if not valid:
            errors.append(f"Unsupported quantum setting or value: {key}.")
    return errors


class CatalogBuilder:
    def __init__(self, project, campaigns):
        self.project, self.campaign_roots = project, campaigns
        self.raw = {}
        self.structures = {}
        self.total_bytes = 0
        self.catalog = {"default_structure_id": None, "structures": [], "calculations": [],
                        "references": [], "campaigns": [], "notices": [], "read_only": True}

    def source(self, root, relative, identifier):
        content = root.read_bytes(relative)
        self.total_bytes += len(content)
        if self.total_bytes > MAX_SNAPSHOT_BYTES:
            raise EvidenceError("Catalog exceeds the 64 MiB snapshot limit")
        mime = "application/json; charset=utf-8" if str(relative).endswith(".json") else "text/plain; charset=utf-8"
        self.raw[identifier] = (content, mime)
        return content, f"/api/evidence?id={quote(identifier)}"

    def record(self, root, relative, identifier):
        try:
            content, href = self.source(root, relative, identifier)
            return decode_json(content), href
        except FileNotFoundError:
            return {"status": "missing", "load_error": "Evidence file is missing"}, None
        except (EvidenceError, OSError) as exc:
            return {"status": "unreadable", "load_error": str(exc)}, None

    def structure(self, identifier, root, relative, design, label, kind, state="initial", recorded_result=None):
        try:
            content, href = self.source(root, relative, f"source-{identifier}")
            symbols, positions = parse_xyz(content)
            count = len(symbols)
            metadata = design.get("metadata", {})
            if not isinstance(metadata, dict):
                raise EvidenceError("Design metadata must be an object")
            for field in ("substrate_indices", "tool_indices"):
                if field in metadata:
                    validated_indices(metadata[field], count, field)
            fixed = validated_indices(design.get("fixed_indices", []), count, "fixed_indices")
            reaction = design.get("hydrogen_transfer") or {}
            if reaction:
                if not isinstance(reaction, dict) or any(key not in reaction for key in ("donor", "hydrogen", "acceptor")):
                    raise EvidenceError("Incomplete hydrogen-transfer identities")
                validated_indices([reaction[key] for key in ("donor", "hydrogen", "acceptor")], count, "hydrogen_transfer")
            bonds, bond_source = visual_bonds(symbols, positions, metadata, state)
            source_name = f"{root.label}/{relative}"
            record = {"id": identifier, "label": label, "kind": kind, "source": source_name,
                      "source_href": href, "source_sha256": hashlib.sha256(content).hexdigest(),
                      "symbols": symbols, "positions": positions, "atom_count": count,
                      "element_counts": dict(Counter(symbols)), "fixed_indices": fixed,
                      "metadata": metadata, "hydrogen_transfer": reaction, "bonds": bonds,
                      "bond_source": bond_source, "units": "angstrom", "frame": 0,
                      "recorded_result": recorded_result}
            self.structures[identifier] = record
            self.catalog["structures"].append({key: record[key] for key in ("id", "label", "kind", "source", "source_href", "atom_count")})
            if self.catalog["default_structure_id"] is None:
                self.catalog["default_structure_id"] = identifier
            return identifier
        except (FileNotFoundError, EvidenceError, OSError) as exc:
            self.catalog["notices"].append(f"{label}: {exc}")
            return None

    def build_archives(self):
        for short, folder, label in (
            ("direct", "h-abstraction-direct-initial", "Direct DFT · 53 atoms"),
            ("df", "h-abstraction-df-initial", "Density-fitting DFT · 53 atoms"),
        ):
            base = PurePosixPath("data/validation") / folder
            identifier = f"calculation-{short}"
            result, href = self.record(self.project, base / "result.json", identifier)
            interpretation, interpretation_href = self.record(self.project, base / "interpretation.json", f"interpretation-{short}")
            design = result.get("design", {})
            if not isinstance(design, dict):
                design = {}
            structure_id = self.structure(f"{short}-structure", self.project, base / "structure.extxyz", design,
                label, "calculation", recorded_result=result)
            final_id = self.structure(f"{short}-final", self.project, base / "input-final.extxyz", design,
                f"{label} · proposed final geometry", "endpoint", state="final")
            self.catalog["calculations"].append({**result, "id": identifier, "label": label,
                "source": str(base / "result.json"), "source_href": href, "structure_id": structure_id,
                "initial_structure_id": structure_id, "final_structure_id": final_id,
                "interpretation": interpretation, "interpretation_href": interpretation_href,
                "path_status": "not computed", "vibrations_status": "not computed"})
        comparison, href = self.record(self.project, "data/validation/h-abstraction-df-initial/comparison.json", "direct-df-comparison")
        self.catalog["direct_df_comparison"] = {**comparison, "source_href": href}
        for suffix, label in (("", "Historical CC reference · guess unrecorded"), ("-atom", "CC reference guess · atom")):
            base = PurePosixPath(f"data/validation/paired-ccpvdz{suffix}")
            identifier = f"reference-{suffix.lstrip('-') or 'minao'}"
            result, href = self.record(self.project, base / "method_comparison.json", identifier)
            annotation, annotation_href = self.record(self.project, base / "interpretation.json", f"{identifier}-annotation")
            annotation_status, annotation_errors = "not available", []
            if not annotation.get("load_error"):
                if annotation.get("schema_version") != 1 or annotation.get("kind") != "retrospective_interpretation":
                    annotation_errors.append("Unsupported retrospective annotation schema or kind.")
                source_bytes = self.raw.get(identifier, (b"", ""))[0]
                if (annotation.get("source") != "method_comparison.json" or result.get("load_error") or
                        annotation.get("source_sha256") != hashlib.sha256(source_bytes).hexdigest()):
                    annotation_errors.append("Retrospective annotation source hash does not match the original record.")
                annotation_status = "invalid source binding" if annotation_errors else "verified source binding"
                recorded_cc = result.get("cc_settings")
                recorded_guess = recorded_cc.get("scf_initial_guess") if isinstance(recorded_cc, dict) else None
                if not suffix and not recorded_guess and not annotation_errors:
                    inferred = annotation.get("inferred_cc_initial_guess")
                    if isinstance(inferred, str):
                        label = f"Historical CC reference · inferred {inferred}"
            elif annotation.get("status") != "missing":
                annotation_status = "invalid source binding"
                annotation_errors.append(annotation["load_error"])
            self.catalog["references"].append({**result, "id": identifier, "label": label, "source_href": href,
                "source": str(base / "method_comparison.json"),
                "retrospective_annotation": annotation if not annotation.get("load_error") else None,
                "retrospective_annotation_href": annotation_href, "retrospective_annotation_status": annotation_status,
                "annotation_errors": annotation_errors})

    def build_campaign(self, root, index):
        identifier = f"campaign-{index}"
        plan, plan_href = self.record(root, "plan.json", f"{identifier}-plan")
        ledger, ledger_href = self.record(root, "campaign.json", f"{identifier}-ledger")
        campaign = {"id": identifier, "label": root.label, "status": "loaded", "entries": [],
                    "source_href": plan_href, "ledger_href": ledger_href,
                    "interpretation": plan.get("interpretation"), "design_validated": False,
                    "automatic_ranking": False, "notices": []}
        self.catalog["campaigns"].append(campaign)
        if plan.get("load_error"):
            campaign.update(status=plan["status"], load_error=plan["load_error"])
            return
        if plan.get("schema_version") != 1:
            campaign.update(status="invalid evidence", load_error="Unsupported campaign plan schema; expected schema_version=1")
            return
        entries = plan.get("designs")
        if not isinstance(entries, list) or len(entries) > MAX_POSES:
            campaign.update(status="unreadable", load_error="Campaign must contain at most 100 design entries")
            return
        expected_plan_hash = ledger.get("plan_sha256")
        actual_hash = hashlib.sha256(self.raw[f"{identifier}-plan"][0]).hexdigest()
        campaign["plan_hash_matches_ledger"] = expected_plan_hash == actual_hash
        campaign_errors = []
        if ledger.get("schema_version") != 1:
            campaign_errors.append("Missing or unsupported campaign ledger schema.")
        if not campaign["plan_hash_matches_ledger"]:
            campaign_errors.append("Plan hash does not match the ledger, or the ledger is unavailable; statuses are unverified.")
        controls, quantum = plan.get("controls"), plan.get("quantum_settings")
        if not isinstance(controls, dict) or not all(key in controls for key in ("stage", "state", "fmax", "steps", "images")):
            campaign_errors.append("Campaign calculation controls are missing or incomplete.")
            controls = {}
        elif (controls["stage"] not in ("singlepoint", "relax", "path") or controls["state"] not in ("initial", "final") or
                type(controls["fmax"]) not in (int, float) or controls["fmax"] <= 0 or
                type(controls["steps"]) is not int or controls["steps"] < 1 or
                type(controls["images"]) is not int or controls["images"] < 5):
            campaign_errors.append("Campaign calculation controls contain unsupported values.")
        if not isinstance(quantum, dict) or not quantum:
            campaign_errors.append("Campaign quantum settings are missing.")
            quantum = {}
        campaign_errors.extend(quantum_settings_errors(quantum))
        campaign["notices"].extend(campaign_errors)
        campaign["integrity_errors"] = campaign_errors
        attempts_by_id = ledger.get("attempts", {})
        if not isinstance(attempts_by_id, dict):
            attempts_by_id = {}
        for number, entry in enumerate(entries):
            row = {"id": f"{identifier}-pose-{number}", "label": f"Pose {number + 1}", "status": "unreadable",
                   "pose": {}, "initial_structure_id": None, "final_structure_id": None,
                   "integrity_errors": list(campaign_errors), "status_verified": False}
            campaign["entries"].append(row)
            try:
                if not isinstance(entry, dict) or not isinstance(entry.get("id"), str):
                    raise EvidenceError("Malformed campaign design entry")
                row["label"], row["source_id"] = entry["id"], entry["id"]
                row["pose"] = entry.get("pose", {})
                relative = entry.get("design")
                SafeRoot.parts(relative)
                design, design_href = self.record(root, relative, f"{row['id']}-design")
                if design.get("load_error"):
                    raise EvidenceError(design["load_error"])
                row["source_href"] = design_href
                if design.get("schema_version") != 1 or design.get("length_unit") != "angstrom":
                    raise EvidenceError("Imported designs require schema_version=1 and explicit length_unit='angstrom'; units are never inferred")
                design_quantum = design.get("quantum")
                row["integrity_errors"].extend(quantum_settings_errors(design_quantum))
                if not isinstance(design_quantum, dict) or any(quantum.get(key) != value for key, value in design_quantum.items()):
                    row["integrity_errors"].append("Design quantum settings disagree with the campaign plan or are missing.")
                actual_hashes = {"design_sha256": hashlib.sha256(self.raw[f"{row['id']}-design"][0]).hexdigest()}
                for state in ("initial", "final"):
                    endpoint = design.get(state)
                    SafeRoot.parts(endpoint)
                    relative_endpoint = PurePosixPath(relative).parent / endpoint
                    row[f"{state}_structure_id"] = self.structure(f"{row['id']}-{state}", root,
                        relative_endpoint, design, f"{row['label']} · {state}", "pose", state=state)
                    endpoint_id = row[f"{state}_structure_id"]
                    actual_hashes[f"{state}_sha256"] = self.structures[endpoint_id]["source_sha256"] if endpoint_id else None
                expected_hashes = entry.get("input_hashes")
                row["input_hashes_match"] = expected_hashes == actual_hashes
                if not row["input_hashes_match"]:
                    row["integrity_errors"].append("Design or coordinate bytes differ from the campaign's recorded input hashes.")
                initial = self.structures.get(row["initial_structure_id"])
                final = self.structures.get(row["final_structure_id"])
                if not initial or not final:
                    row["integrity_errors"].append("Initial or final geometry is unavailable.")
                elif initial["symbols"] != final["symbols"] or any(
                    math.dist(initial["positions"][index], final["positions"][index]) > 1e-8 for index in initial["fixed_indices"]
                ):
                    row["integrity_errors"].append("Endpoints disagree on atom order or fixed-anchor coordinates.")
                attempts = attempts_by_id.get(entry["id"])
                if not isinstance(attempts, list):
                    row["integrity_errors"].append("Attempt ledger is unavailable for this pose.")
                    continue
                row["attempt_count"] = len(attempts)
                if not attempts:
                    row["status"] = row["recorded_status"] = "pending"
                    continue
                latest = attempts[-1]
                if not isinstance(latest, dict):
                    raise EvidenceError("Malformed campaign attempt")
                row["status"] = row["recorded_status"] = latest.get("status", "unknown")
                output = latest.get("output")
                SafeRoot.parts(output)
                result, result_href = self.record(root, PurePosixPath(output) / "result.json", f"{row['id']}-result")
                row["recorded_result"], row["result_href"] = result, result_href
                row["calculation_status"] = result.get("status", "unknown")
                if result.get("load_error"):
                    row["evidence_error"] = result["load_error"]
                    if row["recorded_status"] == "completed" or latest.get("result_sha256"):
                        row["integrity_errors"].append("Recorded completed/result evidence is missing or unreadable.")
                else:
                    result_digest = hashlib.sha256(self.raw[f"{row['id']}-result"][0]).hexdigest()
                    row["result_hash_matches"] = latest.get("result_sha256") == result_digest if latest.get("result_sha256") else None
                    if row["result_hash_matches"] is False:
                        row["integrity_errors"].append("Recorded result bytes differ from the attempt digest.")
                    if row["recorded_status"] == "completed" and row["result_hash_matches"] is not True:
                        row["integrity_errors"].append("Completed attempt lacks a matching recorded result digest.")
                    if result.get("input_hashes") != expected_hashes:
                        row["integrity_errors"].append("Result input hashes disagree with the campaign design.")
                    if result.get("stage") != controls.get("stage") or result.get("state") != controls.get("state"):
                        row["integrity_errors"].append("Result stage/state differs from the plan.")
                    expected_optimization = {"fmax_ev_per_angstrom": controls.get("fmax"), "max_steps_per_stage": controls.get("steps"), "images": controls.get("images")}
                    if result.get("optimization") != expected_optimization:
                        row["integrity_errors"].append("Result optimization controls differ from the plan.")
                    if result.get("quantum_settings") != quantum:
                        row["integrity_errors"].append("Result quantum settings differ from the plan.")
                    if row["recorded_status"] == "completed" and result.get("status") != "completed":
                        row["integrity_errors"].append("Completed ledger entry disagrees with the result status.")
                    if result.get("status") == "completed" and initial:
                        row["integrity_errors"].extend(completed_result_errors(result, controls, initial["atom_count"], initial["fixed_indices"]))
            except (EvidenceError, OSError, TypeError, AttributeError) as exc:
                row["status"], row["load_error"] = "unreadable", str(exc)
                row["integrity_errors"].append(str(exc))
            finally:
                row["status_verified"] = not row["integrity_errors"]
                if row["integrity_errors"]:
                    row["status"] = "invalid evidence"
        if campaign_errors or any(row["integrity_errors"] for row in campaign["entries"]):
            campaign["status"] = "invalid evidence"

    def build(self):
        self.build_archives()
        for index, root in enumerate(self.campaign_roots):
            self.build_campaign(root, index)
        return self

    def scope_ids(self):
        """Make links snapshot-specific; never silently resolve an older ID anew."""
        token = uuid4().hex
        raw_ids = {identifier: f"{token}:{identifier}" for identifier in self.raw}
        structure_ids = {identifier: f"{token}:{identifier}" for identifier in self.structures}

        def links(record):
            for field in ("id", "default_structure_id", "structure_id", "initial_structure_id", "final_structure_id"):
                identifier = record.get(field)
                if isinstance(identifier, str):
                    if identifier in structure_ids:
                        record[field] = structure_ids[identifier]
                    elif field == "id" and identifier in raw_ids:
                        record[field] = raw_ids[identifier]
            for field, value in list(record.items()):
                if field.endswith("_href") and isinstance(value, str) and value.startswith("/api/evidence?"):
                    identifiers = parse_qs(urlsplit(value).query).get("id", [])
                    if len(identifiers) == 1 and identifiers[0] in raw_ids:
                        record[field] = f"/api/evidence?id={quote(raw_ids[identifiers[0]], safe=':')}"

        links(self.catalog)
        for section in ("structures", "calculations", "references"):
            for record in self.catalog[section]:
                links(record)
        links(self.catalog["direct_df_comparison"])
        for campaign in self.catalog["campaigns"]:
            links(campaign)
            for row in campaign["entries"]:
                links(row)
        for record in self.structures.values():
            links(record)
        self.raw = {raw_ids[key]: value for key, value in self.raw.items()}
        self.structures = {structure_ids[key]: value for key, value in self.structures.items()}
        self.catalog["snapshot_id"] = token
        return self


class WorkbenchStore:
    """Publish one consistent read-only snapshot per catalog refresh."""

    def __init__(self, project_root=PROJECT_ROOT, campaign_root=None, campaigns=None):
        self.project = SafeRoot(Path(project_root), "project")
        boundary = SafeRoot(Path(campaign_root or project_root), "campaigns")
        selected = campaigns if campaigns is not None else ["examples/pose-campaign"]
        self.campaigns = []
        self.configuration_notices = []
        for relative in selected:
            try:
                self.campaigns.append(boundary.directory(relative))
            except FileNotFoundError:
                self.configuration_notices.append(f"Campaign directory unavailable: {relative}")
        self._lock = threading.RLock()
        self.snapshot = None
        self.catalog()

    def catalog(self):
        snapshot = CatalogBuilder(self.project, self.campaigns).build().scope_ids()
        snapshot.catalog["notices"].extend(self.configuration_notices)
        with self._lock:
            self.snapshot = snapshot
        return snapshot.catalog

    def structure(self, identifier):
        with self._lock:
            if identifier not in self.snapshot.structures:
                raise KeyError("Unknown or unavailable structure ID")
            return self.snapshot.structures[identifier]

    def evidence(self, identifier):
        with self._lock:
            if identifier not in self.snapshot.raw:
                raise KeyError("Unknown or unavailable evidence ID")
            return self.snapshot.raw[identifier]


def make_server(store, port=8765):
    class Handler(BaseHTTPRequestHandler):
        server_version = "MolecularWorkbench/1.0"

        def send_content(self, status, body, content_type="application/json; charset=utf-8", head=False):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; object-src 'none'; frame-ancestors 'none'; base-uri 'none'")
            self.end_headers()
            if not head:
                self.wfile.write(body)

        def error(self, status, message, head=False):
            self.send_content(status, json.dumps({"error": message}).encode(), head=head)

        def do_HEAD(self):
            self.do_GET(head=True)

        def do_GET(self, head=False):
            host = self.headers.get("Host", "")
            actual_port = self.server.server_port
            if host not in (f"127.0.0.1:{actual_port}", f"localhost:{actual_port}"):
                self.error(403, "Only the local workbench host is allowed", head)
                return
            parsed = urlsplit(self.path)
            query = parse_qs(parsed.query)
            try:
                if parsed.path == "/api/catalog":
                    body = json.dumps(store.catalog(), allow_nan=False).encode()
                    self.send_content(200, body, head=head)
                elif parsed.path in ("/api/structure", "/api/evidence"):
                    identifiers = query.get("id", [])
                    if len(identifiers) != 1 or len(identifiers[0]) > 160:
                        raise EvidenceError("One registered artifact ID is required")
                    if parsed.path == "/api/structure":
                        self.send_content(200, json.dumps(store.structure(identifiers[0]), allow_nan=False).encode(), head=head)
                    else:
                        content, mime = store.evidence(identifiers[0])
                        self.send_content(200, content, mime, head=head)
                else:
                    static = {"/": ("index.html", "text/html; charset=utf-8"),
                              "/index.html": ("index.html", "text/html; charset=utf-8"),
                              "/app.css": ("app.css", "text/css; charset=utf-8"),
                              "/app.js": ("app.js", "text/javascript; charset=utf-8")}
                    if parsed.path not in static:
                        raise KeyError("Not found")
                    filename, mime = static[parsed.path]
                    self.send_content(200, SafeRoot(STATIC_ROOT).read_bytes(filename), mime, head=head)
            except KeyError:
                self.error(404, "Unknown, unavailable, or expired artifact; refresh evidence to load current links", head)
            except (EvidenceError, OSError, ValueError) as exc:
                self.error(400, str(exc), head)

        def read_only(self):
            self.error(405, "This workbench is read-only; only GET and HEAD are supported")

        do_POST = do_PUT = do_PATCH = do_DELETE = do_OPTIONS = read_only

    return ThreadingHTTPServer(("127.0.0.1", port), Handler)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765, help="Loopback port (default: 8765)")
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT, help="Repository evidence boundary")
    parser.add_argument("--campaign-root", type=Path, help="Allowed parent directory for campaign selection; defaults to project root")
    parser.add_argument("--campaign", action="append", help="Relative campaign directory inside campaign root; repeat to load several")
    args = parser.parse_args(argv)
    if not 0 <= args.port <= 65535:
        parser.error("port must be between 0 and 65535")
    try:
        store = WorkbenchStore(args.project_root, args.campaign_root, args.campaign)
    except (EvidenceError, OSError) as exc:
        parser.error(str(exc))
    server = make_server(store, args.port)
    print(f"Molecular workbench: http://127.0.0.1:{server.server_port}", flush=True)
    print("Read-only local evidence. No quantum calculations are launched. Press Ctrl-C to stop.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
