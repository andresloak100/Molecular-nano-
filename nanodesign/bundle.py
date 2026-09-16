"""Bounded, byte-preserving directory exports; no scientific interpretation.

Only Python's standard library is used. A published manifest is a transport
integrity record, not authentication, reference closure or chemical validation.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat


FORMAT = "nanodesign-evidence-bundle"
FORMAT_VERSION = 1
MAX_MANIFEST_BYTES = 16 * 1024**2


class BundleError(ValueError):
    """Unsafe, changed, incomplete or over-budget evidence cannot be exported."""


@dataclass(frozen=True)
class BundleLimits:
    max_files: int = 10_000
    max_file_bytes: int = 32 * 1024**2
    max_total_bytes: int = 256 * 1024**2
    max_entries: int = 20_000
    max_depth: int = 32

    def __post_init__(self):
        for name, value in asdict(self).items():
            if type(value) is not int or value <= 0:
                raise BundleError(f"{name} must be a positive integer")


def _limits(value):
    if value is None:
        return BundleLimits()
    if not isinstance(value, BundleLimits):
        raise BundleError("limits must be a BundleLimits instance")
    return value


def _relative(value, limits):
    if not isinstance(value, str) or not value or value.startswith("/") or "\\" in value or "\x00" in value:
        raise BundleError("Artifact path must be a nonempty relative POSIX path")
    parts = value.split("/")
    if any(p in {"", ".", ".."} for p in parts) or ":" in parts[0]:
        raise BundleError(f"Unsafe artifact path: {value!r}")
    if len(parts) > limits.max_depth or len(value.encode("utf-8")) > 4096:
        raise BundleError("Artifact path exceeds depth or length limit")
    return parts


def _location(value):
    path = Path(value).absolute()
    if ".." in path.parts:
        raise BundleError("Explicit directory paths must not contain '..'")
    return path


def _identity(info):
    return info.st_dev, info.st_ino, stat.S_IFMT(info.st_mode)


def _signature(info):
    return (*_identity(info), info.st_size, info.st_mtime_ns, info.st_ctime_ns)


class _Root:
    """An open directory descriptor, pinned to its original selected location."""

    def __init__(self, path):
        self.path = _location(path)
        if self.path.is_symlink():
            raise BundleError("Selected directory cannot be a symlink")
        self.resolved = self.path.resolve(strict=True)
        info = self.path.lstat()
        if not stat.S_ISDIR(info.st_mode):
            raise BundleError("Selected path must be a directory")
        if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
            raise BundleError("Safe directory export requires POSIX no-follow directory access")
        self.fd = os.open(self.resolved, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        self.identity = _identity(os.fstat(self.fd))
        if self.identity != _identity(info):
            self.close()
            raise BundleError("Selected directory changed while opening")
        try:
            self.check()
        except BaseException:
            self.close()
            raise

    def check(self):
        info = self.path.lstat()
        if stat.S_ISLNK(info.st_mode) or _identity(info) != self.identity or self.path.resolve(strict=True) != self.resolved:
            raise BundleError("Selected directory or ancestor was replaced")

    def close(self):
        if self.fd is not None:
            os.close(self.fd)
            self.fd = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


@contextmanager
def _parent(root_fd, parts):
    fd = os.dup(root_fd)
    try:
        for part in parts[:-1]:
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = child
        yield fd, parts[-1]
    finally:
        os.close(fd)


def _inventory(root, limits):
    files, directories = {}, {}
    total = 0

    def walk(fd, prefix):
        nonlocal total
        # Count before sorting: even an adversarial directory listing is bounded.
        names = []
        with os.scandir(fd) as entries:
            for entry in entries:
                names.append(entry.name)
                if len(names) + len(files) + len(directories) > limits.max_entries:
                    raise BundleError("Directory exceeds max_entries")
        for name in sorted(names):
            relative = f"{prefix}/{name}" if prefix else name
            _relative(relative, limits)
            info = os.stat(name, dir_fd=fd, follow_symlinks=False)
            if stat.S_ISDIR(info.st_mode):
                directories[relative] = _signature(info)
                child = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
                try:
                    if _signature(os.fstat(child)) != _signature(info):
                        raise BundleError("Directory changed during enumeration")
                    walk(child, relative)
                finally:
                    os.close(child)
                if _signature(os.stat(name, dir_fd=fd, follow_symlinks=False)) != _signature(info):
                    raise BundleError("Directory changed during enumeration")
            elif stat.S_ISREG(info.st_mode):
                files[relative] = _signature(info)
                if len(files) > limits.max_files or info.st_size > limits.max_file_bytes:
                    raise BundleError("Evidence exceeds max_files or max_file_bytes")
                total += info.st_size
                if total > limits.max_total_bytes:
                    raise BundleError("Evidence exceeds max_total_bytes")
            else:
                raise BundleError(f"Symlink or nonregular entry is not allowed: {relative}")
            if len(files) + len(directories) > limits.max_entries:
                raise BundleError("Directory exceeds max_entries")

    root.check()
    walk(root.fd, "")
    root.check()
    return files, directories


def _read(root, relative, expected, limit, *, depth=1024):
    parts = _relative(relative, BundleLimits(max_depth=depth))
    root.check()
    with _parent(root.fd, parts) as (parent, name):
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=parent)
        try:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode) or before.st_size > limit:
                raise BundleError(f"Not a bounded regular file: {relative}")
            if expected is not None and _signature(before) != expected:
                raise BundleError(f"Source changed before capture: {relative}")
            with os.fdopen(os.dup(fd), "rb") as stream:
                raw = stream.read(limit + 1)
            if len(raw) > limit or len(raw) != before.st_size or _signature(os.fstat(fd)) != _signature(before):
                raise BundleError(f"File changed during capture: {relative}")
            if _signature(os.stat(name, dir_fd=parent, follow_symlinks=False)) != _signature(before):
                raise BundleError(f"File was replaced during capture: {relative}")
        finally:
            os.close(fd)
    root.check()
    return raw


def _write(root, relative, raw):
    with _parent(root.fd, relative.split("/")) as (parent, name):
        fd = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=parent)
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())


def _mkdir(root, relative):
    with _parent(root.fd, relative.split("/")) as (parent, name):
        os.mkdir(name, 0o700, dir_fd=parent)


def _new_output(destination, source_root):
    # Pin the parent before exclusive creation: a replaced parent must not
    # redirect even an incomplete output into the read-only source tree.
    with _Root(destination.parent) as parent:
        if (parent.resolved / destination.name).is_relative_to(source_root.resolved):
            raise BundleError("Destination must be outside the source tree")
        parent.check()
        os.mkdir(destination.name, 0o700, dir_fd=parent.fd)
        parent.check()
        return _Root(destination)


def _outer_limits(limits):
    return BundleLimits(max_files=limits.max_files + 1,
        max_file_bytes=max(limits.max_file_bytes, MAX_MANIFEST_BYTES),
        max_total_bytes=limits.max_total_bytes + MAX_MANIFEST_BYTES,
        max_entries=limits.max_entries + 2, max_depth=limits.max_depth + 1)


def _check_written_payload(output, records, directories, limits):
    files, dirs = _inventory(output, _outer_limits(limits))
    expected_files = {f"payload/{record['path']}" for record in records}
    expected_dirs = {"payload"} | {f"payload/{relative}" for relative in directories}
    if set(files) != expected_files or set(dirs) != expected_dirs:
        raise BundleError("Output inventory changed before manifest publication")
    for record in records:
        member = f"payload/{record['path']}"
        raw = _read(output, member, files[member], limits.max_file_bytes)
        if len(raw) != record["size_bytes"] or hashlib.sha256(raw).hexdigest() != record["sha256"]:
            raise BundleError("Output bytes changed before manifest publication")
    if _inventory(output, _outer_limits(limits)) != (files, dirs):
        raise BundleError("Output changed during final recheck")


def export_bundle(source, destination, *, limits=None):
    """Copy one explicit tree to a new destination; publish manifest last.

    A failure before manifest publication may leave partial payload bytes but
    no manifest.json. Interruption during publication can leave an extra pending
    file, which verification rejects. No cleanup or retry overwrites evidence.
    Parent directories must already exist.
    """
    limits = _limits(limits)
    try:
        with _Root(source) as source_root:
            destination = _location(destination)
            resolved_destination = destination.parent.resolve(strict=True) / destination.name
            if resolved_destination.is_relative_to(source_root.resolved):
                raise BundleError("Destination must be outside the source tree")
            files, directories = _inventory(source_root, limits)
            # Exclusive creation also refuses existing empty directories/symlinks.
            with _new_output(destination, source_root) as output:
                _mkdir(output, "payload")
                for relative in sorted(directories, key=lambda p: (p.count("/"), p)):
                    _mkdir(output, f"payload/{relative}")
                records = []
                for relative, signature in sorted(files.items()):
                    raw = _read(source_root, relative, signature, limits.max_file_bytes)
                    _write(output, f"payload/{relative}", raw)
                    records.append({"path": relative, "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
                # Recheck complete inventory and each exact source byte stream.
                # This detects mutation; it is not a global filesystem snapshot.
                if _inventory(source_root, limits) != (files, directories):
                    raise BundleError("Source inventory changed during export")
                for record in records:
                    relative = record["path"]
                    raw = _read(source_root, relative, files[relative], limits.max_file_bytes)
                    if hashlib.sha256(raw).hexdigest() != record["sha256"]:
                        raise BundleError(f"Source bytes changed during export: {relative}")
                if _inventory(source_root, limits) != (files, directories):
                    raise BundleError("Source inventory changed during final recheck")
                _check_written_payload(output, records, directories, limits)
                manifest = {
                    "format": FORMAT, "format_version": FORMAT_VERSION,
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "source_name": source_root.resolved.name,
                    "files": records, "directories": sorted(directories),
                    "file_count": len(records), "total_bytes": sum(r["size_bytes"] for r in records),
                    "limits": asdict(limits), "scientific_validation": "not_assessed",
                    "reference_closure_verified": False,
                    "snapshot": {"atomic": False, "source_rechecked": True},
                    "interpretation": "Exact captured bytes with detectable source-change checks. No simultaneous cross-file snapshot, source authentication, reference closure, chemical accuracy or completed calculation is certified. Original records, omissions and failures are preserved without interpretation.",
                }
                raw_manifest = (json.dumps(manifest, indent=2, allow_nan=False) + "\n").encode()
                if len(raw_manifest) > MAX_MANIFEST_BYTES:
                    raise BundleError("Manifest exceeds size limit")
                output.check()
                # Publish complete bytes with a no-overwrite hard link. An
                # interrupted temporary write cannot become a final manifest.
                _write(output, ".manifest.pending", raw_manifest)
                output.check()
                os.link(".manifest.pending", "manifest.json", src_dir_fd=output.fd, dst_dir_fd=output.fd, follow_symlinks=False)
                os.unlink(".manifest.pending", dir_fd=output.fd)
                os.fsync(output.fd)
                return manifest
    except (OSError, UnicodeError) as error:
        raise BundleError(str(error)) from error


def _strict_json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise BundleError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value):
        raise BundleError(f"Non-finite JSON constant: {value}")

    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def _manifest_records(manifest, limits):
    if not isinstance(manifest, dict) or manifest.get("format") != FORMAT or type(manifest.get("format_version")) is not int or manifest["format_version"] != FORMAT_VERSION:
        raise BundleError("Unsupported bundle manifest format/version")
    if manifest.get("scientific_validation") != "not_assessed" or manifest.get("reference_closure_verified") is not False:
        raise BundleError("Manifest cannot certify science or external reference closure")
    if manifest.get("snapshot") != {"atomic": False, "source_rechecked": True}:
        raise BundleError("Unsupported snapshot claim")
    files, dirs = manifest.get("files"), manifest.get("directories")
    if not isinstance(files, list) or not isinstance(dirs, list) or len(files) > limits.max_files or len(files) + len(dirs) > limits.max_entries:
        raise BundleError("Invalid or over-budget manifest inventory")
    records, directories, total = {}, set(), 0
    for directory in dirs:
        _relative(directory, limits)
        if directory in directories:
            raise BundleError("Duplicate manifest directory")
        directories.add(directory)
    for record in files:
        if not isinstance(record, dict) or set(record) != {"path", "size_bytes", "sha256"}:
            raise BundleError("Invalid manifest file entry")
        relative = record["path"]
        _relative(relative, limits)
        if relative in records or relative in directories:
            raise BundleError("Duplicate or conflicting manifest path")
        size = record["size_bytes"]
        if type(size) is not int or size < 0 or size > limits.max_file_bytes:
            raise BundleError("Invalid or over-budget manifest file size")
        if not isinstance(record["sha256"], str) or re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is None:
            raise BundleError("Invalid SHA-256 in manifest")
        total += size
        records[relative] = record
    if total > limits.max_total_bytes:
        raise BundleError("Manifest exceeds max_total_bytes")
    for relative in set(records) | directories:
        parts = relative.split("/")
        if any("/".join(parts[:i]) not in directories for i in range(1, len(parts))):
            raise BundleError("Manifest omits a parent directory")
    if type(manifest.get("file_count")) is not int or manifest["file_count"] != len(records) or type(manifest.get("total_bytes")) is not int or manifest["total_bytes"] != total:
        raise BundleError("Manifest totals do not match its records")
    return records, directories


def verify_bundle(directory, *, limits=None):
    """Read-only, bounded integrity check against the supplied manifest.

    Invalid/incomplete bundles return integrity_verified=False with issues.
    A matching manifest is not authenticated, and no source JSON is interpreted.
    """
    limits = _limits(limits)
    report = {"integrity_verified": False, "scientific_validation": "not_assessed",
              "reference_closure_verified": False, "issues": [], "file_count": 0, "total_bytes": 0}

    def issue(code, path, message):
        report["issues"].append({"code": code, "path": path, "message": message})

    try:
        with _Root(directory) as root:
            raw_manifest = _read(root, "manifest.json", None, MAX_MANIFEST_BYTES)
            manifest = _strict_json(raw_manifest)
            records, directories = _manifest_records(manifest, limits)
            # Account separately for the manifest and payload container.
            outer_limits = _outer_limits(limits)
            files, dirs = _inventory(root, outer_limits)
            expected_files = {"manifest.json"} | {f"payload/{p}" for p in records}
            expected_dirs = {"payload"} | {f"payload/{p}" for p in directories}
            for relative in sorted(expected_files - files.keys()):
                issue("missing_file", relative, "Listed file is missing")
            for relative in sorted(files.keys() - expected_files):
                issue("extra_file", relative, "File is not listed in the manifest")
            for relative in sorted(expected_dirs - dirs.keys()):
                issue("missing_directory", relative, "Listed directory is missing")
            for relative in sorted(dirs.keys() - expected_dirs):
                issue("extra_directory", relative, "Directory is not listed in the manifest")
            for relative, record in records.items():
                member = f"payload/{relative}"
                if member not in files:
                    continue
                raw = _read(root, member, files[member], limits.max_file_bytes)
                if len(raw) != record["size_bytes"]:
                    issue("size_mismatch", member, "Length differs from the manifest")
                if hashlib.sha256(raw).hexdigest() != record["sha256"]:
                    issue("hash_mismatch", member, "SHA-256 differs from the manifest")
                report["file_count"] += 1
                report["total_bytes"] += len(raw)
                if report["total_bytes"] > limits.max_total_bytes:
                    raise BundleError("Actual payload exceeds max_total_bytes")
            if _inventory(root, outer_limits) != (files, dirs) or _read(root, "manifest.json", files["manifest.json"], MAX_MANIFEST_BYTES) != raw_manifest:
                raise BundleError("Bundle changed during verification")
            report["manifest_sha256"] = hashlib.sha256(raw_manifest).hexdigest()
            report["integrity_verified"] = not report["issues"]
    except (BundleError, OSError, UnicodeError, ValueError, RecursionError) as error:
        issue("invalid_bundle", "", str(error))
    return report
