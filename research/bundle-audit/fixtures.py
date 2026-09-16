"""Tiny synthetic file trees for B2; no quantum imports or scientific outputs."""
from __future__ import annotations

import hashlib
from pathlib import Path


def make_source(path: Path) -> dict[str, bytes]:
    """Return the literal bytes the exporter must preserve, including failures."""
    files = {
        "design.json": b'{"fixture":"B2 synthetic","length_unit":"angstrom","design_validated":false}\n',
        "initial.xyz": b"2\nB2 synthetic fixture, not a calculation\nH 0 0 0\nH 0.74 0 0\n",
        "final.xyz": b"2\nB2 synthetic fixture, not a calculation\nH 0 0 0\nH 1.48 0 0\n",
        "attempts/failed/result.json": b'{"status":"failed","error":"synthetic SCF failure","design_validated":false}\n',
        "attempts/pending/note.txt": b"No calculation has been run.\n",
        "empty.txt": b"",
        "notes-Å.txt": "Literal UTF-8 Å and α.\n".encode(),
        ".run-note": b"Hidden files must not be silently excluded.\n",
        "opaque.bin": bytes(range(256)),
    }
    path.mkdir(parents=True, exist_ok=False)
    for relative, content in files.items():
        target = path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    return files


def file_bytes(root: Path) -> dict[str, bytes]:
    """Snapshot a trusted fixture tree; never follow source symlinks."""
    result = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise AssertionError(f"Unexpected symlink in trusted fixture: {path.name}")
        if path.is_file():
            result[path.relative_to(root).as_posix()] = path.read_bytes()
    return result


def fingerprints(root: Path) -> dict[str, dict]:
    return {name: {"size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
            for name, content in file_bytes(root).items()}
