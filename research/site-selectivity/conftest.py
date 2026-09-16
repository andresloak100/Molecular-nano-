"""Make this lane's modules importable by its own tests without installing them."""

from __future__ import annotations

from pathlib import Path
import sys

LANE = Path(__file__).resolve().parent
if str(LANE) not in sys.path:
    sys.path.insert(0, str(LANE))
