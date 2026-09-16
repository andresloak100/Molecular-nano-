"""Run the 53-atom DF single point single-threaded, recording CPU time.

/usr/bin/time is unavailable here, so this wraps the same CLI entry point and
reports wall clock, user+system CPU (children included via os.times), peak RSS,
and page faults from resource.getrusage. This gives A2 a transferable CPU-time
figure that does not depend on the machine being otherwise idle.
"""
import json
import os
import resource
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs/53atom-df-1t"

t = os.times()
wall0 = time.perf_counter()
u0, s0 = t.user + t.children_user, t.system + t.children_system

from nanodesign.cli import main  # noqa: E402

rc = main([
    "calculate", str(HERE / "design-53atom-df-1t.json"),
    "--out", str(OUT), "--stage", "singlepoint", "--state", "initial",
])

wall1 = time.perf_counter()
t = os.times()
u1, s1 = t.user + t.children_user, t.system + t.children_system
ru = resource.getrusage(resource.RUSAGE_SELF)
ruc = resource.getrusage(resource.RUSAGE_CHILDREN)

summary = {
    "cli_return_code": rc,
    "loadavg": os.getloadavg(),
    "wall_seconds": round(wall1 - wall0, 3),
    "user_cpu_seconds": round(u1 - u0, 3),
    "system_cpu_seconds": round(s1 - s0, 3),
    "total_cpu_seconds": round((u1 - u0) + (s1 - s0), 3),
    "peak_rss_mb": round(max(ru.ru_maxrss, ruc.ru_maxrss) / 1024, 1),
    "involuntary_context_switches": ru.ru_nivcsw + ruc.ru_nivcsw,
    "major_page_faults": ru.ru_majflt + ruc.ru_majflt,
    "requested_threads": 1,
    "note": (
        "CPU time is machine-transferable in a way wall clock is not: it is "
        "insensitive to whether the host was otherwise busy. cpu/wall ~1 "
        "confirms the run was compute-bound on one core and uncontended."
    ),
}
Path(HERE / "cputime-summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
