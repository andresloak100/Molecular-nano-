"""Timing instrumentation that survives an oversubscribed host.

Every wall-clock number archived in this repository so far was taken while
several other single-threaded quantum jobs were running. At the time of
writing the one-minute load average on this eight-core host reached 209, so a
wall-clock measurement is a measurement of the scheduler, not of the
calculation.

Process CPU time does not have that problem. PySCF here has no OpenMP and runs
on one thread (see ``threads_honored``), so for a single-threaded process

    process CPU time ~= uncontended wall clock
    wall clock / CPU time = the contention factor actually suffered

Both are recorded. CPU time is the transferable cost of the calculation; wall
clock is what a person waiting on this host today experiences. The report must
not quote one where it means the other.

Caveat worth keeping: CPU time excludes time the process spent blocked on I/O
or swapped out, so it slightly understates true single-user wall clock. On a
memory-pressured host the gap widens. Treat CPU time as a lower bound on
uncontended wall clock, not an exact prediction.
"""

from __future__ import annotations

import os
import resource
import time
from contextlib import contextmanager
from typing import Any


def load_snapshot() -> dict[str, Any]:
    """Load averages and this host's core count at the moment of the call."""
    one, five, fifteen = os.getloadavg()
    cores = os.cpu_count() or 1
    return {
        "load_average_1min": one,
        "load_average_5min": five,
        "load_average_15min": fifteen,
        "logical_cores": cores,
        "oversubscription_1min": one / cores,
    }


def _usage() -> dict[str, float]:
    """Resource usage for this process and its finished children.

    User and system time are kept apart on purpose. Two sessions have argued
    that this host's load is memory-stall rather than CPU queueing. If that is
    right, the cost lands as *system* time (page-fault handling charged to the
    faulting process) and as major faults and involuntary context switches,
    not as user time. Recording the split turns that claim into something this
    lane can check in its own data instead of accepting or rejecting it.
    """
    me = resource.getrusage(resource.RUSAGE_SELF)
    kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "user": me.ru_utime + kids.ru_utime,
        "system": me.ru_stime + kids.ru_stime,
        "major_faults": float(me.ru_majflt + kids.ru_majflt),
        "minor_faults": float(me.ru_minflt + kids.ru_minflt),
        "voluntary_switches": float(me.ru_nvcsw + kids.ru_nvcsw),
        "involuntary_switches": float(me.ru_nivcsw + kids.ru_nivcsw),
    }


def _cpu_seconds() -> float:
    """User + system CPU seconds for this process and its finished children."""
    usage = _usage()
    return usage["user"] + usage["system"]


@contextmanager
def timed(record: dict[str, Any], label: str = "timing"):
    """Record wall clock, CPU time, contention factor and load into ``record``.

    ``record[label]`` is written on the way out, including when the body
    raises: a failed calculation still cost real time and that is a result.
    """
    start_wall = time.monotonic()
    start_usage = _usage()
    start_load = load_snapshot()
    failed = None
    try:
        yield
    except BaseException as exc:  # noqa: BLE001 - recorded, then re-raised
        failed = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        wall = time.monotonic() - start_wall
        end_usage = _usage()
        delta = {key: end_usage[key] - start_usage[key] for key in start_usage}
        cpu = delta["user"] + delta["system"]
        record[label] = {
            "wall_seconds": wall,
            "cpu_seconds": cpu,
            "user_seconds": delta["user"],
            "system_seconds": delta["system"],
            "system_fraction_of_cpu": (delta["system"] / cpu) if cpu > 0 else None,
            "major_page_faults": delta["major_faults"],
            "minor_page_faults": delta["minor_faults"],
            "involuntary_context_switches": delta["involuntary_switches"],
            "voluntary_context_switches": delta["voluntary_switches"],
            "contention_factor_wall_over_cpu": (wall / cpu) if cpu > 0 else None,
            "load_at_start": start_load,
            "load_at_end": load_snapshot(),
            "effective_threads_assumed": 1,
            "interpretation": (
                "user_seconds is the transferable single-thread cost of the "
                "arithmetic. system_seconds plus major_page_faults is the "
                "kernel overhead this process was charged, which is where "
                "memory-stall contention would appear if it is being charged "
                "here at all. wall_seconds is this host under the recorded "
                "load and is not a property of the calculation. None of these "
                "is a validated scientific result."
            ),
            "failed": failed,
        }
