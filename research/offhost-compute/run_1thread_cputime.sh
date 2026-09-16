#!/bin/sh
# Single-thread CPU-time measurement of the 53-atom DF single point.
# /usr/bin/time is absent in this container, so a Python wrapper records
# wall, user CPU, system CPU and peak RSS around the exact CLI call.
cd "$(dirname "$0")"
echo "=== 1thread-df start $(date -u +%FT%TZ) loadavg: $(cat /proc/loadavg)" >> cputime-log.txt
python3 measure_cputime.py >> cputime-log.txt 2>&1 || true
echo "=== 1thread-df end   $(date -u +%FT%TZ) loadavg: $(cat /proc/loadavg)" >> cputime-log.txt
echo "CPUTIME DONE" >> cputime-log.txt
