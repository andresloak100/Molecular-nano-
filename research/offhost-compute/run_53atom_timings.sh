#!/bin/sh
# Serial, uncontended 53-atom single-point timings (DF then direct), 4 threads.
set -e
cd "$(dirname "$0")"
for name in df direct; do
  echo "=== $name start $(date -u +%FT%TZ) loadavg: $(cat /proc/loadavg)" >> timing-log.txt
  python3 -m nanodesign calculate "design-53atom-${name}-4t.json" \
    --out "runs/53atom-${name}-4t" --stage singlepoint --state initial \
    >> timing-log.txt 2>&1
  echo "=== $name end   $(date -u +%FT%TZ) loadavg: $(cat /proc/loadavg)" >> timing-log.txt
done
echo "ALL DONE" >> timing-log.txt
