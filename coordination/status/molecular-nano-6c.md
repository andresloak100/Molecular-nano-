# molecular-nano-6c status (cloud remote-execution session)

Updated: 2026-09-16 ~22:40 UTC. (Prior updates ~21:50, ~22:30 UTC.)

**Alignment with A2's top outstanding item.** A2's REPORT.md names one
uncontended per-evaluation CPU-time measurement as the item its verdict hangs
on (8.6 days vs 20 hours per pose). The 53-atom DF + direct single points now
running here, serially on an idle 4-core container with loadavg logged, are
that measurement's off-host arm: they bound the true per-evaluation cost from
a host with zero contention and honored threading. Numbers land in
`research/offhost-compute/` the moment the runs complete. This serves root's
stated criterion for new support (a concrete requested gap, own paths only,
no duplication of an active lane's science).

**Also produced for G1:** `research/offhost-compute/gpu-readiness-offhost/
offhost-preflight.json` — G1's own inspector run on this container
(Linux/x86_64, PySCF 2.14.0, thread probe requested 2 / effective 2 /
honored TRUE, no CUDA stack). The Linux CPU complement to G1's Darwin
preflight; written into my lane, not G1's.

**Relayed a user directive** on workbench styling ("keep design
apple/chatgpt") to root/V1 with my interpretation flagged:
`coordination/messages/from-molecular-nano-6c-workbench-design-preference.md`.

Identity: Claude Code session in an isolated cloud Linux container (fresh
clone; does **not** share the Mac checkout or its load). Peer name
`molecular-nano-6c [d6a09f]`; not reachable by Mac-side socket messaging —
this repository is my channel, the user relays urgent items. Pushes go to
branch `claude/molecular-nanomachine-design-iu5ik5` (required); `origin/main`
is merged in to stay current.

## Lanes

1. **Org (user-directed, done):** `coordination/ORG.md` — lead assignment
   (Codex root), deputy (original scientific helper), operating rules.
   Revised ~22:20 UTC to defer all lane ownership to root's `ROSTER.md`;
   no dueling boards.
2. **Off-host compute (user-directed, in progress):**
   `research/offhost-compute/` — uncontended 4-thread timings of the exact
   archived 53-atom single points (DF then direct, serial, byte-identical
   input geometries, loadavg logged). Purpose: clean per-gradient cost for
   A2/C2 and a cross-platform reproduction check on the actual candidate.
   Requested roster inclusion from root; no core edits, no writes to
   archived records.
3. **Portability/verification (standing):** merged `origin/main` at
   `fc8103e` into my branch; will re-run the full suite on Linux after the
   quantum timings finish (not during — keeps the measurement uncontended).

## Verified earlier this session (Linux x86_64, 4 cores, 15 GiB RAM)

- PySCF 2.14.0 with working OpenMP; `lib.num_threads()` honors 4.
- 204/204 tests green at my branch base; def2-TZVP benchmark reproduced the
  Mac table to <0.01 kcal/mol (barrier −2.8906, reaction −28.045) in 32.7 s
  wall on 4 threads.
- `65f2dee` gitignore fix confirmed effective on a fresh clone.

## Requests

- **Root:** merge `claude/molecular-nanomachine-design-iu5ik5`; add the
  `offhost-compute` lane to `ROSTER.md` or object here; acknowledge (or
  object to) the user-directed lead assignment in your status file.
- **A2/C2:** off-host timing numbers land in `research/offhost-compute/`
  shortly — treat as different-hardware scenario data, not archive
  corrections.
- **G1:** this container volunteers as the off-host preflight target for
  your GPU/offload equivalence protocol.
- **Any lane:** bounded quantum jobs can be delegated here instead of adding
  to the Mac queue — spec in an addressed message, have the user ping me.

## Jobs running here

`run_53atom_timings.sh` (serial DF → direct 53-atom single points, 4
threads). Logged in `research/offhost-compute/timing-log.txt`. Nothing else.
Container is ephemeral: results are pushed as soon as they exist.
