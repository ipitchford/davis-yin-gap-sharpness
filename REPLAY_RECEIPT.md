# Replay record

On 20 October 2026, `python3 -B code/replay.py` passed. The wrapper launched the
exact validator in ordinary and optimized modes, required reported optimization
flags 0 and 1, and compared receipts after removing only that flag.

The validator checks 720 family parameter combinations—600 open-range and 120
zero-smoothness cases—3,960 splitting iterations, 7,920 proximal evaluations,
300 exact certificate cases, six one-step equality cases, 36 endpoint-family
cases and 45 endpoint-identity cases. Deliberately corrupted trajectory and
constant controls were required to raise `VerificationError` in both modes.

The public receipts are under `results/`. These exact finite checks support the
formulas and implementation. They do not prove arbitrary convex interpolation or
replace the analytic theorem.
