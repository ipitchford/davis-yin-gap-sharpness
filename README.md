# Exact finite-horizon gap constants and non-attainment for Davis–Yin splitting

Anonymous. Version 1.0.0-candidate, revised 20 October 2026.

Version DOI: [10.5281/zenodo.23505130](https://doi.org/10.5281/zenodo.23505130).

This paper determines exact worst-case coefficients for an averaged primal–dual
comparison gap of Davis–Yin three-operator splitting throughout `0 < alpha L < 2`.
It proves the sharp denominator `K + 1 - alpha L`, supplies a two-dimensional
family approaching both natural metric bounds, and classifies attainment: one
step can attain equality, whereas no positive-denominator finite-dimensional
instance does so after two or more steps.

The revision also proves the critical endpoint `alpha L = 2` for `K >= 2`, with
exact sharp suprema and a restricted non-attainment theorem for strictly convex
`f`. General endpoint attainment remains open. Smooth strongly convex
approximations preserve the suprema without a common conditioning lower bound.

## Start here

- `paper.pdf` and `paper.tex`: complete definitions, proofs and limits.
- `CLAIMS.md` / `CLAIMS.json`: claim-to-evidence map.
- `code/replay.py`: exact ordinary/optimized replay and negative controls.
- `results/`: regenerated receipts.
- `AI_INDEX.md`: agent navigation and non-claims.

## Replay

```sh
python3 -B code/replay.py
```

The wrapper verifies actual optimization flags, compares semantic receipts, and
requires deliberate trajectory and constant corruptions to raise
`VerificationError` in ordinary and optimized modes. The finite checks support
the formulas; the analytic argument establishes the universal result.

Text and generated data: CC0-1.0. Original code: MIT.
