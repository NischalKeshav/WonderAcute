# proofs/

The theory track of the research plan (`project.md` Step 1, plus the closed forms behind the Step 3–5 data analysis).

| File | What it is |
|---|---|
| `formal-model.md` | Step 1. Address table and flag rules, Lemmas 1–3, Theorems 1–3 (Clifford circuits, Clifford+T, the exact class of free classical functions), the boundary table for non-Clifford gates, caveats for later steps, sources. |
| `propagation-tables.json` | 54 signed rows *G P G*† = ±*Q* for H, S, S†, CNOT, CZ, SWAP. Read by `tests/test_frame.py`. |
| `predictions.md` | Closed forms for gates saved per circuit family (Step 3), lookup costs and what the pass does to them (Step 4), the shared-oracle saving and break-even (Step 5), a noise bound (Step 6), and an audit of the working paper's own numbers. |
| `check_proofs.py` | Checks every identity in `formal-model.md` against the matrices; `--write` regenerates the JSON. |
| `check_logic.py` | Checks the reasoning: the invariant after every instruction of random circuits in every branch, Lemma 2 for general π, the compile-time form, Theorem 3 by enumeration, every boundary row, the QFT claim. |
| `check_predictions.py` | Recomputes every formula in `predictions.md` and compares it row by row with the benchmark CSVs in `results/` (found automatically in the repo or the code worktree). |
| `sources/` | arXiv copies of the papers checked in `formal-model.md` §9. |

Run everything:

```
python3 proofs/check_proofs.py && python3 proofs/check_logic.py && python3 proofs/check_predictions.py
```

Only NumPy is needed. Still unverified against the original text: Riesebos 2017, Nielsen & Chuang section numbers, Babbush 2018, Low–Kliuchnikov–Schaeffer 2018, Zalka 1999 (all marked in the documents).
