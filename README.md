# Optimizations to Quantum Computing Architectures Utilizing Classical Structures

Working folder for the ISEF project by Nischal Keshav.

| File | What it is |
|---|---|
| `research-plan.md` | ISEF Research Plan / Project Summary (current version, Markdown) |
| `research-plan.docx` | Same plan as a two-page Word file for submission |
| `research-directions.md` | Seven candidate research projects with the papers and textbooks behind each |
| `whitepaper.html` | Whitepaper on the original architecture, with diagrams and interactive charts. Opens in any browser. |
| `source/` | The original working paper, "Modern Architecture for Quantum Computers" (PDF) |

| `proofs/` | Formal model and propagation tables (step 1) |
| `src/qarch/` | Address-table compiler pass, circuit generators, QRAM circuits, noise model (steps 2–6) |
| `benchmarks/` | Scripts that produce `results/` (gate savings, oracle cost, shared oracle, noise) |
| `results/` | CSV/PNG outputs with a README summarising the findings |

## Running the code

```
uv sync                                        # Python 3.12 venv with Qiskit, Aer, Stim
uv run pytest                                  # rule tests, equivalence checks, QRAM checks
uv run python benchmarks/gate_savings.py --quick
uv run python benchmarks/oracle_cost.py
uv run python benchmarks/shared_oracle.py
uv run python benchmarks/noise.py --quick
```
