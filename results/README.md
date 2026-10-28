# Results

Produced by the scripts in `benchmarks/`; regenerate with `uv run python benchmarks/<name>.py`.
Every run records its inputs (seeds, gate set, basis, coupling map, commit) next to its output.
Code: `src/qarch/` (tests: `uv run pytest`, 259 tests).

## Step 2 — compiler pass (`src/qarch/passes.py`)

`compress()` absorbs X, Y, Z and SWAP into an address table (a qubit permutation plus a Pauli
frame), propagates the frame through Clifford gates, rewrites T/Tdg/RZ/P under a pending X,
flushes only the X part before diagonal non-Clifford gates, and flushes fully before anything
else. Verified against the proofs agent's propagation tables, exact operators (n ≤ 10), Stim
tableaux (n = 64) and an exact Stim lock-step measurement check. Not supported: control flow.

## Step 3 — gate savings (`gate_savings_{all,line}*`)

Columns: `raw`, `pass`, `o3` (Qiskit optimization_level=3), `o3+pass`, `pass+o3`.
Gates saved per qubit (slope of a linear fit, n = 4…64 step 4, 3 seeds), all-to-all coupling:

| family | pass vs raw | O3+pass vs O3 | pass+O3 vs O3 |
|---|---|---|---|
| random layered, depth 32 (`DEFAULT_LAYER_PROBS`) | 8.11 (R² 0.99) | 2.67 (0.95) | 3.33 (0.96) |
| Grover, Clifford+T MCZ, 1 iteration | 4.00 (1.00) | 2.00 (1.00) | 6.00 (1.00) |
| random Clifford (Qiskit synthesis) | 0.98 (0.98) | 1.02 (0.96) | 0.78 (0.96) |
| QFT | 0.50 (1.00) = ⌊n/2⌋ SWAPs | 0 | 0 |
| Bernstein–Vazirani | constant 1 (the ancilla X) | 1 | 1 |
| GHZ, repetition-code syndrome extraction | 0 | 0 | 0 |
| Shor-15 (order-aware compiled circuit, not general Shor) | 0.5 per counting qubit | 0 | 0 |

Readings.
- Pauli removal never saves a two-qubit gate; only SWAP removal does (`cx_equiv` column).
- Qiskit's O3 already removes SWAPs (`ElidePermutations`), so the QFT and Shor SWAP savings
  overlap with the baseline. The whitepaper's "1 × n" QFT figure is ⌊n/2⌋ and already standard.
- What the pass adds on top of O3 is Pauli absorption through Cliffords and the T rewrite:
  about 2 gates per qubit on Grover, 1–3 per qubit on random Clifford circuits, depending on
  the Pauli density of the input (stated in `*_params.json`).
- On the line coupling map, `o3+pass` deletes routing SWAPs by relabeling and therefore no
  longer obeys the hardware graph: it is the "free relabeling" assumption of the architecture,
  not a hardware saving. `pass+o3` is the routed, hardware-fair column and tracks the
  all-to-all numbers.
- Grover with a native multi-controlled Z is reported for n ≤ 10 only and is noisy on top of
  O3 (O3 resynthesizes the MCZ differently at each n). Its qubit counts are not comparable to
  the Clifford+T variant, which uses n − 3 clean ancillas.

## Step 4 — oracle lookup cost (`oracle_cost.*`)

Per lookup, d = 1 data bit, built circuits (closed forms in `qarch.qram` agree at every N):

| N | bucket brigade qubits / CSWAP / depth | table lookup qubits / Toffoli / depth |
|---|---|---|
| 8 | 35 / 72 / 65 | 6 / 12 / 23 |
| 64 | 262 / 732 / 125 | 12 / 124 / 223 |
| 1024 | 4106 / 12236 / 205 | 20 / 2044 / 3586 |

The table lookup is the plain compute/uncompute unary iteration (2(N−2) Toffolis); Babbush et
al.'s N−1 figure is included in the CSV. The select-swap row is a closed form written from
memory and not yet checked against Low–Kliuchnikov–Schaeffer. The lookup's own X gates are
mostly not absorbable: in the table lookup every absorbed X is flushed again before the next
Toffoli; in the bucket brigade about 28 % survive. The pass does not make lookups cheaper.

## Step 5 — shared oracle (`shared_oracle.*`)

k Grover instances, one queued oracle vs k copies, N = 1024, k = 3:
bucket brigade 12,318 → 4,126 qubits (66.5 % fewer; the whitepaper's 66 % holds for qubit
count) but depth × 2.81, so qubit·depth × 0.94; table lookup 33 % fewer qubits, depth × 2.99,
qubit·depth × 1.99. Sharing pays only when the oracle dominates the qubit count, and it always
costs depth because the memory answers one query at a time.

## Step 6 — noise (`noise.*`)

Amplitude damping + dephasing after every gate (p1 = 1e-3, p2 = 1e-2) and on every idle qubit
per tick (2e-3), 1000 shots per seed. Quick run (3 seeds), unbounded ops per tick:
Grover-3 success 0.634 → 0.690, Grover-4 0.127 → 0.195 (compressed higher, CIs disjoint),
BV-6 0.959 → 0.962 (no difference). Packing more operations into a tick always helps in this
model because it has no crosstalk term; a tick-count optimum needs a crosstalk model, which
is future work.
