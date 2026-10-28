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
- CP/RZZ rewrite (`rewrite_cp`, **off by default**; `--cp-rewrite` turns it on): a pending X on a controlled-phase gate becomes
  CP(−θ) plus a P(θ) on the other qubit instead of a flushed X. On O3 output this *lowers* the
  total-gate saving (Clifford+T Grover: 2 → 1 per qubit; native-MCZ Grover: +3 → −27 per
  qubit, because every CP in the synthesized MCZ gains a P) but *raises* the virtual-Z saving
  (`gates_vz`, diagonal one-qubit gates free: 7 and +26 per qubit). It is the right choice on
  virtual-Z hardware and the wrong one when every gate costs; QFT and Shor are unaffected.
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
per tick (p_idle = 2e-3); 1000 shots × 10 seeds; 95 % CI = 1.96 × max(across-seed standard
error, binomial error of the pooled shots). Rate rationale: Google's 2024 surface-code
processor (arXiv 2408.13687) reports, from memory and still to be confirmed against the paper,
roughly 3–4 × 10⁻⁴ per one-qubit gate and 3 × 10⁻³ per CZ; the values here are about three times worse, a generic near-term setting chosen so that
3–5-qubit Grover instances land at informative success probabilities instead of near 0 or 1.
`noise_idle_sweep.*` varies p_idle from 0 to 8 × 10⁻³ so the conclusion does not rest on one
guess. Two pairs are compared with identical seeds and shots: raw vs pass, and Qiskit O3 vs
pass followed by O3 (the step-3 baseline). Success probability, original → compressed,
unbounded ops per tick (clock ticks in brackets):

| circuit | original | compressed |
|---|---|---|
| Grover n=3 (2 iterations) | 0.621 ± 0.011 (57) | 0.680 ± 0.012 (53) |
| Grover n=4 (3 iterations) | 0.131 ± 0.007 (193) | 0.197 ± 0.009 (187) |
| Grover n=5 (4 iterations) | 0.021 ± 0.004 (417) | 0.047 ± 0.004 (409) |
| Bernstein–Vazirani n=6 | 0.955 ± 0.004 (6) | 0.966 ± 0.004 (5) |
| Bernstein–Vazirani n=10 | 0.853 (10) | 0.916 (9) |

The compressed circuit succeeds more often in every case; the gap is similar at one operation
per tick (Grover-4: 0.080 → 0.136). Packing more operations into a tick always helps in this
model because it has no crosstalk term; the "best number of operations per tick" asked for in
the plan needs a crosstalk penalty, which is future work.

`shared_oracle_noise.*` repeats the step-5 comparison under this noise model at N = 8 with the
table-lookup oracle (the bucket brigade's 35 qubits are out of statevector range), k = 1…3
instances, sweeping the idle rate: in the shared design the queued instances take idle error
while they wait. See the CSV for the numbers.
