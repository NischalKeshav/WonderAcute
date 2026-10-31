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
  all-to-all numbers on average, but routing is heuristic: on the line map `pass+o3` is worse
  than `o3` in CX-equivalents on 24 of 48 `random_layered` rows (worst +402), never on the
  structured families (`proofs/predictions.md` §6). Removing SWAPs before routing changes the
  router's starting point, and a different layout can cost more than the removed SWAPs saved.
- CP/RZZ rewrite (`rewrite_cp`, **off by default**; `--cp-rewrite` turns it on): a pending X on a controlled-phase gate becomes
  CP(−θ) plus a P(θ) on the other qubit instead of a flushed X. On O3 output this *lowers* the
  total-gate saving (Clifford+T Grover: 2 → 1 per qubit; native-MCZ Grover: +3 → −27 per
  qubit, because every CP in the synthesized MCZ gains a P) but *raises* the virtual-Z saving
  (`gates_vz`, diagonal one-qubit gates free: 7 and +26 per qubit). It is the right choice on
  virtual-Z hardware and the wrong one when every gate costs; QFT and Shor are unaffected.
- Exact Grover savings (derived in `proofs/predictions.md`, checked against the CSV): with z the
  number of zero bits in the marked item, one iteration saves exactly n + z gates for the
  native-MCZ variant and 2n + 2z for the Clifford+T variant vs raw. z depends on the seed, which
  is why those fits have R² < 1 while the per-qubit slope is exact.
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
Toffoli; in the bucket brigade about 28 % survive. Exactly (`proofs/predictions.md` §2a,
checked on every built row): flushed X = 4(N−1−k) + 4d(N−1) + 2dw with w the set bits of the
table, so gates removed = 6N + 2k − 2 − w at d = 1 and the flushed fraction tends to 2/3; the
only X pairs that cancel are a descent's trailing X against the ascent's leading X on the same
router. No CSWAP is ever removed, so Toffoli and T counts are unchanged. The pass does not make
lookups cheaper.

## Step 5 — shared oracle (`shared_oracle.*`)

k Grover instances, one queued oracle vs k copies, N = 1024, k = 3:
bucket brigade 12,318 → 4,126 qubits (66.5 % fewer; the whitepaper's 66 % holds for qubit
count) but depth × 2.81, so qubit·depth × 0.94; table lookup 33 % fewer qubits, depth × 2.99,
qubit·depth × 1.99. In closed form (`proofs/predictions.md`): with a address qubits per
instance, Q oracle qubits, lookup depth L and diffusion depth D, sharing saves a fraction
(1 − 1/k)·Q/(a + Q) of the qubits, so the whitepaper's 66 % is the k = 3 ceiling Q ≫ a, and
qubit·depth improves only while k < Q·D/(a·L). The whitepaper's own table (61/64/65/66 % for
its four cases) omits the address registers on the duplicated side; with them included the
values are 63.2/64.9/65.8/66.2 % (`proofs/predictions.md` §5). Sharing pays only when the oracle dominates the
qubit count, and it always costs depth because the memory answers one query at a time.

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

| circuit | raw | pass | O3 | pass + O3 |
|---|---|---|---|---|
| Grover n=3 (2 iterations) | 0.621 ± 0.011 (57) | 0.680 ± 0.012 (53) | 0.661 ± 0.011 (49) | 0.706 ± 0.012 (47) |
| Grover n=4 (3 iterations) | 0.131 ± 0.007 (193) | 0.197 ± 0.009 (187) | 0.153 ± 0.008 (188) | 0.215 ± 0.008 (190) |
| Grover n=5 (4 iterations) | 0.021 ± 0.004 (417) | 0.047 ± 0.004 (409) | 0.024 ± 0.003 (425) | 0.049 ± 0.005 (413) |
| Bernstein–Vazirani n=6 | 0.955 ± 0.004 (6) | 0.966 ± 0.004 (5) | 0.957 ± 0.004 (6) | 0.968 ± 0.004 (5) |
| Bernstein–Vazirani n=10 | 0.853 ± 0.007 (10) | 0.916 ± 0.005 (9) | 0.852 ± 0.007 (10) | 0.915 ± 0.005 (9) |

The compressed circuit succeeds more often in every pair, against the raw circuit and against
Qiskit's O3; the gap is similar at one operation per tick (Grover-4 raw: 0.080 → 0.136). The
idle sweep (`noise_idle_sweep.*`, O3 vs pass + O3, unbounded ops per tick) shows the advantage
growing with the idle rate: Grover-4 0.295 → 0.344 at p_idle = 0 and 0.046 → 0.087 at 8e-3;
BV-10 0.925 → 0.937 at 0 and 0.684 → 0.857 at 8e-3.

Where the BV-10 gain comes from (diagnostic at p_idle = 8e-3, O3 → pass+O3): re-inserting
the saved tick as an all-qubit identity gives 0.815, so the tick accounts for ~4 points; with
dephasing only the pair is 0.906 → 0.929, with amplitude damping only 0.745 → 0.915. Most of
the gap is therefore a **T1 asymmetry**: in the compressed circuit the answer lives in the
classical frame and the physical control qubits wait for readout in |0⟩, which amplitude
damping cannot touch, while in the O3 circuit the controls holding a 1 sit in |1⟩ and decay.
This is a real property of frame tracking under biased noise, but it is specific to amplitude
damping and to a model without readout error; `noise.py --channel dephasing|damping` separates
the two.

**Dephasing-only sweep (`noise_idle_sweep_dephasing.*`, O3 vs pass + O3).** With amplitude
damping switched off the compressed and baseline circuits are equal within CI for Grover at
every idle rate (n=4: 0.664 vs 0.664 at p_idle = 0, 0.341 vs 0.337 at 8e-3; n=5: 0.458 vs
0.460 and 0.068 vs 0.070; n=3: 0.860 vs 0.862 and 0.802 vs 0.803), and only BV keeps a small
tick-effect gain (BV-10: 0.963 → 0.972 at 0, 0.906 → 0.929 at 8e-3; BV-6: 0.974 → 0.978 at
8e-3). So the success-probability gains under the combined model are the amplitude-damping
asymmetry above, not the removed gates: the pass removes one-qubit gates at p1 = 1e-3, which is
negligible next to the untouched two-qubit gates at p2 = 1e-2. The reason is structural:
dephasing and depolarizing channels are invariant under Pauli conjugation, so moving a Pauli
into the frame cannot change what they do to the state; amplitude damping is not, and a frame
that keeps the physical register near |0…0⟩ (BV's controls, Grover's register during the
MCZ decomposition once the diffusion X layer is absorbed) is protected from it.
`noise_idle_sweep_depolarizing.*` (the Pauli-symmetric control, 1q/2q depolarizing at the
same rates) confirms it: Grover equal within CI at every idle rate (n=4: 0.449 vs 0.449 at 0,
0.259 vs 0.257 at 2e-3, 0.092 vs 0.090 at 8e-3; n=5: 0.208 vs 0.211, 0.055 vs 0.055, 0.030 vs
0.030) and BV-10 showing only the tick effect (0.939 → 0.951 at 0, 0.668 → 0.710 at 8e-3).
`noise_idle_sweep_damping.*` (damping only) is the mirror image. This is
Lemma 4 of `proofs/predictions.md` §4: dephasing and depolarizing commute with every Pauli
flag, amplitude damping commutes with Z only, and an X flag swaps which logical state is the
damping fixed point. The same lemma says the effect is circuit-dependent and can go against
the pass: a logically |0⟩ qubit carrying an X flag sits physically in |1⟩ and decays. The
defensible claims are
(i) frame tracking shortens the circuit by a tick wherever an absorbed gate was alone in its
layer, and (ii) under T1-dominated noise it protects the answer by keeping the physical
register in |0⟩; it does not make Grover more robust to dephasing.

Packing more operations into a tick always helps in
this model because it has no crosstalk term; the "best number of operations per tick" asked
for in the plan needs a crosstalk penalty, which is future work.

`shared_oracle_noise.*` repeats the step-5 comparison under this noise model at N = 8 with the
table-lookup oracle (the bucket brigade's 35 qubits are out of statevector range), k = 1…3
instances, sweeping the idle rate: in the shared design the queued instances take idle error
while they wait. See the CSV for the numbers.
