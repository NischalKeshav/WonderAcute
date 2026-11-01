# Research Summary

Maintained by the review agent. It summarizes research done so far and what still needs checking. Last update: 2026-10-31.

## 1. Where the research stands

| Source | Status |
|---|---|
| `source/…Modern Architecture for Quantum Computers.pdf` | Original working paper: oracle ROM sector, classical controller doing X/SWAP via an address table, packed clock. Its claims were calculated by hand, and it treats oracle lookups as free. |
| `research-plan.md` | ISEF plan; the de facto spec. Scope is research-directions ideas 1, 3, 4, 5. |
| `research-directions.md` | Seven directions with a literature table. 27 papers confirmed against arXiv (2026-10-01); 4 cited from memory. |
| `whitepaper.html` | Diagrams and interactive charts for the original architecture. |
| `project.md` | Added 14:46; identical to `research-plan.md`. This is the spec. |
| `proofs/` (14:50) | **Step 1 complete.** Formal model, 3 lemmas, 3 theorems, 54-row propagation table, and a checker. Independently verified by review. See §2a. |
| Code (10-31) | Steps 2–6 done (§2b–2e), 10 commits, latest 8687387. Dephasing-only sweep reviewed: **the noise advantage is entirely a T1 effect** (§2e); README corrected. Shared-oracle-under-noise run restarted. |
| `proofs/predictions.md` (10-31) | Closed forms for Steps 3–5, matched to every benchmark row by `check_predictions.py`; independently re-derived by review (§2f). |

## 2. Key findings from the literature (already in the plan)
- **Lookups are not free.** QRAM/QROM lookups cost O(N) gates, or O(N) qubits with O(log N) depth (Jaques & Rattew 2023; Giovannetti et al. 2008; Babbush et al. 2018). This is the largest gap in the original paper.
- **One query at a time.** A bucket-brigade QRAM serves one query at a time; Fat-Tree QRAM (Xu, Lu, Ding 2025) addresses parallel queries. Zalka (1999): parallel Grover instances give only a √k speedup.
- **Pauli frame = the "flag" idea.** Knill 2005 introduced it; Riesebos 2017 describes it at the architecture level. Qubit relabeling for SWAP is standard in mapping/routing (Li, Ding, Xie 2019, SABRE).
- **Zoned hardware** (trapped-ion QCCD, neutral atoms) physically matches the active/inactive split (Thaker 2006; Pino 2020; Bluvstein 2021/2023). This is out of scope for the plan.

## 2a. Step 1 results (proofs agent; verified by review)
- **Theorem 1:** in any circuit of Cliffords, Paulis, SWAPs, measurements (including mid-circuit), resets and feed-forward, every X/Y/Z/SWAP can be removed. Every other gate runs once, unchanged, on relabeled lines, and outcome distributions are identical. This is the plan's expected "proof" outcome, confirmed.
- **Theorem 2:** this extends to Clifford+T and Rz(θ) by substituting T↔T† (Rz(θ)↔Rz(−θ)) when the X flag is set. The T-count is unchanged. With feed-forward, the T/T† choice is made at run time.
- **Theorem 3, the answer to the research question:** the gate-free operations are exactly Paulis times permutations (4ⁿ·n! up to phase). A classical reversible function is gate-free **iff** it is a bit permutation followed by a fixed bit flip: 2ⁿ·n! of the (2ⁿ)! functions (48 of 40,320 for n = 3). CNOT and Toffoli aren't in this class. This is the precise "yes, for a specific class" boundary the plan expected.
- **Where it stops:** only native multi-qubit non-Clifford gates. Toffoli/CCZ with X on a control leaves a CNOT/CZ residue; Toffoli with Z on the target leaves an extra CZ (frame form Z_t·CZ).
- **Native controlled-phase (added 16:00):** X_a·CP(θ)·X_a = CP(−θ)·P_b(θ), so an X entering a native CP leaves a single-qubit phase on the other qubit. That's free on virtual-Z hardware, otherwise one gate. This is relevant to Shor, where the QFT follows X-containing arithmetic. If CP is decomposed to Rz + CNOT, there's no residue (Theorem 2).
- **Rz frame-update residue** is X·Rz(−2θ) (S† is the θ = π/4 case, i.e. T), not "S-type" in general. This is a wording fix.
- **Logic checker** `proofs/check_logic.py` checks the invariant step by step over ~1,000 random branches, including reset and feed-forward, and enumerates Theorem 3 exactly (8/24 maps for n = 2, 48/40,320 for n = 3).
- **QFT prediction:** savings are exactly ⌊n/2⌋ SWAPs (3/2 CNOT per qubit), so linear, with slope ½ SWAP/qubit, not the working paper's n.
- **Correction to the working paper's proof [1]:** "UX is unitary" doesn't show the merged gate is free. Lemma 1 supplies the missing step.
- **Caveats for Steps 2–6:** connectivity (all-to-all assumed), Pauli removal saves no 2-qubit gates, and unmeasured outputs need the frame applied or handed on.

## 2b. Step 3 results: full all-to-all sweep (n = 4…64 step 4, seeds 0–2, 1 Grover iteration)
Source: `results/gate_savings_all*.csv` in the code worktree (written 15:05, **not yet committed**). Review re-derived the key rows from the raw CSV, and they match the code agent's table.

Gates saved per added qubit [R²], with mean saved in parentheses. *pass+o3* (compress, then O3) is the fair comparison; *o3+pass* runs the pass on O3's output.

| Family | pass vs raw | o3+pass vs O3 | pass+o3 vs O3 |
|---|---|---|---|
| Grover, Clifford+T | 4.00 [1.00] (135) | 2.00 [1.00] (66) | **6.00 [1.00] (188)**, exactly 6n − 16 |
| random_layered (**rerun at fixed depth 32**) | 8.11 [0.994] (saved/gates ≈ 0.34) | 2.67 [0.95] | 3.33 [0.96] (saved/gates ≈ 0.21) |
| random_clifford_synth | 0.98 [0.98] (32) | 1.02 [0.96] (31) | 0.78 [0.96] (24) |
| QFT | 0.50 [1.00] (17) = ⌊n/2⌋ SWAPs | **0** | **0** |
| Bernstein–Vazirani | 0 (constant 1) | 0 (1) | 0 (1) |
| GHZ, repetition syndrome | 0 | 0 | 0 |
| Grover, native MCZ (n = 4–10 only) | 2.00 [0.99] (13) | −0.92 [0.61], noise | 0.06 [0.01], noise |
| Shor-15 (compiled, order-aware; n_count 8/12/16) | 0.5/n_count (inverse-QFT SWAPs) | 0 | 0 |

**Two-qubit gates:** compared with O3, at most 5 CX-equivalents are saved in any row. Pauli removal saves essentially no two-qubit gates, and SWAP savings are already captured by O3's ElidePermutations.

**What the data says about the plan's "measured savings" outcome:**
1. **Supported, in a narrower form.** Savings over Qiskit O3 grow linearly with n for Grover (6 gates/qubit, R² = 1.00) and for random Clifford circuits (~0.8–1 gate/qubit). They come from Pauli absorption through Cliffords plus the T→T† rewrite, which O3 doesn't do. They are **single-qubit gates and depth**, not two-qubit gates.
2. **Not supported** for GHZ, syndrome extraction, or BV (nothing to remove), or for QFT once compared with O3.
3. **The working paper's QFT claim is overstated twice over:** the saving is n/2 SWAPs, not n, and O3 already gets it.
4. **random_layered** was quadratic in the first sweep (depth = 2n). It was rerun at fixed depth 32 and is now linear (row above). The pass removes about 21% of the gates left after O3.
5. Native-MCZ Grover on top of O3 is noise: O3 resynthesizes the MCZ differently for each n, and the pass must flush before `cp` gates.
6. **Line coupling (hardware-fair: compress, then route with O3):** Grover-Clifford+T saves 6.49 gates/qubit (R² 0.997), random_layered 3.52 (R² 0.69), random_clifford_synth 0.79 (R² 0.96), native-MCZ Grover 1.17 (R² 0.76); QFT, Shor, GHZ, BV and syndrome save 0. The savings survive routing. The "o3+pass on line" column deletes routing SWAPs (up to 260 gates/qubit) and is **not** a hardware saving.

## 2c. Step 4 results: oracle lookup cost per query (d = 1 data bit; `results/oracle_cost.csv`)

| N | Bucket brigade: qubits / CSWAP / depth | Table lookup (unary iteration): qubits / Toffoli / depth |
|---|---|---|
| 8 | 35 / 72 / 65 | 6 / 12 / 23 |
| 64 | 262 / 732 / 125 | 12 / 124 / 223 |
| 1024 | 4106 / 12,236 / 205 | 20 / 2044 / 3586 |

- **Bucket brigade:** qubits ≈ 4N + k; CSWAP = 12N − 4k − 12 (each CSWAP = 1 Toffoli = 7 T); depth ≈ 20·log₂N. O(N) qubits buy O(log N) depth.
- **Table lookup:** qubits = 2k (k − 1 ancillas); Toffoli = 2(N − 2) as built (Babbush's measurement-uncompute version: N − 1, i.e. 4(N − 1) T); depth O(N).
- **Select-swap** (closed form, unverified): Toffoli ⌈N/λ⌉ + d(λ − 1), λ ≈ √(N/d). This interpolates between the two designs.
- **This confirms the plan's central concern: the oracle isn't free.** At N = 1024 one lookup costs either ~4,100 qubits (bucket brigade) or ~2,000 Toffolis ≈ 14,000 T gates with depth ~3,600 (table lookup). The working paper treated this as zero.
- **Architecture finding:** the pass doesn't make lookups cheaper. In the table lookup every absorbed X is re-emitted before the next Toffoli; in the bucket brigade ~28% of X gates survive absorption (12,754 absorbed, 9,180 flushed at N = 1024), and those are pairs that cancel between the descent and the ascent. This is consistent with proofs §7: an X on a Toffoli control leaves a residue.

## 2d. Step 5 results: shared vs duplicated oracle (`results/shared_oracle.csv`)
Model: k Grover instances, each running ⌊π/4·√N⌋ iterations (25 at N = 1024). Duplicated = one oracle per instance, running in parallel. Shared = one oracle, with lookups queued.

| Oracle, N = 1024 | k | Qubits: duplicated → shared | Depth ratio | Qubit·depth ratio |
|---|---|---|---|---|
| Bucket brigade | 3 | 12,318 → 4,126 (−66.5%) | ×2.81 | **0.94** |
| Bucket brigade | 8 | 32,848 → 4,176 (−87%) | ×7.35 | 0.93 |
| Table lookup | 3 | 60 → 40 (−33%) | ×2.99 | **1.99** |
| Table lookup | 8 | 160 → 90 (−44%) | ×7.96 | 4.48 |

Break-even (review re-derived it from the CSV; "wins" means shared qubit·depth < duplicated):
- **Bucket brigade:** never at N = 8. k = 2–4 at N = 64. All k = 2–8 at N ≥ 256. The gain is at most ~6–7%.
- **Table lookup:** never, for any N or k (qubit·depth 1.7–2× worse at k = 3), because the oracle is only half the qubits.

What this means for the working paper's shared-oracle claim:
- **The "66% fewer qubits at k = 3" claim holds for qubit count**, but only with a qubit-heavy oracle (bucket brigade, N ≥ 64).
- **It costs time:** queued lookups make the run k× longer (×2.8 at k = 3). In qubit·depth terms sharing is roughly break-even for the bucket brigade and a clear loss for the table lookup.
- **Idle error decides it.** Waiting instances sit idle for ~k× longer, so the Step 6 noise model (idle error on waiting qubits) is what determines whether the qubit saving is real. This is consistent with the bucket-brigade bottleneck noted by Jaques & Rattew 2023 and Xu, Lu, Ding 2025 (Fat-Tree QRAM, which exists to serve parallel queries).

**Shared oracle under noise (10-31, `results/shared_oracle_noise.csv`; verified by review).** N = 8 table-lookup phase oracle, k Grover instances, 2 iterations, exact Clifford+T; p1 = 1e-3, p2 = 1e-2; 1000 shots × 5 seeds; CIs ±0.009–0.012. Per-instance success at p_idle = 0 / 1e-3 / 2e-3 / 4e-3:

| k | Design | Qubits | Ticks | Success |
|---|---|---|---|---|
| 1 | (same circuit) | 6 | 291 | 0.230 / 0.192 / 0.174 / 0.148 |
| 2 | duplicated | 12 | 291 | 0.251 / 0.210 / 0.194 / 0.169 |
| 2 | shared | 9 | 530 | 0.233 / 0.169 / 0.128 / 0.111 |
| 3 | duplicated | 18 | 291 | 0.243 / 0.198 / 0.171 / 0.156 |
| 3 | shared | 12 | 786 | 0.218 / 0.147 / 0.127 / 0.109 |

Reading: with this small oracle sharing saves 25–33% of the qubits and costs the queued instances about a third of their success probability at p_idle = 2e-3, more as the idle rate rises; even at p_idle = 0 the shared design is slightly worse, because the shared output/ancilla qubits take gate noise from every instance's lookup. **Answer to the spec's question:** a shared oracle still saves qubits once lookup cost is counted, at a depth and success-probability price that grows linearly with k. The bucket-brigade case (where sharing saves the most qubits) is out of noisy-simulation range and is listed as future work.

## 2e. Step 6 results: success probability under noise (`results/noise.csv`)
Noise model: amplitude damping + dephasing; p1 = 1e-3 (1q gate), p2 = 1e-2 (2q gate), p_idle = 2e-3 (idle qubit per tick); no crosstalk. 1000 shots × 10 seeds, ±95% CI. "Original" is transpiled at **O0**.

| Circuit | Original → compressed, unbounded ops/tick | Ticks | At 1 op/tick |
|---|---|---|---|
| Grover n=3 | 0.621 ±0.011 → **0.680 ±0.012** | 57 → 53 | |
| Grover n=4 | 0.131 ±0.007 → **0.197 ±0.009** | 193 → 187 | 0.080 → 0.136 |
| Grover n=5 | 0.021 ±0.004 → **0.047 ±0.004** | 417 → 409 | |
| BV n=6 | 0.955 ±0.004 → 0.966 ±0.004 (CIs touch) | 6 → 5 | |
| BV n=10 | 0.853 → **0.916** | 10 → 9 | |

- **Supports the plan's "system effects" outcome:** compressed circuits succeed more often, with non-overlapping CIs in every case except BV-6. The relative gain grows as noise dominates (Grover-5 more than doubles, 0.021 → 0.047).
- **Caveat 1, resolved:** fair O3 pair (same seeds and shots), success probability ±95% CI:

| Circuit | O3 | pass+O3 | Gates (O3 → pass+O3) | Ticks |
|---|---|---|---|---|
| Grover n=3 | 0.661 ±0.011 | **0.706 ±0.012** | | 49 → 47 |
| Grover n=4 | 0.153 ±0.008 | **0.215 ±0.008** | 288 → 283 | 188 → 190 |
| Grover n=5 | 0.024 ±0.003 | **0.049 ±0.005** | 629 → 594 | 425 → 413 |
| BV n=6 | 0.957 ±0.004 | 0.968 ±0.004 (CIs touch) | | 6 → 5 |
| BV n=10 | 0.852 ±0.007 | **0.915 ±0.005** | 23 → 22 | 10 → 9 |

**The noise advantage survives the fair baseline.** Idle sweep (O3 → pass+O3, p_idle = 0 / 5e-4 / 1e-3 / 2e-3 / 4e-3 / 8e-3):
- Grover-4: 0.295→0.344, 0.248→0.304, 0.211→0.269, 0.153→0.215, 0.093→0.146, 0.046→0.087
- Grover-5: 0.099→0.135, 0.059→0.095, 0.040→0.074, 0.024→0.049, 0.012→0.031, 0.009→0.024
- BV-10: 0.925→0.937, 0.902→0.928, 0.878→0.924, 0.852→0.915, 0.797→0.900, 0.684→0.857

At p_idle = 0 the gain comes from removed gates. It grows with idle error.

**BV-10 diagnostic (resolved, 16:45).** At p_idle = 8e-3: O3 0.684, pass+O3 0.857, pass+O3 with the removed tick re-inserted 0.815. **The clock tick is worth ~4 pp.** Channel split: dephasing only 0.906 → 0.929 (the tick effect); amplitude damping only 0.745 → 0.915. Cause: the pass turns the ancilla X into a Z flag, the CNOTs copy it onto the controls, and the final H turns it into X flags. **The secret is carried entirely by the classical frame,** the physical data qubits finish in |0…0⟩, and amplitude damping can't touch them.

**Interpretation (review):** for fully Clifford circuits with a compile-time frame, the controller computes the measured answer classically (Gottesman–Knill; consistent with Theorem 1). BV's noise advantage therefore means less quantum work is left, not that the same quantum work is more robust. It also relies on T1-biased noise and the absence of readout error. **Grover (Clifford+T), whose answer stays in the physical state, is the meaningful noise result:** +0.05 at p_idle = 0, rising with idle error.

Simulator check: MPS and statevector agree exactly for the shared-oracle k = 1 configuration (0.174 ±0.011).

**Dephasing-only sweep (reviewed 10-31, `results/noise_idle_sweep_dephasing.csv`; same p1, p2, shots and seeds, channel = dephasing).** O3 vs pass+O3, p_idle = 0 / 2e-3 / 8e-3:
- Grover-4: 0.664 vs 0.664 / 0.557 vs 0.555 / 0.341 vs 0.337
- Grover-5: 0.458 vs 0.460 / 0.247 vs 0.253 / 0.068 vs 0.070
- BV-10: 0.963 vs 0.972 / 0.945 vs 0.954 / 0.906 vs 0.929

**Revised interpretation.** Under dephasing the compressed circuit is indistinguishable from the baseline for Grover and gains ≤ 2 pp for BV (the one removed tick). Therefore the whole advantage in the combined model, for Grover as well as BV, comes from amplitude damping. The reason is general: dephasing and depolarizing channels are invariant under Pauli conjugation (X·Z·X = −Z leaves the channel unchanged), so a Pauli frame cannot change their effect; only removed gates and ticks can, and 5 of 288 gates is too few to see. Amplitude damping is not Pauli-symmetric: it drives |1⟩ → |0⟩. In Grover's diffusion step the original circuit's X layer holds the register near |1…1⟩ during the long multi-controlled-Z decomposition; after the pass absorbs those X gates the physical register sits near |0…0⟩, where damping does nothing. The earlier statement that Grover's gain is "gate-count driven" was wrong.

**Channel split (10-31, `results/noise_idle_sweep_{damping,depolarizing}.csv`, O3 → pass+O3, p_idle = 0 / 2e-3 / 8e-3; same rates, 1000 × 10; * = CIs disjoint):**

| Circuit | Depolarizing (Pauli-symmetric control) | Amplitude damping only |
|---|---|---|
| Grover-4 | 0.449→0.449, 0.259→0.257, 0.092→0.090 | 0.401→0.472*, 0.224→0.326*, 0.058→0.141* |
| Grover-5 | 0.208→0.211, 0.055→0.055, 0.030→0.030 | 0.169→0.238*, 0.032→0.090*, 0.005→0.026* |
| BV-10 | 0.939→0.951*, 0.853→0.878*, 0.668→0.710* (tick only) | 0.956→0.972*, 0.898→0.946*, 0.745→0.915* |

**Settled conclusion for the plan's "system effects" outcome.** Under Pauli-symmetric noise (dephasing, depolarizing) the compressed circuit performs the same as the O3 baseline, apart from BV's one saved clock tick; under amplitude damping it performs better, often by a factor of 2–3 at high idle rates, because the frame keeps the physical register near |0⟩ (Lemma 4, `predictions.md` §4). This is a genuine property of frame tracking under T1-biased noise, but it is circuit-dependent (an X flag left on a logical |0⟩ parks it in |1⟩, where it decays) and the model has no readout error. `results/README.md` step 6 now says exactly this.
- **Caveat 2:** with no crosstalk, more ops per tick is always better (Grover-4 rises from 0.080 at 1 op/tick to 0.131 unbounded). The model can't identify an optimal ops per tick.
- **Caveat 3:** the noise rates are deliberately ~3× worse than Google 2024 (arXiv 2408.13687: ~3–4e-4 per 1q gate, ~3e-3 per CZ, according to the code agent; not checked by review), and a p_idle sweep (0 … 8e-3) is running.

**CP rewrite trade-off (commit 912fe69):** a pending X on a controlled-phase gate can become CP(−θ) plus P(θ) on the other qubit. On O3 output for native-MCZ Grover this costs **26.7 gates/qubit** in the plain count but saves **25.8/qubit** when diagonal 1q gates are free (virtual-Z, McKay 2017). Grover-Clifford+T: 1.0 vs 7.0/qubit. The fair pass+O3 column is unaffected. So the architecture's value depends on whether Z rotations are free on the hardware. That's worth one sentence in the write-up.

## 2f. Closed-form predictions (proofs agent, 10-31; verified by review)
`proofs/predictions.md` derives the Step 3–5 numbers from the flag rules, and `check_predictions.py` matches them to the CSVs row by row (168 gate-savings rows, 48 closed-form and 16 built lookup rows, 64 shared-oracle rows; all pass on review's run). Review separately re-derived the Grover and shared-oracle formulas from the raw CSVs with its own code: 0 mismatches.

- **Step 3, exact savings (pass vs raw):** GHZ, syndrome 0; BV 1; QFT and Shor-15 ⌊n/2⌋ SWAPs; Grover native MCZ **n + z**, Grover Clifford+T **2n + 2z**, where z = number of zero bits of the marked item. The regression's R² < 1 for Grover comes only from z varying with the seed (z = n, n−1, n−1 for seeds 0–2), so the "4 gates/qubit" slope is exact at fixed z. O3 removes the QFT/Shor SWAPs itself, so the pass adds 0 there and the Grover savings survive O3.
- **Step 4:** bucket brigade CSWAP = 4(N−1−k) + 8d(N−1), qubits = k + d + 2(N−1) + N + Nd + d; the pass flushes 4(N−1−k) + 4d(N−1) + 2dw X gates (w = set table bits), tending to 2/3 of them, and removes **no CSWAP**, so the Toffoli/T cost is untouched. Table lookup: every X is flushed, the pass saves nothing. Select-swap: N/λ + dλ is minimized at λ = √(N/d) with value 2√(Nd) (the Toffoli formula itself is still from memory).
- **Step 5:** qubits saved = **(1 − 1/k)·Q/(a + Q)** (Q = oracle qubits, a = address bits), so 2/3 is a ceiling at k = 3 and the working paper's 66% is that ceiling for a bucket brigade; the checker also reproduces the whitepaper's 61/64/65/66% table as 1 − (N + 3 log N)/(3N). Depth ratio (kL + D)/(L + D) ≈ k. **Shared wins in qubit·depth iff k < Q·D/(a·L):** bucket brigade N = 64 gives k* = 4.4, N = 8 gives 1.15 (never), table lookup k* < 1 (never). Zalka: k parallel copies of the *same* search need √(Nk) total queries, so sharing is only sensible for independent searches.
- **Step 6:** removing g of G gates can at most multiply success by ≈ e^(gε); for Grover that is ~6% of gates, too small to see at p1 = 1e-3. This is consistent with the dephasing-only result in §2e.
- **Lemma 4 (added 10-31, verified by review):** for every Pauli P, P·𝒩·P = 𝒩 when 𝒩 is dephasing or depolarizing, so a Pauli frame cannot change their effect; amplitude damping satisfies this only for P = Z, and X·𝒩_AD·X has fixed point |1⟩. This is the proof behind the §2e noise result and shows the T1 gain is a fixed-point asymmetry that depends on which flag the frame leaves on a line.
- **Depth (added 10-31, 165 rows):** raw − pass depth is exactly BV 1, QFT 1, GHZ/syndrome 0, Grover native 2, Grover Clifford+T 4 if qubit 1 is X-wrapped else 2. Review confirmed from the CSV.
- **Connectivity (added 10-31, from `gate_savings_line.csv`):** with an optimal router the pass can never hurt, since the original SWAPs are one valid routing of the relabeled circuit; with SABRE on a line it can. pass+O3 vs O3 in CNOT-equivalents: never worse on the 216 structured rows; worse on 24/48 `random_layered` rows (worst +402 of ~11,800). Total gates are lower in 263 of 264 rows (exception: `random_clifford_synth` n = 4, seed 1, 31 → 33 gates, equal CNOTs). This makes `formal-model.md` §8 caveat 1 quantitative.

## 3. Review notes affecting methodology
These come from the review agent's background knowledge and still need to be checked against the installed package versions and the cited papers.
- **Clifford conjugation is exact:** for P ∈ {X, Y, Z} and Clifford C, C†PC is again a Pauli (up to sign), which is what makes the Step 1 proof work. Standard rules: H swaps X↔Z; S maps X→Y and leaves Z fixed; CNOT copies X from control to target and Z from target to control.
- **T with the gate kept fixed:** T†XT = (X − Y)/√2 is not a Pauli, so the frame would need an S-type Clifford correction.
- **T with a rewrite (correction from the setup agent, checked by hand):** X·T = e^{iπ/4}·T†·X. Because the frame is known at compile time, the compiler can replace T with T† (RZ(θ) with RZ(−θ)) and keep the frame Pauli; only a global phase changes. Z flags commute with any diagonal gate. **The real limit is therefore multi-qubit and non-diagonal non-Clifford gates** (e.g. X on a CCZ control leaves CZ, which is Clifford but not Pauli). This is one degree more permissive than the plan's "Clifford" wording.
- **The baseline already elides SWAPs: confirmed.** Qiskit 2.6.1 O3 runs `ElidePermutations`. The quick run shows QFT savings vs O3 = 0 (see §2b).
- **Stim tracks Pauli frames natively** (its frame simulator), so it's a good independent check for Step 2.
- **Grover/QFT/Shor extension.** QFT's terminal SWAPs are absorbable. Grover's diffusion X gates surround a multi-controlled Z, so they aren't Pauli-absorbable. Shor's modular-arithmetic X gates are mixed. Expect partial savings.

## 4. Setup agent's implementation plan (reviewed 14:50, matches spec)
- `src/qarch/frame.py`: address table (permutation plus x/z bits) with Clifford update rules (H: x↔z; S: z^=x; CNOT: x_t^=x_c, z_c^=z_t; CZ, SX, Sdg).
- `passes.py`: `AddressTablePass` removes X/Y/Z/SWAP and re-targets gates; for non-Cliffords, uses the `flush` or `rewrite` policy; stores measurement flips in `property_set`.
- `verify.py`: Operator/Statevector at n ≤ 16, Stim tableaux up to n = 64; ≥200 random circuits.
- `circuits.py`: BV, GHZ, repetition/surface-code syndrome (`stim.Circuit.generated`), `random_clifford`, QFT, Grover, Shor N = 15/21.
- `benchmarks/`: `gate_savings.py`, `oracle_cost.py`, `shared_oracle.py`, `noise.py` writing CSV and PNG to `results/`.
- Built-in sanity check: QFT saves exactly ⌊n/2⌋ SWAPs.
- Review note: `qiskit.quantum_info.random_clifford` synthesizes to H/S/CX, so it will contain almost no explicit Paulis or SWAPs. Savings on that family will be near zero unless the generator also samples Pauli/SWAP layers. Flagged in TODO.

## 5. Still to research
- [x] Proofs §9 sources checked against arXiv (15:05, copies in `proofs/sources/`): Gottesman 1998 (Table 1 = our H/S/CNOT rows), Aaronson–Gottesman 2004 (§III bit rules = `frame.py`), Knill quant-ph/0410199 (defines the Pauli frame; notes that non-Clifford gates may need "non-Pauli compensation"), McKay 2017 (virtual Z by phase offset: related to, but not the same mechanism as, Theorem 2's Rz rule), Li–Ding–Xie 2019 (SWAP = 3 CNOT; connectivity). Review spot-checked these.
- [ ] Still unchecked: Riesebos DAC '17, N&C §4.3 and §10.5.2, Kielpinski 2002, Van Meter & Horsman 2013. The student needs to get copies.
- [ ] Knill title mismatch: the arXiv title is "Very Noisy Devices" but the plan cites the Nature title "Realistically Noisy Devices". Likely the same work; unconfirmed.
- [ ] Journal references for Aaronson–Gottesman and McKay come from memory (not in the arXiv copies).
- [ ] Verify Smolin, Smith & Vargo, "Oversimplifying quantum factoring", Nature 499 (2013), arXiv 1301.7007 (from memory), if the Shor caveat is cited.
- [x] ElidePermutations in O3 confirmed (Qiskit 2.6.1).
- [x] Babbush 2018 unary-iteration count (N − 1 Toffoli, 4(N − 1) T) carried in the CSV.
- [ ] Low, Kliuchnikov, Schaeffer 2018 (arXiv 1812.00954): verify the select-swap formula. The review agent's download was blocked; the student needs to check it.
- [ ] Source the Step 6 parameters (p1 = 1e-3, p2 = 1e-2, p_idle = 2e-3) against published hardware figures (e.g. Google 2024, arXiv 2408.13687), or present them as a sweep.
- [ ] Results from the proofs and compiler agents get summarized here as they land.
