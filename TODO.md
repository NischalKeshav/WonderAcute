# TODO

Spec: `project.md` (byte-for-byte copy of `research-plan.md`). Only the student changes the spec; agents must not edit it.
Maintained by the review agent. Last review: 2026-10-31. History of each step is in `research.md` §2a–2e.

Agents:
- **Proofs** ("project.md proofs kickoff"): Step 1, `proofs/`. Finished.
- **Code** ("project.md setup"): `src/qarch`, `tests/`, `benchmarks/`, `results/`, on branch `worktree-code-compiler` (worktree `.claude/worktrees/code-compiler`, 16 commits, latest fad563e). Not merged into `master`.
- **Review** ("spec compliance review"): keeps this file and `research.md`. Note: a second session named "spec compliance review (2)" appeared on 2026-10-31; only one should edit these files.

## Status

| Step | Status | Spec compliance |
|---|---|---|
| 0. Environment | Done (uv, Python 3.12, Qiskit 2.6.1, Aer 0.18, Stim 1.16; 264 tests) | OK |
| 1. Formal model and proofs | Done; re-reviewed at the student's request; two checkers pass. **New 10-31: `proofs/predictions.md` + `check_predictions.py`** derive closed forms for Steps 3–5 and match every benchmark row (168 + 48 + 16 + 64) | **Compliant, independently verified** (review re-derived Grover n+z / 2n+2z on 63 rows and the break-even rule on 64 rows) |
| 2. Compiler pass | Done | **Compliant** (feed-forward unsupported; must be stated) |
| 3. Gate-savings benchmark | Done: all-to-all and line sweeps, n = 4…64, 3 seeds | **Compliant** |
| 4. Oracle lookup cost | Done, N = 8…1024 | **Compliant** |
| 5. Shared-oracle test | Done, k = 1…8, N = 8…1024, break-even found | **Compliant** |
| 6. Noise simulation | Done, incl. O3-fair pair and idle sweeps over four channels (combined, dephasing, depolarizing, damping) | **Compliant**; interpretation settled (T1 fixed-point effect, Lemma 4) |
| Data analysis and figures | CSV + PNG per step in `results/`, regression slopes/R², 95% CIs | Done |
| Shared oracle under noise (extra, review-requested) | **Done** (`results/shared_oracle_noise.csv`, verified by review) | Answers the spec's "does a shared oracle still save qubits once lookup cost is counted" under noise |

## Open items for the code agent
- [x] **Noise interpretation (10-31).** The code agent reached the same conclusion independently and `results/README.md` step 6 now states it: under dephasing pass+O3 = O3 within CI for every circuit including Grover; the entire combined-model advantage is the amplitude-damping (T1) asymmetry. `predictions.md` §4 predicted the gate-count effect would be too small to see.
- [x] Damping-only and depolarizing idle sweeps done and verified by review against the CSVs (research.md §2e). They split exactly as Lemma 4 predicts: under depolarizing, Grover is equal within CI at every rate and BV-10 keeps only the tick effect; under damping alone every pair separates (Grover-4 0.224 → 0.326 at p_idle = 2e-3). README step 6 states this with Lemma 4 cited.
- [x] Damping sweep committed (1d36c5b); worktree clean.
- [x] Dephasing sweep files committed (54105dd).
- [ ] Note: sessions other than the proofs agent have edited `proofs/predictions.md` and `check_predictions.py` (§5 working-paper claims, bucket-brigade flush form). The proofs agent kept them. Fine, but one owner per file is safer; student's call.
- [x] Shared-oracle-under-noise done and verified (research.md §2d). Sharing saves 25–33% of qubits with the N = 8 table-lookup oracle and costs queued instances about a third of their success probability at p_idle = 2e-3 (k = 3: 0.171 duplicated vs 0.127 shared), slightly worse even at p_idle = 0 because the shared ancillas take every instance's gate noise.
- [x] Final commit confirmed: fad563e; worktree clean; 16 commits on the branch.

**Code agent status:** every spec item (Steps 2–6) is delivered on `worktree-code-compiler`, 264 tests pass. Future work it lists in `results/README.md`: crosstalk term (for an ops-per-tick optimum), bucket brigade under noise (out of simulation range), select-swap formula and Google 2024 rates to be checked against the papers.
- [x] README line 68 states the shared-oracle costs use the built 2(N−2)-Toffoli table lookup (Babbush's N−1 is in the CSV).
- [ ] **Readout error** (requested by the student 10-31): add Aer ReadoutError (symmetric p_ro, default 1e-2), rerun O3 vs pass+O3 for Grover 3–5 and BV 6/10 at p_idle = 2e-3, plus a p_ro sweep (0, 5e-3, 1e-2, 2e-2). Report whether the T1 advantage survives.
- [ ] **Crosstalk term** (requested by the student 10-31): per tick, each qubit in a 2-qubit gate takes extra depolarizing error p_xt·(m−1), m = number of 2-qubit gates in that tick. Sweep p_xt ∈ {0, 1e-3, 3e-3, 1e-2} × ops-per-tick ∈ {1, 2, 4, 8, unbounded} for Grover-4/5; report the optimal ops-per-tick per p_xt, original vs compressed. Answers research-directions idea 5.

## Decisions for the student
- [ ] **How to present the noise result.** Now backed by a proof: `predictions.md` §4 Lemma 4 (verified by review on the Kraus operators) shows a Pauli frame commutes with dephasing and depolarizing but amplitude damping commutes only with Z; with an X flag on the line the channel's fixed point moves to logical |1⟩, so the T1 effect **can go either way** depending on what the frame leaves on each line. Recommended wording: "Under Pauli-symmetric noise (dephasing, depolarizing) the compressed circuit performs the same as the baseline; the removed single-qubit gates are too few to matter. Under amplitude damping it performs better, because the frame keeps the physical qubits closer to |0⟩, which T1 decay cannot touch (Grover-4: 0.153 → 0.215; BV-10: 0.852 → 0.915 at p_idle = 2e-3)." Do not describe it as a clock-tick effect. For BV add that the answer is computed classically (Gottesman–Knill), so less quantum work is left.
- [ ] **How to restate "measured savings."** Exact per-family formulas are in predictions.md §1 (Grover: n+z native, 2n+2z Clifford+T, z = zero bits of the marked item, which is why R² < 1; QFT and Shor-15: ⌊n/2⌋ SWAPs). Over Qiskit O3: linear savings for Grover (6 gates/qubit, R² 1.00) and random Clifford circuits (~0.8–1/qubit), in single-qubit gates and depth only; zero two-qubit savings; zero for GHZ, syndrome, BV; QFT's ⌊n/2⌋ SWAPs are already removed by O3. The savings survive line-coupling routing (Grover 6.49/qubit).
- [ ] **How to restate the shared-oracle claim.** Closed forms (predictions.md §3, verified): qubits saved = (1 − 1/k)·Q/(a+Q), so 2/3 is a *ceiling* at k = 3 and the working paper's 66% is that ceiling for a bucket-brigade oracle; depth ≈ ×k; sharing wins in qubit·depth iff k < Q·D/(a·L) (bucket brigade N = 64: k < 4.4; table lookup: never). Zalka: a shared oracle only makes sense for k *independent* searches, since k parallel copies of the same search cost √(Nk) total queries. Cite Zalka and Fat-Tree QRAM.
- [ ] **`rewrite_cp` default.** Currently off (plain gate count is the spec's metric). Rewrite-on results kept as `results/gate_savings_*_cprewrite_*` for the virtual-Z-hardware story (+25.8 gates/qubit on native-MCZ Grover). Confirm or overrule.
- [ ] **Merge `worktree-code-compiler` into `master`** when satisfied (your call; the review agent won't merge).
- [ ] **Write-up caveats to carry:** feed-forward unsupported in code; removing a SWAP assumes free relabeling, now quantified in `predictions.md` §6 (with an optimal router the pass can never hurt; with SABRE on a line, pass+O3 is never worse than O3 in CNOT-equivalents on the 216 structured rows but worse on 24/48 `random_layered` rows, worst +402 of ~11,800, and has more total gates in 1 of 264 rows); the "o3+pass on line" column is not a hardware saving; Shor-15 is a compiled, order-aware circuit (Smolin, Smith & Vargo 2013); Grover-Clifford+T uses n−3 ancillas; residual flushes before non-diagonal non-Cliffords can cost more than they save (`flushed` column); the "clear boundary" is native multi-qubit non-Clifford gates, slightly broader than the plan's "Clifford" wording.
- [ ] **Novelty framing.** Pauli frames (Knill 2005; Riesebos 2017) and SWAP elision by relabeling are established and partly implemented in Qiskit/Stim; frame the contribution as combining and evaluating them in this architecture, as the plan's "Gap" paragraph says.

## Sources the student must verify (agents could not)
- [ ] Riesebos et al., "Pauli frames for quantum computer architectures," DAC '17 (no copy).
- [ ] Nielsen & Chuang §4.3 (Toffoli in Clifford+T) and §10.5.2 (Clifford group): section numbers from memory.
- [ ] Knill: the plan cites *Nature* 434, 39–44 (2005), "realistically noisy devices"; the arXiv copy (quant-ph/0410199) is titled "Very Noisy Devices." Confirm they are the same work and cite the version you read.
- [ ] Journal refs not in arXiv copies: Aaronson–Gottesman *PRA* 70, 052328; McKay *PRA* 96, 022330.
- [ ] Low, Kliuchnikov & Schaeffer 2018 (arXiv 1812.00954): the select-swap formula Toffoli = ⌈N/λ⌉ + d(λ−1) is from memory.
- [ ] Google 2024 (arXiv 2408.13687) error rates quoted in `results/README.md` (~3–4e-4 per 1q gate, ~3e-3 per CZ): from memory.
- [ ] Smolin, Smith & Vargo, "Oversimplifying quantum factoring," *Nature* 499 (2013), arXiv 1301.7007: from memory.
- [ ] Still "cited from memory" in `research-directions.md`: Kielpinski 2002; Van Meter & Horsman 2013.
- [ ] If the proofs are cited, add Gottesman 1998 and Aaronson–Gottesman 2004 to the plan's bibliography (spec change, student's call).
