# Research Summary

Maintained by the review agent. It summarizes research done so far and what still needs checking. Last update: 2026-10-28 14:50.

## 1. Where the research stands

| Source | Status |
|---|---|
| `source/…Modern Architecture for Quantum Computers.pdf` | Original working paper: oracle ROM sector, classical controller doing X/SWAP via an address table, packed clock. Its claims were calculated by hand, and it treats oracle lookups as free. |
| `research-plan.md` | ISEF plan; the de facto spec. Scope is research-directions ideas 1, 3, 4, 5. |
| `research-directions.md` | Seven directions with a literature table. 27 papers confirmed against arXiv (2026-10-01); 4 cited from memory. |
| `whitepaper.html` | Diagrams and interactive charts for the original architecture. |
| `project.md` | Added 14:46; identical to `research-plan.md`. This is the spec. |
| Agent work (2026-10-28) | No results yet. The setup agent has a full implementation plan (below); the proofs agent has no files yet. |

## 2. Key findings from the literature (already in the plan)
- **Lookups are not free.** QRAM/QROM lookups cost O(N) gates, or O(N) qubits with O(log N) depth (Jaques & Rattew 2023; Giovannetti et al. 2008; Babbush et al. 2018). This is the largest gap in the original paper.
- **One query at a time.** A bucket-brigade QRAM serves one query at a time; Fat-Tree QRAM (Xu, Lu, Ding 2025) addresses parallel queries. Zalka (1999): parallel Grover instances give only a √k speedup.
- **Pauli frame = the "flag" idea.** Knill 2005 introduced it; Riesebos 2017 describes it at the architecture level. Qubit relabeling for SWAP is standard in mapping/routing (Li, Ding, Xie 2019, SABRE).
- **Zoned hardware** (trapped-ion QCCD, neutral atoms) physically matches the active/inactive split (Thaker 2006; Pino 2020; Bluvstein 2021/2023). This is out of scope for the plan.

## 3. Review notes affecting methodology
These come from the review agent's background knowledge and still need to be checked against the installed package versions and the cited papers.
- **Clifford conjugation is exact:** for P ∈ {X, Y, Z} and Clifford C, C†PC is again a Pauli (up to sign), which is what makes the Step 1 proof work. Standard rules: H swaps X↔Z; S maps X→Y and leaves Z fixed; CNOT copies X from control to target and Z from target to control.
- **T with the gate kept fixed:** T†XT = (X − Y)/√2 is not a Pauli, so the frame would need an S-type Clifford correction.
- **T with a rewrite (correction from the setup agent, checked by hand):** X·T = e^{iπ/4}·T†·X. Because the frame is known at compile time, the compiler can replace T with T† (RZ(θ) with RZ(−θ)) and keep the frame Pauli; only a global phase changes. Z flags commute with any diagonal gate. **The real limit is therefore multi-qubit and non-diagonal non-Clifford gates** (e.g. X on a CCZ control leaves CZ, which is Clifford but not Pauli). This is one degree more permissive than the plan's "Clifford" wording.
- **The baseline may already elide SWAPs.** Recent Qiskit versions (1.x/2.x) include an `ElidePermutations` pass in the preset pass managers at higher optimization levels; it removes SWAPs by relabeling. If it runs at `optimization_level=3`, the SWAP savings over the baseline will be close to zero. That's still an honest result, but it has to be disclosed. **To verify** once `.venv` is installed.
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
- [ ] Verify the 4 "from memory" citations (Knill 2005, Riesebos 2017, Kielpinski 2002, Van Meter & Horsman 2013).
- [ ] Confirm the Qiskit preset passes (ElidePermutations, Optimize1qGatesDecomposition, InverseCancellation) and what they already remove. **Blocked at 14:50: `.venv` still has no packages.**
- [ ] QROM / select-swap cost formulas (Babbush 2018; Low, Kliuchnikov, Schaeffer 2018) as targets for the Step 4 circuits.
- [ ] Amplitude-damping/dephasing parameters for Step 6 (typical T1/T2 values for superconducting hardware) so the noise levels are defensible.
- [ ] Results from the proofs and compiler agents get summarized here as they land.
