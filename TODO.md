# TODO

Spec: `project.md`, added 14:46 as a byte-for-byte copy of `research-plan.md`. Only the student changes the spec; agents must not edit it.
Maintained by the review agent. Last review: 2026-10-28 14:52.

Agents:
- **Proofs** ("project.md proofs kickoff"): Step 1, the formal model.
- **Setup/compiler** ("project.md setup"): environment, `src/qarch`, Step 2 onward.
- **Review** ("spec compliance review"): checks both against the spec and keeps this file and `research.md` current.

## Status at last review

| Step | Owner | Status | Spec compliance |
|---|---|---|---|
| 0. Environment | Setup | In progress: `git init`, uv project `qarch`, Python 3.12, empty `src/ tests/ benchmarks/`. **Packages still not installed in `.venv`** as of 14:50. Implementation plan written (`~/.claude/plans/sunny-zooming-quail.md`), now including the review notes. **Waiting for the student to approve the plan (plan mode)** before any code. | OK. Deps match the plan exactly. |
| 1. Formal model and proofs | Proofs | Working. No `proofs/` files yet. | Can't assess yet |
| 2. Compiler pass | Setup | Planned (`frame.py`, `passes.py`, `verify.py`) | Plan matches spec |
| 3. Gate-savings benchmark | Setup | Planned (`circuits.py`, `metrics.py`, `benchmarks/gate_savings.py`) | Plan matches spec; baseline O3 and O3+pass |
| 4. Oracle lookup cost | Setup (scaffold only this session) | Planned | Matches spec; adds T-count (fine) |
| 5. Shared-oracle test | Setup (scaffold only) | Planned | Matches spec |
| 6. Noise simulation | Setup (scaffold only) | Planned | Matches spec |
| Data analysis and figures | Setup (CSV and PNG to `results/`) | Planned | Regression is planned; CIs planned for Step 6 |

## 0. Environment (Setup)
- [x] uv project with exactly the plan's tool list
- [ ] `uv sync` finishes and `import qiskit, qiskit_aer, stim` works
- [ ] Record the installed Qiskit version in `research.md` (it decides which baseline passes exist; see research.md §3)
- [ ] Fix all random seeds in generators and simulators
- [ ] Make a first commit. The repo has no commits, so no agent can use worktree isolation and they all write into the shared folder.

## 1. Formal model (Proofs)
- [ ] Define the address table (permutation π, logical→physical) and the flag vector (an X/Z bit per qubit, i.e. a Pauli frame).
- [ ] Prove: a Pauli flag P, or an address swap, followed by any unitary *U* equals *U* followed by the updated frame (P·U = U·(U†PU)). This is valid for every U, but U†PU is a Pauli only when U is Clifford.
- [ ] Propagation tables for **H, S, CNOT** (also S†, CZ, SWAP, X, Y, Z for completeness).
- [ ] Handle Y and Z flags, not just X. Show global phase can be ignored and that a Z flag doesn't change a Z-basis measurement.
- [ ] Measurement rule: an X (or Y) flag flips the classical outcome bit. π relabels which bit is which.
- [ ] **T / diagonal-gate boundary.** Kept unchanged, T turns an X flag into a non-Pauli ((X+Y)/√2, i.e. an S-type correction). But X·T = e^{iπ/4}·T†·X, so a compiler that **rewrites T→T†** (RZ(θ)→RZ(−θ)) keeps the frame Pauli. Prove both and state the class exactly: single-qubit diagonal gates are absorbable with a rewrite; non-diagonal non-Cliffords and multi-qubit ones (CCZ, MCZ) are not, because X on a CCZ control leaves a CZ correction.
- [ ] Cover mid-circuit measurement and classically controlled gates, or explicitly exclude them.
- [ ] Entangling limit: CNOT can't be removed. A SWAP (= 3 CNOT) can be removed by relabeling.
- [ ] Deliverable: files under `proofs/` whose propagation tables are written so Step 2's unit tests can check them line by line.

## 2. Compiler pass (Setup)
- [ ] Qiskit `TransformationPass` that strips X/Y/Z/SWAP and outputs (compressed circuit, π, frame).
- [ ] Post-processing that applies the frame and π to measurement counts.
- [ ] Verification: Stim at large *n*, statevector at *n* ≤ ~16, property tests on random Cliffords.
- [ ] Unit tests generated from the Step 1 propagation tables. Setup's `frame.py` is its "single source of truth"; it must agree with `proofs/`, so reconcile once proofs land.
- [ ] Implement the planned `flush` and `rewrite` policies for non-Clifford gates, and report results for both.

## 3. Gate-savings benchmark (Setup)
- [ ] BV, GHZ, syndrome extraction, and random Clifford circuits at *n* = 4…64. Then QFT, Grover, and Shor (Shor as gate counts only).
- [ ] Measure total gates, 2-qubit gates, and depth, before and after.
- [ ] **Baseline is `transpile(optimization_level=3)`.** Report both "O3" and "O3 + our pass". Disclose any SWAP elision Qiskit already performs (research.md §3).
- [ ] Linear regression of gates saved vs *n* per family (slope, R²).

## 4–6 (need an owner once Step 2 runs)
- [ ] Step 4: bucket-brigade and QROM / table-lookup circuits, *N* = 8…1024. Qubits, gates, and depth per lookup.
- [ ] Step 5: *k* = 1…8 Grover instances, duplicated oracle vs shared queued oracle. Find the break-even point. Cite Zalka.
- [ ] Step 6: amplitude damping and dephasing plus idle error. Vary ops per tick. ≥1000 shots, several seeds, 95% CIs.

## Decisions for the student (flagged by review)
- [ ] **"Savings grow linearly with *n*" is a hypothesis.** BV, GHZ, and syndrome circuits contain few X/SWAP gates. Random Cliffords decomposed into H/S/CX have no explicit Paulis at all unless the generator emits them. Decide whether the Clifford generator draws from a gate set that includes Paulis and SWAPs, and report that choice.
- [ ] **Grover's X gates surround a multi-controlled Z** (non-Clifford for n ≥ 3). Pushing an X through it leaves CZ-type Clifford corrections, so those X's are not removable as Pauli flags. QFT's final SWAP layer *is* removable (⌊n/2⌋ SWAPs, not the whitepaper's n). T/RZ gates *are* passable with the T→T† rewrite. The "clear boundary" sits at multi-qubit and non-diagonal non-Clifford gates, which is broader than the plan's wording ("Clifford") and may be worth a sentence in the final write-up. The spec itself stays unchanged; this is your call.
- [ ] **Novelty framing.** Pauli frame tracking (Knill 2005; Riesebos 2017) and SWAP elision by relabeling are established techniques, and Qiskit and Stim implement versions of them. The contribution should be framed as combining and evaluating them inside this architecture, as the plan's "Gap" paragraph says, not as inventing them.
- [ ] Bibliography: Riesebos 2017, Knill 2005, Kielpinski 2002, and Van Meter & Horsman 2013 are still "cited from memory" in research-directions.md. Verify before submission.
