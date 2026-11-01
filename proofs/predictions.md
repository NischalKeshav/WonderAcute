# Closed-Form Predictions for Steps 3–5

The research plan's data analysis fits gates saved against *n* and reports slopes. This document derives those slopes (and the Step 4 lookup costs and the Step 5 break-even point) from the circuit definitions, so the benchmark numbers can be checked against a prediction rather than read off a plot. Everything here follows from Theorems 1–2 and the boundary table in `formal-model.md` §7, applied to the circuits in `src/qarch/circuits.py` and `src/qarch/qram.py`.

`proofs/check_predictions.py` recomputes every formula below from scratch and compares it with the benchmark CSVs row by row. As of 2026-10-31 every row matches (168 gate-savings rows, 48 closed-form lookup rows, 16 built lookup rows including the bucket-brigade flush count, 64 shared-oracle rows), and §5 checks the working paper's own numbers.

Notation: *n* qubits, *z* = number of qubits the Grover oracle wraps in X (the zero bits of the marked item), *N* = 2ᵏ oracle entries of *d* bits, *k* also denotes the number of Grover instances in §3 (as in the plan).

---

## 1. Step 3: gates removed by the pass, family by family

"Pass vs raw" means the address-table pass applied to the textbook circuit, all-to-all coupling, no other optimization. Measurements are not counted as gates.

| Family | Removable gates | 2-qubit gates removed | Slope in *n* | Why |
|---|---|---|---|---|
| GHZ | 0 | 0 | 0 | H and CNOT only; nothing to absorb |
| Repetition-code syndrome | 0 | 0 | 0 | CNOT, measure, reset only |
| Bernstein–Vazirani | 1 | 0 | 0 | the ancilla's X (prepares \|−⟩); the H after it turns the X flag into a Z flag, which the following CNOT targets ignore (Theorem 1) |
| QFT | ⌊*n*/2⌋ | ⌊*n*/2⌋ SWAPs | ½ gate, ½ SWAP, 3/2 CNOT-equivalent | final bit-reversal SWAP network absorbed by π (Lemma 2); CP(θ) gates carry no flags |
| Shor-15, *n* counting qubits | ⌊*n*/2⌋ | ⌊*n*/2⌋ SWAPs | ½ | inverse-QFT SWAPs. The work register's X is flushed: it meets a controlled-SWAP whose target it is, which is non-Clifford (§7) |
| Grover, 1 iteration, MCZ decomposed to Clifford+T | 2*z* + 2*n* | 0 | 4 (for *z* = *n*) | every X passes the H/CNOT/T gates of the Toffoli chain (Theorem 2); T becomes T† where needed |
| Grover, 1 iteration, MCZ native | *n* + *z* | 0 | 2 (for *z* = *n*) | see below |

**Grover with a native MCZ** is the one case that needs the boundary table. One iteration is

  H⊗ⁿ · [X^z · MCZ · X^z] · [H⊗ⁿ · X⊗ⁿ · MCZ · X⊗ⁿ · H⊗ⁿ] · measure.

Walk the flags through it:
1. Oracle, X before MCZ (*z* gates): absorbed, then flushed in front of the MCZ because an X on an MCZ qubit leaves a non-Pauli residue (§7). *z* gates out, *z* gates back: saves 0.
2. Oracle, X after MCZ (*z* gates): absorbed. The next H turns them into Z flags.
3. Diffusion, X before MCZ (*n* gates): absorbed. On the *z* qubits that already carry Z, the flag is now Y; on the other *n*−*z* it is X. The MCZ forces a flush: **one** gate per qubit (Y or X), *n* gates in all, where the original circuit had *z* + *n*. Saves *z*.
4. Diffusion, X after MCZ (*n* gates): absorbed; the final H turns them into Z flags, and Z flags have no effect on measurement (Lemma 3). Saves *n*.

Total saved *n* + *z*, with 2*z* + 2*n* absorbed and *z* + *n* flushed. Depth goes from 9 layers to 7. With marked item 0 (*z* = *n*) that is 2*n* out of 7*n* + 2 gates.

**What this means for the plan's regression.** Every family's savings are exactly linear with the slopes above and R² = 1 when *z* is fixed; the benchmark's R² < 1 for Grover comes only from *z* varying with the seed (seeds 0, 1, 2 give *z* = *n*, *n*−1, *n*−1, so the mean saved for Grover Clifford+T is 4*n* − 4/3, which is what the regression file reports). The savings are all single-qubit gates except the SWAPs in QFT and Shor, which matches caveat 2 of `formal-model.md` §8.

**Depth.** The pass removes whole layers only where the removed gates sit on the critical path: BV −1 (the ancilla X), QFT −1 (the SWAP layer), GHZ and syndrome 0, Grover native-MCZ −2 (9 layers to 7: two flushes replace four X layers), Grover Clifford+T −4 when the first Toffoli's inner control (qubit 1) is wrapped in X and −2 when it is not, because the chain's critical path enters through that qubit (`cx(b, c)` is the Toffoli's first two-qubit gate). Shor-15's measured saving (*n* − 2) depends on how the SWAP layer synchronised the two registers and is not derived here.

**Against Qiskit's optimizing compiler (the plan's fair baseline).** Qiskit O3 already removes the QFT SWAP network by relabeling, so the predicted *additional* saving of the pass after O3 is 0 for QFT and Shor-15, and the benchmark finds exactly that. O3 does not track a Pauli frame, so the Grover savings survive.

## 2. Step 4: memory-lookup cost

### 2a. Bucket brigade (qubit-router model)

A tree with *N* − 1 internal nodes, each holding a router qubit and an input qubit; *N* leaves hold the *d*-bit cells. Level *j* (root = level 0) has 2ʲ nodes.

- *Routing the address.* Address bit *j* must travel from the root to level *j*, through the 2ʲ − 1 nodes above that level, with 2 CSWAPs per node (one per child, selected by an X pair on the router). Summed over *j* = 0…*k*−1 that is 2·Σ(2ʲ − 1) = **2(N − 1 − k)** CSWAPs down, and the same again to uncompute: 4(N − 1 − k).
- *Reading the data.* A one-qubit bus descends through all *N* − 1 nodes (2 CSWAPs each), the addressed cell is copied with CNOTs, the bus ascends, and this is done twice per data bit (once to copy, once to clean up): 4·*d*·2(N − 1) = 8*d*(N − 1) CSWAPs.
- *Totals.* CSWAP = **4(N − 1 − k) + 8d(N − 1)**, which is 12*N* − 4*k* − 12 for *d* = 1. One Toffoli per CSWAP, so T = 7·CSWAP with the exact Toffoli. Data CNOTs = 2*N d*. Qubits = *k* + *d* + 2(*N* − 1) + *N* + *N d* + *d* (address, output, routers and inputs, leaf selects, cells, bus). At *N* = 1024, *d* = 1: 12,264 CSWAPs and 4,106 qubits per lookup.

**What the pass does to it.** Every CSWAP comes with an X on its router, so there are as many structural X gates as CSWAPs, plus *w* data-preparation X gates (*w* = set bits in the table). By §7 an X on a CSWAP's control is not absorbable past it. Follow one router at depth δ through the circuit: it sees *m* = *k* − 1 − δ address passes in phase 1 and again in phase 3, and *d* descent/ascent cycles in phase 2. Each descent is `C X C X` and each ascent `X C X C` (C = CSWAP on this router, X = X on it). Scanning that string, an X is flushed when the next event on the router is a C, and two adjacent X's cancel. The only adjacent pairs are where a descent's trailing X meets an ascent's leading X and vice versa, 4 per data bit. So the router flushes 4*m* + 4*d* of its 4*m* + 8*d* X gates. Summing 2^δ routers per depth and using Σ 2^δ(*k* − 1 − δ) = *N* − 1 − *k*:

  **flushed X = 4(N − 1 − k) + 4d(N − 1) + 2dw**,  **removed gates = X + SWAP − flushed = 4d(N − 1) + 2(N − 1 + k) + 4d + w − 2dw**,

where the 2(*N* − 1 + *k*) + 4*d* SWAPs (address loading, router setting, bus loading) are all absorbed by π, and each set table bit costs 2*d* flushes because the CNOT copies its X flag onto the leaf, which then meets a CSWAP on the way up. For *d* = 1 the pass removes 6*N* + 2*k* − 2 − *w* of about 28*N* gates, and the flushed fraction of X gates tends to 2/3 (the benchmark's 67–73%). None of the removed gates is a CSWAP, so the Toffoli and T counts, which dominate the lookup cost, are unchanged. This is the exact form of the "the pass does not make lookups cheaper" finding. The checker verifies both formulas on every built row.

### 2b. Table lookup (unary iteration)

The recursion visits every internal node below the root once (*N* − 2 of them), doing an X pair and a compute/uncompute Toffoli pair, plus an X pair at the root: **2(N − 2) Toffolis** and **2(N − 1) X gates**, with *k* − 1 ancillas and up to *N d* data CNOTs. (With measurement-based uncomputation, the literature's count is *N* − 1 Toffolis and 4(*N* − 1) T gates; that figure is quoted from memory of Babbush et al. 2018 and is **not** verified here.)

**Prediction: the pass saves nothing on a table lookup.** Each X in the recursion is immediately before or after a Toffoli that uses the same qubit as a control. An X flag on a Toffoli control is not absorbable (§7), so every one of the 2(*N* − 1) X gates is flushed. The benchmark confirms this at every *N*: gates after the pass equal gates before.

### 2c. Select-swap

With chunk size λ (a power of two), a select over ⌈*N*/λ⌉ chunks plus a swap network over λ registers of *d* bits costs **⌈N/λ⌉ + d(λ − 1) Toffolis** and *k* + *d*λ + (*k* − 1) qubits. The formula is from memory of Low, Kliuchnikov & Schaeffer 2018 and is listed in `TODO.md` for the student to verify. What *is* proved here: *N*/λ + *d*λ is minimized at **λ\* = √(N/d)** with minimum **2√(Nd)**, so select-swap interpolates between the table lookup (λ = 1, cost ≈ *N*) and a pure swap network (λ = *N*, cost ≈ *dN*), and its Toffoli count grows as √*N* instead of *N* at the price of √(*N d*) extra qubits.

## 3. Step 5: one shared oracle versus *k* copies

Let *a* = log₂*N* (address bits per instance), *Q* = oracle qubits beyond the address register, *L* = depth of one lookup, *D* = 2*a* + 1 = depth of one diffusion step, and *I* = ⌊(π/4)√*N*⌋ iterations. The benchmark's model is

  duplicated: qubits *k*(*a* + *Q*), depth *I*(*L* + *D*);  shared: qubits *ka* + *Q*, depth *I*(*kL* + *D*).

**Qubits saved.** 1 − (*ka* + *Q*) / (*k*(*a* + *Q*)) = **(1 − 1/k) · Q/(a + Q)**. So sharing can never save more than 1 − 1/*k* of the qubits (2/3 for *k* = 3), and it approaches that ceiling only when the oracle dwarfs the address register. For the bucket brigade at *N* = 1024, *Q* = 4096 and *a* = 10, so the saving is (2/3)(4096/4106) = 66.5%: this is where the working paper's "66% for three instances" comes from, and it is the ceiling, not a typical value. For a table lookup *Q* = *a* = 10, and the same three instances save only (2/3)(1/2) = **33%**.

**Depth paid.** shared/duplicated = (*kL* + *D*)/(*L* + *D*), which is close to *k* whenever the lookup is deeper than the diffusion step, which it is for every design and *N* in the benchmark.

**Break-even in qubit·depth.** Shared wins iff (*ka* + *Q*)(*kL* + *D*) < *k*(*a* + *Q*)(*L* + *D*). Expanding and cancelling the common terms *kaD* and *kQL*:

  *k*²*aL* + *QD* < *kaL* + *kQD*  ⟺  *k*(*k* − 1)*aL* < (*k* − 1)*QD*  ⟺  **k < Q·D / (a·L)**  (for *k* ≥ 2).

So there is a single threshold *k*\* = *QD*/(*aL*). The bucket brigade at *N* = 64 has *Q* = 256, *D* = 13, *a* = 6, *L* = 125, so *k*\* = 4.4: sharing wins for 2–4 instances and loses from 5 on, exactly as the CSV shows. At *N* = 8 *k*\* = 1.15, so sharing never wins; for the table lookup *Q* = *a* and *L* ≫ *D*, so *k*\* < 1 and sharing never wins in qubit·depth. The checker verifies the rule on all 64 rows.

**Zalka's bound (literature context, from memory).** Zalka (1999) shows Grover's √*N* query count is optimal and that *k* processors searching the same problem in parallel need Θ(√(*N*/*k*)) queries each, so the total number of oracle queries is √(*Nk*), more than one instance uses. Therefore a shared oracle serving *k* copies of the *same* search is strictly worse than running one instance, and sharing only makes sense for *k* *independent* searches, which is the case modeled above. This should be stated in the Step 5 write-up. (Riesebos and de Wolf's notes also cover this; neither copy is in `proofs/sources/`.)

## 4. Step 6: what the frame does under noise

**Lemma 4 (the frame is invisible to dephasing, not to amplitude damping).** Let 𝒩 be a one-qubit channel acting on an idle physical qubit that carries flag *P* ∈ {X, Y, Z}. The logical qubit then experiences the channel *P*𝒩*P*. For dephasing (Kraus operators √(1−p)·I, √p·Z) and for depolarizing, *P*𝒩*P* = 𝒩 for every Pauli *P*, because conjugating a Kraus operator by a Pauli changes it by at most a sign. For amplitude damping (K₀ = diag(1, √(1−γ)), K₁ = √γ |0⟩⟨1|), Z·AD·Z = AD but X·AD·X ≠ AD: the conjugated channel damps toward |1⟩. So an idle qubit that is logically |1⟩ but physically |0⟩ (an X flag) is a fixed point of AD and takes **no** idle error, while stored directly it decays with probability γ per tick; a logically |0⟩ qubit carrying an X flag is exposed instead.

*Proof.* Direct computation on the Kraus operators; `check_predictions.py` verifies all three channels and the fixed-point statement. ∎

What this says about the Step 6 results:
- Under **dephasing only** the compressed and original circuits differ only by the removed gates' own error (bound below). With one-qubit error 10⁻³ and the two-qubit gates untouched, that difference is below the benchmark's confidence intervals, which is what `noise_idle_sweep_dephasing.csv` shows: equal within CI for Grover at every idle rate.
- Under **amplitude damping** any difference is the fixed-point asymmetry of Lemma 4, not "fewer gates". It is physical (T₁ decay is amplitude damping), but it depends on which basis state the frame leaves on the physical line and can go either way. In Grover the diffusion's X⊗ⁿ leaves every line carrying an X flag during the MCZ, so the gain appears only with damping switched on.
- BV's small dephasing-only gain is the clock-tick effect: one layer fewer, so one idle tick fewer on every other qubit.

**Bound from gate count alone.** If each gate fails independently with probability ε and the circuit succeeds only when no gate fails, removing *g* of *G* gates multiplies the success probability by (1 − ε)^(−*g*) ≈ e^(*gε*). For Grover Clifford+T the pass removes 4*n* of 67*n* − 146 gates, about 6% of them for large *n*, so the predicted gain is a factor of about e^(0.06·*Gε*): small per gate, and only visible when *Gε* is of order 1. Depth matters separately through idle error, and the pass changes depth by a constant (4 layers for Grover, 1 for QFT), not by a multiple. This is a bound on what the noise simulation can show, not a prediction of its numbers: the simulated channels are not of this one-parameter form, and the checker does not test this section.

## 5. The working paper's claims, checked against the proofs

The working paper (`source/`) and the whitepaper (`whitepaper.html`) make five quantitative claims. Here is what the proofs and the closed forms say about each.

| Claim | Status | What is actually true |
|---|---|---|
| QFT: "1·*n* transforms removed" | **Overstated by 2×** | The bit-reversal network has ⌊*n*/2⌋ SWAPs, so ½ gate per qubit (3/2 CNOT-equivalents). Confirmed exactly for *n* = 4…64. Against Qiskit O3, which already relabels those SWAPs, the extra saving is 0. |
| Shor: "2·*n*, further multiplied by the error rate" | **Overstated; the error-rate factor has no basis in the gate model** | In the textbook *N* = 15 circuit the removable gates are the ⌊*n*/2⌋ inverse-QFT SWAPs (½ per counting qubit). The work register's X is flushed at a controlled-SWAP. Fewer gates do raise the success probability, but by at most about e^(*gε*) (§4), which is not a multiplier on the gate count. |
| Grover: "linear, constant depends on how the oracle and diffusion are compiled" | **Confirmed, and now quantified** | *n* + *z* gates per iteration with a native MCZ, 2*n* + 2*z* with the MCZ in Clifford+T (*z* = zero bits of the marked item). The paper was right that the constant depends on compilation. |
| "66% fewer qubits for 3 instances on a 1,024-line oracle" | **Confirmed as a ceiling, with a bookkeeping slip** | The saving is (1 − 1/*k*)·*Q*/(*a* + *Q*) ≤ 1 − 1/*k* = 2/3. The whitepaper's table (61/64/65/66% for *N* = 128…1024) uses 3*N* conventional versus *N* + 3⌈log₂*N*⌉ shared, counting the address registers on one side only; counting them on both gives 63.2/64.9/65.8/66.2%. With a real bucket brigade (*Q* ≈ 4*N*) the figure is 66.5%; with a table lookup (*Q* ≈ log *N*) it is 33%. Sharing also multiplies depth by about *k*, and wins in qubit·depth only while *k* < *QD*/(*aL*) (§3). |
| Proof [1]: "*UX* is unitary, so X need not be a separate step" | **True but incomplete** | Unitarity shows the merged gate exists, not that it is free. Lemma 1 of `formal-model.md` supplies the missing step: free when *U* is Clifford (frame update) or when *X U X* is a native gate of equal cost (substitution). Theorem 3 gives the exact class. |

Two further points about the whitepaper's figures. Figure 5 (fidelity against operations per tick, 1% per tick) illustrates the clock argument, but the pass changes depth by a constant (1 layer for QFT, 2–4 for Grover), so its benefit is in gate count, not ticks. Figure 7 (*k* instances, 128-line oracle) follows the same formula as the table and has the same slip; the corrected curve is (1 − 1/*k*)·128/135, which approaches 94.8% of 1 − 1/*k*.

## 6. Connectivity: when removing a SWAP costs CNOTs

`formal-model.md` §8 caveat 1 says SWAP removal assumes all-to-all coupling. The line-coupling benchmark (`gate_savings_line.csv`) makes this precise.

**With an optimal router the pass never hurts.** The original circuit, SWAPs included, is one valid routing of the relabeled circuit: a router may reinsert exactly the SWAPs the pass removed. So the minimum routed cost of the pass's output is at most the routed cost of the original, and `pass → route` is never worse than `route` alone.

**With a heuristic router it can.** SABRE (inside Qiskit O3) is not optimal, and the pass changes the circuit it sees. On the line, the structured families (BV, GHZ, syndrome, QFT, Shor-15, both Grovers, Qiskit's random Cliffords) are never worse in CNOT-equivalents in any of 216 rows, but random layered circuits are worse in 24 of 48 rows, better in 17 and equal in 7, with a worst case of +402 CNOT-equivalents on a circuit whose routed cost is about 11,800. Total gate count is still lower in 263 of the 264 rows, because the removed one-qubit gates outnumber the added routing CNOTs; the one exception (random Clifford, *n* = 4, seed 1: 31 gates to 33 with the same CNOT count) comes from O3's one-qubit resynthesis landing differently on the relabeled circuit. So on constrained hardware the SWAP savings can be partly or wholly undone by routing, and a fair benchmark must route *after* the pass, as the `pass+o3` column does. The `o3+pass` column on a coupling map deletes the router's SWAPs by relabeling and does not describe a hardware circuit.

## 7. Sources for this document

- The circuit definitions: `src/qarch/circuits.py`, `src/qarch/qram.py` on branch `worktree-code-compiler` (commit 8687387).
- The flag rules: `proofs/formal-model.md` §4–§7.
- From memory, not verified: Babbush et al. 2018 (*N* − 1 Toffoli, 4(*N* − 1) T for QROM), Low–Kliuchnikov–Schaeffer 2018 (select-swap count), Zalka 1999 (parallel Grover bound). The downloads were blocked in this session; see `TODO.md`.
