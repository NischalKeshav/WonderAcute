# Step 1: Formal Model, Flag Rules, and Proofs

This is research plan Step 1, "Formal model" (`research-plan.md` §3). It defines the address table and flags, proves that X/Y/Z flags and address swaps reproduce the original circuit, derives how flags move through H, S, S†, CNOT, CZ and SWAP, and works out what happens at T and other non-Clifford gates.

Two scripts check this document:

- `proofs/check_proofs.py` checks every identity quoted below against the matrices, and regenerates `proofs/propagation-tables.json` (the flag tables, in a form the Step 2 compiler's unit tests can check row by row).
- `proofs/check_logic.py` checks the chain of reasoning: Lemma 2 with a general address table, the invariant of §2 after every instruction of random circuits in every measurement branch (the induction step of Theorems 1–2, including reset and feed-forward), the compile-time form, Theorem 3 by exhaustive enumeration, every row of the boundary table in §7, and the QFT claim.

Run both with `python3 proofs/check_proofs.py && python3 proofs/check_logic.py`.

---

## 1. Notation

- *n* logical qubits are labeled 1…*n*. Each is stored on a physical line.
- Pauli matrices: *I, X, Y, Z*, with *Y = iXZ*. For bit vectors *a, b* ∈ {0,1}ⁿ, write *P(a,b)* = ⊗ᵢ Xᵃⁱ Zᵇⁱ. Every *n*-qubit Pauli operator equals *P(a,b)* up to a phase in {±1, ±i}.
- ≅ means "equal up to a global phase". A global phase cannot be observed, so the proofs drop it freely.
- A **Clifford** gate is a unitary *C* that maps every Pauli to a Pauli under conjugation: *CPC*† ≅ Pauli for every Pauli *P*. H, S, S†, CNOT, CZ and SWAP are Clifford. T is not.
- Circuits are read left to right. "Apply *A*, then *B*" is the operator *BA*.

## 2. The architecture as a mathematical object

The working paper's classical controller holds two kinds of classical data:

1. **Address table** π: a permutation of {1…*n*}. Logical qubit *i* lives on physical line π(*i*).
2. **Flags** (*x*, *z*) ∈ {0,1}ⁿ × {0,1}ⁿ: one X-flag bit and one Z-flag bit per logical qubit. Both flags set means a Y flag.

Let *R*_π be the unitary that moves line π(*i*) into position *i* for every *i*. If |φ⟩ is the state of the physical lines, the controller's **claim** is

> **(Invariant)** |ψ_ideal⟩ ≅ *P(x,z)* · *R*_π · |φ⟩,

where |ψ_ideal⟩ is the state the original circuit would have produced at this point. At the start π = identity, *x = z* = 0, |φ⟩ = |ψ_ideal⟩ = |0…0⟩, so the invariant holds.

Every rule below keeps the invariant true. Theorem 1 follows by induction over the gates.

## 3. Three lemmas

**Lemma 1 (pushing a flag through any gate).** For every unitary *U* and every Pauli *P*,

  *UP* = (*UPU*†) *U*  and also  *UP* = *P* (*P*†*UP*).

*Proof.* Insert *U*†*U* = *I* or *PP*† = *I*. ∎

These two identities give the two ways to deal with a flag in front of a gate:
- **Frame update:** apply *U* unchanged and replace the flag *P* by *UPU*†. This keeps the flag a Pauli exactly when *UPU*† is a Pauli. That holds for every *P* if and only if *U* is Clifford (the definition).
- **Gate substitution:** keep the flag *P* and apply *P*†*UP* in place of *U*. This works for every *U*. The only question is whether *P*†*UP* costs as much as *U* does on the hardware.

**Lemma 2 (relabeling).** Let *U* act on logical qubits *i₁…i_k*, and let *U*^π be the same matrix applied to lines π(*i₁*)…π(*i_k*). Then *U R*_π = *R*_π *U*^π. In particular a logical SWAP(*i,j*) satisfies SWAPᵢⱼ *R*_π = *R*_π′, where π′ is π with the entries for *i* and *j* exchanged.

*Proof.* *R*_π only renames tensor factors, and *U* acts on the factors named *i₁…i_k*. After the renaming those factors are lines π(*i₁*)…π(*i_k*). For SWAP, exchanging two factors and then renaming by π is the same as renaming by π′. ∎

**Lemma 3 (measurement through a flag).** Measure every qubit of *P(a,b)*|ψ⟩ in the computational basis. Outcome *m* occurs with probability Pr_ψ(*m* ⊕ *a*). The Z flags *b* have no effect.

*Proof.* *Zᵇ*|*y*⟩ = (−1)^{*b·y*}|*y*⟩ and *Xᵃ*|*y*⟩ = |*y*⊕*a*⟩, so ⟨*m*|*XᵃZᵇ*|ψ⟩ = ±⟨*m*⊕*a*|ψ⟩. Squaring removes the sign. The same argument applies to measuring a single qubit *i*: the outcome bit is flipped exactly when *aᵢ* = 1. ∎

## 4. The controller's rules

| Instruction in the original circuit | What the controller does | Physical operation |
|---|---|---|
| **X, Y or Z** on *i* | *xᵢ* ^= 1 (X), *zᵢ* ^= 1 (Z), or both (Y) | **none** |
| **SWAP**(*i,j*) | swap π(*i*) ↔ π(*j*), swap (*xᵢ,zᵢ*) ↔ (*xⱼ,zⱼ*) | **none** |
| Clifford *C* on *i₁…i_k* | apply *C* on lines π(*i₁*)…π(*i_k*), then update flags by the table in §5 | *C* (same gate) |
| **T** (or T†) on *i* | if *xᵢ* = 1, apply T† (or T) instead; flags unchanged | one T-type gate |
| Rz(θ) on *i* | if *xᵢ* = 1, apply Rz(−θ) instead; flags unchanged | one Rz |
| Measure *i* (Z basis) | measure line π(*i*) and get *b*; **report *m* = *b* ⊕ *xᵢ***; set *zᵢ* = 0 | measurement |
| Reset *i* to \|0⟩ | reset line π(*i*); set *xᵢ = zᵢ* = 0 | reset |
| Gate conditioned on an earlier outcome *m* | use the reported (corrected) *m* as the condition | as above |
| End of circuit, all qubits measured | nothing more | none |
| End of circuit, state handed on unmeasured | either apply Xˣⁱ Zᶻⁱ on line π(*i*) for every *i* (at most *n* single-qubit gates, depth 1) and hand π on, or hand (π, *x*, *z*) on with the state. The relabeling π is never undone physically; the next routine reads it. | at most *n* 1-qubit gates |

Why each rule keeps the invariant:

- *Pauli Q.* The new ideal state is *Q P(x,z) R*_π|φ⟩ ≅ *P(x′,z′) R*_π|φ⟩, because the product of two Paulis is a Pauli up to phase, with bits added mod 2.
- *SWAP.* SWAPᵢⱼ *P(x,z)* = *P(x′,z′)* SWAPᵢⱼ, where the flags of *i* and *j* are exchanged. Then Lemma 2 gives SWAPᵢⱼ *R*_π = *R*_π′.
- *Clifford C.* *C P R*_π|φ⟩ = (*CPC*†) *C R*_π|φ⟩ = *P′ R*_π *C*^π|φ⟩, by Lemma 1 (frame update) and then Lemma 2.
- *T.* Lemma 1 (substitution) with T X = e^{iπ/4} X T†, T Y = e^{iπ/4} Y T† and T Z = Z T. Only the X part of the flag matters (§6).
- *Measurement.* Lemma 3 shows the reported bit has the right distribution. After the outcome is fixed, qubit *i* is in a basis state, so *Zᵢ* acts on it as ±1, a global phase, and *zᵢ* can be cleared. *xᵢ* stays set because the line holds |*b*⟩ while the ideal qubit holds |*m*⟩ = Xˣⁱ|*b*⟩.
- *Reset.* Reset to |0⟩ is "measure, then apply X if the result is 1". Under the rules above that is: measure line π(*i*) (get *b*, reported *m* = *b* ⊕ *xᵢ*), then the conditional X becomes *xᵢ* ^= *m*, giving *xᵢ* = *b* with the line in |*b*⟩. Physically resetting the line to |0⟩ instead and setting *xᵢ* = *zᵢ* = 0 describes the same ideal state. Either way the flags on the other qubits are untouched, since the Paulis on them commute with anything done to qubit *i*. The invariant is kept branch by branch, which is all that outcome distributions depend on.

## 5. Flag propagation tables

**Bit form** (signs dropped, since a frame sign is a global phase). This is the form `frame.py` uses. Each rule takes the flags in front of the gate and gives the flags behind it.

| Gate | Update |
|---|---|
| H | swap *x* ↔ *z* |
| S, S† | *z* ^= *x* |
| CNOT(c, t) | *x_t* ^= *x_c*; *z_c* ^= *z_t* |
| CZ(a, b) | *z_a* ^= *x_b*; *z_b* ^= *x_a* |
| SWAP(a, b) | swap (*x_a, z_a*) ↔ (*x_b, z_b*), and swap π(a) ↔ π(b) |
| X, Y, Z | flags of a Pauli gate are flipped (§4), and existing flags pass through unchanged |

**Signed form** (exact matrices). `propagation-tables.json` lists *G P G*† = ±*Q* for every non-identity Pauli *P* on the gate's qubits: 3 rows for each 1-qubit gate and 15 rows for each 2-qubit gate, 54 rows in total. The letters follow gate-argument order, with control first. The generators:

| Gate | X ↦ | Y ↦ | Z ↦ |
|---|---|---|---|
| H | Z | −Y | X |
| S | Y | −X | Z |
| S† | −Y | X | Z |

| Gate | XI ↦ | IX ↦ | ZI ↦ | IZ ↦ |
|---|---|---|---|---|
| CNOT(c,t) | XX | IX | ZI | ZZ |
| CZ | XZ | ZX | ZI | IZ |
| SWAP | IX | XI | IZ | ZI |

All other rows are products of these. Both forms are checked against the matrices by the script.

## 6. Theorems

**Theorem 1 (Clifford circuits).** Take any circuit made of Clifford gates, Pauli gates, SWAPs, computational-basis measurements (at the end or in the middle), resets, and gates conditioned on measurement outcomes. Run it under the controller rules of §4. Then:

(a) no physical X, Y, Z or SWAP gate is ever applied;
(b) every other gate is applied exactly once, unchanged, on relabeled lines;
(c) the joint distribution of all reported measurement outcomes equals that of the original circuit.

*Proof.* Use induction over the instructions. The invariant of §2 holds at the start and is preserved by every rule (§4). So at each measurement the ideal state is *P(x,z)* *R*_π|φ⟩, and Lemma 3 (single-qubit version) shows the reported bit has the correct conditional distribution given all earlier reported bits. Conditioned gates use the reported bits, so every branch follows the original circuit's branch. Multiplying the conditional probabilities gives (c). Parts (a) and (b) can be read off the table in §4. ∎

When the circuit has no mid-circuit measurements, the flags are fixed at compile time. The compiler can then delete the Paulis and SWAPs, relabel the remaining gates, and post-process the outcome bits with a fixed XOR mask and bit permutation. Nothing has to happen at run time.

**Theorem 2 (extension to Clifford+T, and to Rz(θ)).** Add T, T† and Rz(θ) gates to the circuits of Theorem 1, and handle them with the substitution rule (T ↔ T† and Rz(θ) ↔ Rz(−θ) when *xᵢ* = 1). Then (a)–(c) still hold, and the number of T-type gates is unchanged.

*Proof.* The rule preserves the invariant by Lemma 1 (substitution) and the identities

  T X = e^{iπ/4} X T†,  T Y = e^{iπ/4} Y T†,  T Z = Z T,  X Rz(θ) X = Y Rz(θ) Y = Rz(−θ),  Z Rz(θ) Z = Rz(θ).

The rest of the proof is the same as for Theorem 1. ∎

*Two versions of the T rule (both proved).*
1. *Substitution (used above):* the flag stays a Pauli and the gate becomes T†. This costs nothing extra if the controller can choose between T and T† when it issues the gate. With no mid-circuit measurements that choice is made at compile time. With feed-forward it becomes a run-time choice, which is the same kind of decision as an ordinary classically conditioned gate.
2. *Frame update (gate left unchanged):* T X T† = (X+Y)/√2 = e^{iπ/4} X S†, and in the other direction T† X T = (X−Y)/√2. Neither is a Pauli, so if T is left as is, the flag picks up an **S-type Clifford correction**. A table of bits can no longer represent the flag. Either the correction is applied as a physical S†, or the frame must be extended to track Cliffords. Applying the S† and fusing it with the T gives T S† = T†, which is version 1.

**Theorem 3 (exactly which classical reversible functions are free).** This answers the research question. The operations the controller can perform with no physical gate are exactly the group

  𝒢ₙ = { phase · *P(a,b)* · *R*_π }  (Paulis combined with qubit permutations; 4ⁿ·*n*! elements up to phase).

On classical inputs |*y*⟩ these act as *y* ↦ σ(*y*) ⊕ *a*, where σ permutes bit positions. So **a classical reversible function on *n* bits can be done with no quantum gates if and only if it is a bit permutation followed by a fixed bit flip.** There are 2ⁿ·*n*! such functions out of (2ⁿ)! reversible functions. For *n* = 3 that is 48 out of 40,320. CNOT and Toffoli are not in the class.

*Proof.* Table-only rules can only change π and (*x,z*), so every reachable operation is in 𝒢ₙ, and every element of 𝒢ₙ is reachable. On basis states, *Zᵇ* contributes only a sign and *Xᵃ R*_π maps |*y*⟩ ↦ |σ(*y*)⊕*a*⟩. Each output bit of such a map depends on exactly one input bit. CNOT's target output *y_c ⊕ y_t* depends on two input bits, so CNOT is not in the class. There is also a quantum reason: every element of 𝒢ₙ is a product of single-qubit unitaries followed by a permutation, so it maps product states to product states, while CNOT(|+⟩|0⟩) is a Bell state. ∎

## 7. Where the method stops working

Lemma 1 always allows a substitution, so the real question is whether *P*†*UP* costs as much as *U*.

| Gate *U* | Flag on its qubits | Frame-update form *UPU*† | Substitution form *P*†*UP* | Free? |
|---|---|---|---|---|
| Clifford | any | Pauli | — | **yes** |
| T, Rz(θ) | X or Y | Pauli × Rz(−2θ), a Z-rotation residue (S† for T) | T†, Rz(−θ) | **yes, by substitution** |
| T, Rz(θ) | Z | commutes | — | **yes** |
| CP(θ) (controlled phase) | Z on either qubit | commutes | — | **yes** |
| CP(θ) | X on one qubit | Pauli × 2-qubit non-Clifford residue | CP(−θ) · P(θ) on the *other* qubit | **almost**: one extra single-qubit phase gate, which is free on hardware with virtual Z (McKay et al.) |
| Toffoli (c₁,c₂,t) | X on t, or Z on c₁/c₂ | commutes | — | **yes** |
| Toffoli | X on c₁ | X_{c₁} · CNOT(c₂→t), a Clifford residue | negative-control Toffoli | only if the negative-control variant is native |
| Toffoli | Z on t | Z_t · CZ(c₁,c₂) | Tof · CZ(c₁,c₂) | **no**: one extra CZ either way |
| CCZ / multi-controlled Z | X on a control | X · CZ (or C^{k−1}Z) on the others | negative-control MCZ | only if the variant is native |

Here P(θ) = diag(1, e^{iθ}) ≅ Rz(θ). The CP(θ) rows matter for QFT and Shor: CP(θ) = P(θ/2)_a P(θ/2)_b · CNOT · P(−θ/2)_b · CNOT, so when it is decomposed the flags pass through gate by gate with Rz angles negated (Theorem 2), and when it is kept native the substitution form above applies.

So:

- **On a Clifford+T machine**, where Toffoli and multi-controlled Z are broken down into Clifford+T, Theorems 1–2 apply gate by gate. Every Pauli and SWAP can be removed with no added gates, including the X gates around the multi-controlled Z in Grover's diffusion operator. The cost is that the controller must choose T versus T† when it issues each gate. Treating that decomposition as exact relies on the standard Clifford+T decompositions of Toffoli (Nielsen & Chuang §4.3). This document does not verify those decompositions, and Step 2 should test them on the actual Qiskit decomposition.
- **On a machine where Toffoli or CCZ is a native gate**, the boundary is real. An X flag on a control leaves a CNOT or CZ residue (or needs a negative-control variant), and a Z flag on a Toffoli target leaves a CZ. These are multi-qubit, non-Clifford native gates, which is where the "yes, for a specific class" boundary from the plan's expected outcomes lies.
- **QFT:** the SWAPs at the end of the circuit are absorbed by π (Lemma 2), and the controlled-phase gates contain no Paulis. The savings are therefore exactly the ⌊*n*/2⌋ final SWAPs. That matches the working paper's "linear" entry, but the slope is 1/2 SWAP per qubit, which is 3/2 CNOTs per qubit when each SWAP is counted as 3 CNOTs. In Shor's circuit the QFT follows modular arithmetic that contains X gates, so flags do enter the CP(θ) gates there; by the CP row above they pass through at the cost of single-qubit phase gates, which are free on virtual-Z hardware and otherwise count against the savings.

## 8. Caveats Steps 2–6 should carry forward

1. **Connectivity.** Lemma 2 assumes any two lines can interact. On hardware with limited coupling, a relabeled gate *U*^π can land on lines that are not connected. Removing a SWAP is then only valid if the router does not need it. Routing SWAPs (Li, Ding & Xie 2019) exist for exactly this reason, so the Step 3 baseline should state its coupling map.
2. **Paulis are cheap gates.** Removing X/Y/Z saves 1-qubit gates and depth, but no 2-qubit gates. The 2-qubit savings come only from SWAPs (3 CNOTs each). Step 3's "two-qubit gates" column will show this.
3. **Unmeasured outputs.** If a state is handed on instead of measured (for example, GHZ preparation feeding another routine), the frame has to be applied or handed on with the state (last row of the §4 table).
4. **Measurement basis.** Lemma 3 covers Z-basis measurement. Measuring in another basis means applying a Clifford first and then measuring in Z, which the rules already handle.
5. **Relation to the working paper's proof [1].** The paper argues that *UX* is unitary and can therefore be applied as a single gate. That is true, but it does not show the merged gate is free. Lemma 1 adds what was missing: the merge costs nothing when *U* is Clifford (frame update) or when *P*†*UP* is a native gate of equal cost (substitution).

## 9. Sources

The results above are standard in the Pauli-frame and stabilizer literature. This document restates them in the architecture's terms and proves them directly. The arXiv copies are saved in `proofs/sources/`.

**Checked against the arXiv text (2026-10-28):**

- Gottesman, D. (1998). The Heisenberg representation of quantum computers. arXiv:quant-ph/9807006. §3 "The Clifford Group" defines the Clifford group as the normalizer of the Pauli group. Table 1 (p. 4) gives exactly the H, S and CNOT rows of our §5: R: X→Z, Z→X; P: X→Y, Z→Z; CNOT: XI→XX, IX→IX, ZI→ZI, IZ→ZZ. Gottesman writes R for H and P for S. His Y is the standard σ_y, the same as ours; a footnote notes that earlier papers used Y = −iσ_y. p. 17 says that Clifford gates conditioned on measurement results keep the circuit efficiently trackable, which supports the feed-forward part of Theorem 1.
- Aaronson, S., & Gottesman, D. (2004). Improved simulation of stabilizer circuits. arXiv:quant-ph/0406196. §III "Efficient simulation of stabilizer circuits" gives the bit rules word for word as in §5 (and `frame.py`): CNOT a→b sets x_b ⊕= x_a and z_a ⊕= z_b; Hadamard swaps x_a and z_a; Phase sets z_a ⊕= x_a. The paper also updates a sign bit r, which the frame drops as a global phase. The journal reference (*Phys. Rev. A* 70, 052328) does not appear in the arXiv copy and is cited from memory.
- Knill, E. arXiv:quant-ph/0410199, "Quantum Computing with Very Noisy Devices" (v2, Nov 2004). This defines the "Pauli frame" as a Pauli product that the controller tracks instead of applying (p. 3). It also states that Pauli corrections "can be absorbed into the Pauli frame and do not need to be known immediately", while for non-Clifford gates "a non-Pauli compensation may be required". That is the frame-update side of the T boundary in §6–7. **Note:** the arXiv title differs from the Nature title in the plan's bibliography ("Quantum computing with realistically noisy devices", *Nature* 434, 39–44, 2005). I believe they are the same work, but the arXiv copy does not mention Nature, so this link is unconfirmed.
- McKay, D. C., Wood, C. J., Sheldon, S., Chow, J. M., & Gambetta, J. M. Efficient Z-gates for quantum computing. arXiv:1612.00858 (v2, 2017). A "virtual Z-gate" is implemented by adding a phase offset to all later drive pulses. The classical control hardware does this, and it takes zero time (p. 1). This is hardware precedent for Z rotations handled purely classically. It is close to Theorem 2's Rz rule but not the same mechanism: McKay changes the frame of later gates, while Theorem 2 negates the angle of the gate itself. The journal reference (*Phys. Rev. A* 96, 022330) does not appear in the arXiv copy and is cited from memory.
- Li, G., Ding, Y., & Xie, Y. Tackling the qubit mapping problem for NISQ-era quantum devices. arXiv:1809.02573 (v2, May 2019). SABRE inserts SWAPs because physical qubits have limited connections (abstract, Fig. 3), and "one SWAP even requires 3 CNOT gates" (p. 3). These support caveats 1 and 2 in §8. ASPLOS '19 is the venue given in the plan's bibliography.

**Not yet checked (no copy in `proofs/sources/`):**

- Riesebos, L., et al. (2017). Pauli frames for quantum computer architectures. *DAC '17*. Cited from memory in `research-directions.md`.
- Nielsen, M. A., & Chuang, I. L. (2010). §4.3 (Toffoli built from Clifford+T gates) and §10.5.2 (the Clifford group). Section numbers are from memory.
