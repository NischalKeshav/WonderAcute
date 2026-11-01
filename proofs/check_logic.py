"""Checks the chain of reasoning in formal-model.md, step by step.

check_proofs.py checks the *identities* the proofs quote. This file checks that the
*steps* of each proof lead to its conclusion:

  L2    Lemma 2 with a general address table pi on n qubits (check_proofs.py only
        tests one SWAP on two qubits), including the SWAP(i,j) R_pi = R_pi' claim.
  INV   The induction step of Theorems 1 and 2: the invariant
            |psi_ideal>  =  P(x,z) R_pi |phi>   (up to global phase)
        holds after EVERY instruction of random circuits, in EVERY measurement
        branch. Instructions: H S Sdg X Y Z T Tdg Rz(theta) CNOT CZ SWAP, measure,
        reset, and any gate conditioned on an earlier outcome.
  MEAS  The measurement rule: matched branch probabilities and matched
        post-measurement states (that is what lets the induction continue).
  COMP  The compile-time claim after Theorem 1: with no mid-circuit measurement,
        run the stripped circuit once and fix the bits with a mask and a permutation.
  T3    Theorem 3 by exhaustive enumeration for n = 2, 3 (48 of 40,320 for n = 3).
  BND   Every row of the section-7 boundary table, including CP(theta).
  QFT   The QFT claim: drop the floor(n/2) final SWAPs, relabel, same unitary.

Run:  python3 proofs/check_logic.py
"""

import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_proofs import (CNOT, CZ, H, I2, S, SDG, SWAP, T, TDG, X, Y, Z, apply,  # noqa: E402
                          as_signed_pauli, close, controlled, equal_up_to_phase, kron,
                          measure_branches, pauli)

rng = np.random.default_rng(2026)


def haar(d):
    q, r = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return q * (np.diag(r) / abs(np.diag(r)))


def random_state(n):
    v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    return v / np.linalg.norm(v)


def rz(th):
    return np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])


def phase_gate(th):
    return np.diag([1, np.exp(1j * th)])


def cp(th):
    return np.diag([1, 1, 1, np.exp(1j * th)])


def relabel(state, n, pi):
    """R_pi |phi>: logical position i receives the content of physical line pi[i]."""
    return np.transpose(state.reshape([2] * n), axes=list(pi)).reshape(-1)


def frame_op(x, z):
    return kron(*(np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b) for a, b in zip(x, z)))


# ---------------------------------------------------------------- L2

def check_lemma2():
    for _ in range(200):
        n = int(rng.integers(2, 5))
        pi = list(rng.permutation(n))
        phi = random_state(n)
        k = int(rng.integers(1, n + 1))
        qs = [int(q) for q in rng.choice(n, k, replace=False)]
        u = haar(2 ** k)
        lhs = apply(relabel(phi, n, pi), n, u, qs)                  # U R_pi |phi>
        rhs = relabel(apply(phi, n, u, [pi[q] for q in qs]), n, pi)  # R_pi U^pi |phi>
        assert close(lhs, rhs)
        i, j = [int(q) for q in rng.choice(n, 2, replace=False)]
        pi2 = list(pi)
        pi2[i], pi2[j] = pi2[j], pi2[i]
        assert close(apply(relabel(phi, n, pi), n, SWAP, [i, j]), relabel(phi, n, pi2))
    print("[ok] L2   Lemma 2: U R_pi = R_pi U^pi for random pi, random k-qubit U (k=1..n), "
          "and SWAP(i,j) R_pi = R_pi' on 2-4 qubits")


# ---------------------------------------------------------------- INV / MEAS

GATES_1Q = {"H": H, "S": S, "Sdg": SDG, "X": X, "Y": Y, "Z": Z, "T": T, "Tdg": TDG}
GATES_2Q = {"CNOT": CNOT, "CZ": CZ, "SWAP": SWAP}
PAULI_BITS = {"X": (1, 0), "Y": (1, 1), "Z": (0, 1)}
CLIFFORD = {"H", "S", "Sdg", "CNOT", "CZ"}


def bit_rule(name, flags):
    """Section 5 bit form. flags = [(x,z), ...] in gate-argument order."""
    f = [list(v) for v in flags]
    if name == "H":
        f[0] = [f[0][1], f[0][0]]
    elif name in ("S", "Sdg"):
        f[0][1] ^= f[0][0]
    elif name == "CNOT":
        f[1][0] ^= f[0][0]
        f[0][1] ^= f[1][1]
    elif name == "CZ":
        xa, xb = f[0][0], f[1][0]
        f[0][1] ^= xb
        f[1][1] ^= xa
    return [tuple(v) for v in f]


def gate_matrix(op):
    name = op[0]
    if name == "Rz":
        return rz(op[2]), [op[1]]
    if name in GATES_1Q:
        return GATES_1Q[name], [op[1]]
    return GATES_2Q[name], [op[1], op[2]]


class Stats:
    branches = 0
    steps = 0
    free_gates = 0
    phys_gates = 0


def run_side_by_side(ops, n, stats):
    """Step the original circuit (psi) and the controller (phi, pi, x, z) together.

    Asserts the invariant after every instruction in every branch, and that each
    measurement/reset branch of the controller pairs with an ideal branch of the same
    probability. Returns the distribution of reported outcomes.
    """
    dist = {}

    def invariant(psi, phi, pi, x, z):
        assert equal_up_to_phase(psi, frame_op(x, z) @ relabel(phi, n, pi))
        stats.steps += 1

    def go(i, psi, phi, pi, x, z, bits, pr):
        invariant(psi, phi, pi, x, z)
        if i == len(ops):
            dist[tuple(bits)] = dist.get(tuple(bits), 0) + pr
            stats.branches += 1
            return
        op = ops[i]
        if op[0] in ("measure", "reset"):
            q = op[1]
            ideal = {m: (p, s) for m, p, s in measure_branches(psi, n, q)}
            for b, p, phi2 in measure_branches(phi, n, pi[q]):
                m = b ^ x[q]
                assert m in ideal and abs(ideal[m][0] - p) < 1e-9, "MEAS: branch probabilities differ"
                psi2 = ideal[m][1]
                x2, z2 = list(x), list(z)
                z2[q] = 0                                  # Z flag on a measured qubit is a phase
                if op[0] == "reset":                       # reset = measure, then flip to |0>
                    if m:
                        psi2 = apply(psi2, n, X, [q])
                    if b:
                        phi2 = apply(phi2, n, X, [pi[q]])
                    x2[q] = 0
                    go(i + 1, psi2, phi2, pi, x2, z2, bits, pr * p)
                else:
                    go(i + 1, psi2, phi2, pi, x2, z2, bits + [m], pr * p)
            return
        if op[0] == "if":                                  # ("if", bit_index, gate_op)
            if bits[op[1]] == 0:
                return go(i + 1, psi, phi, pi, x, z, bits, pr)
            op = op[2]
        name = op[0]
        u, qs = gate_matrix(op)
        psi2 = apply(psi, n, u, qs)                        # the original circuit applies the gate
        x, z, pi = list(x), list(z), list(pi)
        if name in PAULI_BITS:                             # controller: flip flags only
            x[qs[0]] ^= PAULI_BITS[name][0]
            z[qs[0]] ^= PAULI_BITS[name][1]
            stats.free_gates += 1
            return go(i + 1, psi2, phi, pi, x, z, bits, pr)
        if name == "SWAP":                                 # controller: edit the table only
            a, b = qs
            pi[a], pi[b] = pi[b], pi[a]
            x[a], x[b] = x[b], x[a]
            z[a], z[b] = z[b], z[a]
            stats.free_gates += 1
            return go(i + 1, psi2, phi, pi, x, z, bits, pr)
        lines = [pi[q] for q in qs]
        stats.phys_gates += 1
        if name in ("T", "Tdg", "Rz"):                     # substitution rule, flags unchanged
            if name == "Rz":
                u2 = rz(-op[2]) if x[qs[0]] else u
            else:
                u2 = {"T": TDG, "Tdg": T}[name] if x[qs[0]] else u
            return go(i + 1, psi2, apply(phi, n, u2, lines), pi, x, z, bits, pr)
        assert name in CLIFFORD
        new = bit_rule(name, [(x[q], z[q]) for q in qs])   # frame update by the section-5 table
        for q, (xq, zq) in zip(qs, new):
            x[q], z[q] = xq, zq
        go(i + 1, psi2, apply(phi, n, u, lines), pi, x, z, bits, pr)

    go(0, np.eye(2 ** n, dtype=complex)[0], np.eye(2 ** n, dtype=complex)[0],
       list(range(n)), [0] * n, [0] * n, [], 1.0)
    assert abs(sum(dist.values()) - 1) < 1e-9
    return dist


def random_op(n, ones, n_bits):
    r = rng.random()
    if r < 0.5 or n < 2:
        g = str(rng.choice(ones))
        q = int(rng.integers(n))
        return ("Rz", q, float(rng.uniform(0, 2 * np.pi))) if g == "Rz" else (g, q)
    a, b = (int(v) for v in rng.choice(n, 2, replace=False))
    return (str(rng.choice(["CNOT", "CZ", "SWAP"])), a, b)


def random_circuit(n, depth, use_t, mid):
    ones = ["H", "S", "Sdg", "X", "Y", "Z"] + (["T", "Tdg", "Rz"] if use_t else [])
    ops, n_bits = [], 0
    for _ in range(depth):
        r = rng.random()
        if mid and r < 0.12:
            ops.append(("measure", int(rng.integers(n))))
            n_bits += 1
        elif mid and r < 0.20:
            ops.append(("reset", int(rng.integers(n))))
        elif mid and r < 0.32 and n_bits:
            ops.append(("if", int(rng.integers(n_bits)), random_op(n, ones, n_bits)))
        else:
            ops.append(random_op(n, ones, n_bits))
    for q in range(n):
        ops.append(("measure", q))
    return ops


def check_invariant():
    cases = [("Clifford, final measurement only", False, False),
             ("Clifford + measure/reset/feed-forward", False, True),
             ("Clifford+T+Rz, final measurement only", True, False),
             ("Clifford+T+Rz + measure/reset/feed-forward", True, True)]
    for label, use_t, mid in cases:
        stats = Stats()
        for _ in range(120):
            n = int(rng.integers(1, 5))
            run_side_by_side(random_circuit(n, int(rng.integers(5, 25)), use_t, mid), n, stats)
        print(f"[ok] INV  {label}: invariant held at {stats.steps} steps over {stats.branches} "
              f"branches; {stats.free_gates} Pauli/SWAP handled with no physical gate")
    print("[ok] MEAS measurement and reset branches paired with equal probabilities every time")


# ---------------------------------------------------------------- COMP

def check_compile_time():
    """No mid-circuit measurement: strip Paulis/SWAPs, run once, fix bits classically."""
    for _ in range(150):
        n = int(rng.integers(1, 5))
        ops = random_circuit(n, int(rng.integers(5, 25)), True, False)
        gates = [op for op in ops if op[0] != "measure"]
        # Compile: a single pass over the gates with the controller rules (no state).
        pi, x, z, stripped = list(range(n)), [0] * n, [0] * n, []
        for op in gates:
            name = op[0]
            _, qs = gate_matrix(op)
            if name in PAULI_BITS:
                x[qs[0]] ^= PAULI_BITS[name][0]
                z[qs[0]] ^= PAULI_BITS[name][1]
            elif name == "SWAP":
                a, b = qs
                pi[a], pi[b] = pi[b], pi[a]
                x[a], x[b] = x[b], x[a]
                z[a], z[b] = z[b], z[a]
            elif name in ("T", "Tdg", "Rz"):
                lines = [pi[q] for q in qs]
                if name == "Rz":
                    stripped.append((rz(-op[2]) if x[qs[0]] else rz(op[2]), lines))
                else:
                    m = {"T": T, "Tdg": TDG}[name]
                    stripped.append((m.conj().T if x[qs[0]] else m, lines))
            else:
                lines = [pi[q] for q in qs]
                stripped.append((gate_matrix(op)[0], lines))
                for q, (xq, zq) in zip(qs, bit_rule(name, [(x[q], z[q]) for q in qs])):
                    x[q], z[q] = xq, zq
        # Run original and stripped circuits once each.
        psi = np.eye(2 ** n, dtype=complex)[0]
        for op in gates:
            u, qs = gate_matrix(op)
            psi = apply(psi, n, u, qs)
        phi = np.eye(2 ** n, dtype=complex)[0]
        for u, lines in stripped:
            phi = apply(phi, n, u, lines)
        p_ideal = np.abs(psi) ** 2
        p_phys = np.abs(phi) ** 2
        # Post-process: logical bit i = physical bit on line pi[i], XOR x[i].
        fixed = np.zeros_like(p_ideal)
        for s in range(2 ** n):
            bits = [(s >> (n - 1 - line)) & 1 for line in range(n)]
            logical = [bits[pi[i]] ^ x[i] for i in range(n)]
            fixed[int("".join(map(str, logical)), 2)] += p_phys[s]
        assert close(fixed, p_ideal)
        assert len(stripped) == sum(op[0] not in ("X", "Y", "Z", "SWAP") for op in gates)
    print("[ok] COMP compile-time form: stripped circuit + XOR mask + bit permutation reproduces "
          "the ideal distribution (150 random Clifford+T+Rz circuits)")


# ---------------------------------------------------------------- T3

def check_theorem3():
    for n in (2, 3):
        N = 2 ** n
        gens = []
        for i in range(n):
            gens.append(tuple(y ^ (1 << i) for y in range(N)))                      # X_i
        for i, j in itertools.combinations(range(n), 2):
            def sw(y, i=i, j=j):
                bi, bj = (y >> i) & 1, (y >> j) & 1
                return y ^ ((bi ^ bj) << i) ^ ((bi ^ bj) << j)
            gens.append(tuple(sw(y) for y in range(N)))                             # SWAP_ij
        ident = tuple(range(N))
        seen, frontier = {ident}, [ident]
        while frontier:                                                             # closure
            nxt = []
            for f in frontier:
                for g in gens:
                    h = tuple(g[f[y]] for y in range(N))
                    if h not in seen:
                        seen.add(h)
                        nxt.append(h)
            frontier = nxt
        assert len(seen) == 2 ** n * math.factorial(n)
        for f in seen:
            assert sorted(f) == list(range(N))                                      # bijection
            a = f[0]
            lin = [f[y] ^ a for y in range(N)]                                      # strip the flip
            # each output bit depends on exactly one input bit, and distinct bits -> permutation
            src = []
            for out in range(n):
                deps = [k for k in range(n) if (lin[1 << k] >> out) & 1]
                assert len(deps) == 1
                src.append(deps[0])
            assert sorted(src) == list(range(n))
            assert all(lin[y] == sum(((y >> src[o]) & 1) << o for o in range(n)) for y in range(N))
        cnot = tuple(y ^ (((y >> 1) & 1) << 0) for y in range(N))                  # CNOT 1 -> 0
        assert cnot not in seen
        if n == 3:
            tof = tuple(y ^ ((((y >> 2) & (y >> 1)) & 1) << 0) for y in range(N))
            assert tof not in seen
            assert len(seen) == 48 and math.factorial(N) == 40320
        print(f"[ok] T3   n={n}: X_i and SWAP_ij generate exactly {len(seen)} = 2^n n! classical maps, "
              f"all of the form y -> sigma(y) xor a; CNOT{' and Toffoli' if n == 3 else ''} not among them "
              f"(of {math.factorial(N)} reversible maps)")


# ---------------------------------------------------------------- BND

def check_boundary():
    tof = controlled(X, 2)
    ccz = controlled(Z, 2)
    dag = lambda m: m.conj().T  # noqa: E731
    neg = lambda u, k: pauli("X" + "I" * (k - 1)) @ u @ pauli("X" + "I" * (k - 1))  # noqa: E731

    # Clifford: residue always Pauli (54 rows already checked; spot check with random Clifford words)
    for _ in range(50):
        n = 3
        word = np.eye(8, dtype=complex)
        for _ in range(int(rng.integers(1, 8))):
            g = str(rng.choice(["H", "S", "CNOT", "CZ", "SWAP"]))
            if g in ("H", "S"):
                word = apply_mat(word, n, GATES_1Q[g], [int(rng.integers(3))])
            else:
                a, b = (int(v) for v in rng.choice(3, 2, replace=False))
                word = apply_mat(word, n, GATES_2Q[g], [a, b])
        p = "".join(rng.choice(list("IXYZ"), 3))
        if set(p) != {"I"}:
            assert as_signed_pauli(word @ pauli(p) @ dag(word), 3) is not None
    # T / Rz with X or Y: substitution gives T^dag / Rz(-theta); frame form residue is a Z rotation
    for th in rng.uniform(0, 2 * np.pi, 5):
        for p in (X, Y):
            assert close(p @ rz(th) @ p, rz(-th))
            assert equal_up_to_phase(rz(th) @ p @ dag(rz(th)), p @ rz(-2 * th))
        assert close(Z @ rz(th) @ Z, rz(th))
    assert equal_up_to_phase(T @ X @ dag(T), X @ SDG)               # theta = pi/4: residue is S^dag
    assert equal_up_to_phase(Y @ T @ Y, TDG)
    # Toffoli rows
    assert close(tof @ pauli("IIX"), pauli("IIX") @ tof)
    assert close(tof @ pauli("ZII"), pauli("ZII") @ tof)
    assert close(tof @ pauli("IZI"), pauli("IZI") @ tof)
    assert close(tof @ pauli("XII") @ tof, pauli("XII") @ kron(I2, CNOT))     # frame form
    assert close(neg(tof, 3), tof @ kron(I2, CNOT))                           # = negative-control Tof
    assert close(pauli("IIZ") @ tof @ pauli("IIZ"), tof @ kron(CZ, I2))       # substitution form
    assert close(tof @ pauli("IIZ") @ tof, pauli("IIZ") @ kron(CZ, I2))       # frame form
    assert as_signed_pauli(tof @ pauli("IIZ") @ tof, 3) is None
    # CCZ rows
    assert close(ccz @ pauli("XII") @ dag(ccz), pauli("XII") @ kron(I2, CZ))
    assert close(neg(ccz, 3), ccz @ kron(I2, CZ))
    for p in ("ZII", "IZI", "IIZ"):
        assert close(ccz @ pauli(p), pauli(p) @ ccz)
    # CP(theta) rows (new): X on either qubit -> CP(-theta) plus a phase gate on the other qubit
    for th in rng.uniform(0, 2 * np.pi, 5):
        assert close(pauli("XI") @ cp(th) @ pauli("XI"), cp(-th) @ kron(I2, phase_gate(th)))
        assert close(pauli("IX") @ cp(th) @ pauli("IX"), cp(-th) @ kron(phase_gate(th), I2))
        assert close(pauli("ZI") @ cp(th) @ pauli("ZI"), cp(th))
        assert as_signed_pauli(cp(th) @ pauli("XI") @ dag(cp(th)), 2) is None
        # and CP decomposes into Rz + CNOT, where the flags pass gate by gate
        dec = kron(phase_gate(th / 2), phase_gate(th / 2)) @ CNOT @ kron(I2, phase_gate(-th / 2)) @ CNOT
        assert close(dec, cp(th))
        assert equal_up_to_phase(phase_gate(th), rz(th))
    print("[ok] BND  every row of the section-7 table holds: Clifford residues are Pauli; T/Rz residue "
          "is X Rz(-2t) (S^dag for T) or substitute Rz(-t); Toffoli X-on-control -> CNOT residue / "
          "negative control; Toffoli Z-on-target -> CZ; CCZ X-on-control -> CZ; "
          "CP(t) X-on-either -> CP(-t) + P(t) on the other qubit; CP = Rz+CNOT decomposition")


def apply_mat(m, n, u, qs):
    """Left-multiply matrix m by u acting on qubits qs."""
    return np.column_stack([apply(m[:, c], n, u, qs) for c in range(m.shape[1])])


# ---------------------------------------------------------------- QFT

def check_qft():
    for n in range(2, 6):
        full = np.eye(2 ** n, dtype=complex)
        noswap = np.eye(2 ** n, dtype=complex)
        n_cp = 0
        for j in range(n):
            full = apply_mat(full, n, H, [j])
            noswap = apply_mat(noswap, n, H, [j])
            for k in range(j + 1, n):
                th = np.pi / 2 ** (k - j)
                full = apply_mat(full, n, cp(th), [k, j])
                noswap = apply_mat(noswap, n, cp(th), [k, j])
                n_cp += 1
        n_swaps = 0
        for j in range(n // 2):
            full = apply_mat(full, n, SWAP, [j, n - 1 - j])
            n_swaps += 1
        assert n_swaps == n // 2
        # Reference: the QFT matrix.
        w = np.exp(2j * np.pi / 2 ** n)
        qft = np.array([[w ** (a * b) for b in range(2 ** n)] for a in range(2 ** n)]) / np.sqrt(2 ** n)
        assert close(full, qft)
        # Address table after absorbing the SWAPs: pi = bit reversal. R_pi applied to every column.
        pi = list(range(n))[::-1]
        relabeled = np.column_stack([relabel(noswap[:, c], n, pi) for c in range(2 ** n)])
        assert close(relabeled, qft)
    print("[ok] QFT  n=2..5: QFT with its floor(n/2) final SWAPs removed, relabeled by bit reversal, "
          "equals the QFT matrix; the CP gates carry no Paulis")


if __name__ == "__main__":
    check_lemma2()
    check_invariant()
    check_compile_time()
    check_theorem3()
    check_boundary()
    check_qft()
    print("all logic checks passed")
