"""Numerical checks for every identity used in proofs/formal-model.md.

Run:  python3 proofs/check_proofs.py            (check everything)
      python3 proofs/check_proofs.py --write    (regenerate propagation-tables.json, then check)

Conventions (same as formal-model.md):
  * A Pauli string such as "XZ" lists one letter per gate argument, in argument order.
  * The first argument is the most significant qubit in the matrix (kron(first, second)).
  * CNOT and CZ take (control, target).
  * A table row "C: P -> s Q" means  C P C^dagger = s Q  exactly, with s in {+1, -1}.
    This is the rule for moving a Pauli flag from before C to after C.
"""

import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TABLES = HERE / "propagation-tables.json"

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j])
SDG = S.conj().T
T = np.diag([1, np.exp(1j * np.pi / 4)])
TDG = T.conj().T
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def pauli(s):
    return kron(*(PAULI[c] for c in s))


def controlled(u, n_controls=1):
    """Apply u to the last qubit when all n_controls leading qubits are 1."""
    d = 2 ** (n_controls + 1)
    out = np.eye(d, dtype=complex)
    out[d - 2:, d - 2:] = u
    return out


CNOT = controlled(X)
CZ = controlled(Z)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)

CLIFFORD_GATES = {"H": H, "S": S, "Sdg": SDG, "CNOT": CNOT, "CZ": CZ, "SWAP": SWAP}


def close(a, b):
    return np.allclose(a, b, atol=1e-10)


def equal_up_to_phase(a, b):
    """True if a = c * b for some |c| = 1."""
    idx = np.unravel_index(np.argmax(np.abs(b)), b.shape)
    if abs(b[idx]) < 1e-12:
        return False
    c = a[idx] / b[idx]
    return abs(abs(c) - 1) < 1e-10 and close(a, c * b)


def as_signed_pauli(m, n):
    """Return (sign, string) with m = sign * pauli(string), or None if m is not a signed Pauli."""
    for s in itertools.product("IXYZ", repeat=n):
        p = pauli("".join(s))
        for sign in (1, -1):
            if close(m, sign * p):
                return sign, "".join(s)
    return None


# ---------------------------------------------------------------- tables

def build_tables():
    tables = {}
    for name, g in CLIFFORD_GATES.items():
        n = int(np.log2(g.shape[0]))
        rows = []
        for s in itertools.product("IXYZ", repeat=n):
            p = "".join(s)
            if set(p) == {"I"}:
                continue
            sign, q = as_signed_pauli(g @ pauli(p) @ g.conj().T, n)
            rows.append({"in": p, "sign": sign, "out": q})
        tables[name] = {"qubits": n, "rows": rows}
    return {
        "convention": (
            "Row {in: P, sign: s, out: Q} means G P G^dagger = s Q. Letters are in gate-argument "
            "order; CNOT and CZ are (control, target). Move a flag P from before G to after G by "
            "replacing it with Q. The sign is a global phase and can be dropped by the frame."
        ),
        "gates": tables,
    }


def check_tables(data):
    n_rows = 0
    for name, entry in data["gates"].items():
        g = CLIFFORD_GATES[name]
        n = entry["qubits"]
        assert len(entry["rows"]) == 4 ** n - 1, name
        for row in entry["rows"]:
            lhs = g @ pauli(row["in"]) @ g.conj().T
            assert close(lhs, row["sign"] * pauli(row["out"])), (name, row)
            n_rows += 1
    print(f"[ok] propagation tables: {len(data['gates'])} gates, {n_rows} rows match matrices")


# ---------------------------------------------------------------- lemmas

def check_lemmas():
    rng = np.random.default_rng(1)

    def haar(d):
        q, r = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
        return q * (np.diag(r) / abs(np.diag(r)))

    # Lemma 1: P U = U (U^dag P U) for every unitary U (pure algebra; sampled here).
    for _ in range(20):
        u = haar(4)
        for p in ("XI", "IZ", "YX"):
            assert close(pauli(p) @ u, u @ (u.conj().T @ pauli(p) @ u))
    print("[ok] Lemma 1: P U = U (U^dag P U) on random 2-qubit unitaries")

    # Lemma 2 (relabeling): SWAP (A x B) SWAP = B x A, so SWAP then U = U' then SWAP.
    for _ in range(20):
        a, b = haar(2), haar(2)
        assert close(SWAP @ kron(a, b) @ SWAP, kron(b, a))
        u = haar(4)
        assert close(u @ SWAP, SWAP @ (SWAP @ u @ SWAP))
    print("[ok] Lemma 2: SWAP relabels any gate (U SWAP = SWAP U', U' = SWAP U SWAP)")

    # T gate: conjugating a Pauli by T does not give a Pauli.
    txt = T @ X @ TDG
    assert close(txt, (X + Y) / np.sqrt(2))
    assert close(TDG @ X @ T, (X - Y) / np.sqrt(2))
    assert as_signed_pauli(txt, 1) is None
    # ...but it is a Pauli times a Clifford: T X T^dag = e^{i pi/4} X S^dag.
    assert close(txt, np.exp(1j * np.pi / 4) * X @ SDG)
    # Substitution form: T X = e^{i pi/4} X T^dag, T Y = e^{i pi/4} Y T^dag, T Z = Z T.
    assert close(T @ X, np.exp(1j * np.pi / 4) * X @ TDG)
    assert close(T @ Y, np.exp(1j * np.pi / 4) * Y @ TDG)
    assert close(T @ Z, Z @ T)
    assert close(TDG @ X, np.exp(-1j * np.pi / 4) * X @ T)
    print("[ok] T: T X T^dag = (X+Y)/sqrt2 = e^{i pi/4} X S^dag (not Pauli); "
          "T X = e^{i pi/4} X T^dag (substitution)")

    # Rz(theta): X Rz(theta) X = Rz(-theta).
    for th in rng.uniform(0, 2 * np.pi, 10):
        rz = np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
        rzm = np.diag([np.exp(1j * th / 2), np.exp(-1j * th / 2)])
        assert close(X @ rz @ X, rzm)
        assert close(Z @ rz @ Z, rz)
    print("[ok] Rz: X Rz(t) X = Rz(-t), Z Rz(t) Z = Rz(t)")

    # Toffoli (controls 0,1; target 2) and CCZ.
    tof = controlled(X, 2)
    ccz = controlled(Z, 2)
    neg0_tof = kron(X, I2, I2) @ tof @ kron(X, I2, I2)
    # X on a control: frame-update form produces a CNOT (Clifford, not Pauli).
    assert close(tof @ pauli("XII") @ tof.conj().T, pauli("XII") @ kron(I2, CNOT))
    assert as_signed_pauli(tof @ pauli("XII") @ tof.conj().T, 3) is None
    # Substitution form: Tof X_0 = X_0 (negative-control-0 Toffoli).
    assert close(tof @ pauli("XII"), pauli("XII") @ neg0_tof)
    # Commuting cases.
    assert close(tof @ pauli("IIX"), pauli("IIX") @ tof)
    assert close(tof @ pauli("ZII"), pauli("ZII") @ tof)
    # Z on the target: Z_t Tof Z_t = Tof CZ_{01}  (extra CZ needed).
    assert close(pauli("IIZ") @ tof @ pauli("IIZ"), tof @ kron(CZ, I2))
    # CCZ with X on a control: CCZ X_0 CCZ^dag = X_0 CZ_{12}.
    assert close(ccz @ pauli("XII") @ ccz.conj().T, pauli("XII") @ kron(I2, CZ))
    print("[ok] Toffoli/CCZ: X on control -> CNOT/CZ residue (frame form) or "
          "negative-control variant (substitution form); Z on Tof target -> extra CZ")

    # Measurement lemma: Pr_{X^a Z^b psi}(x) = Pr_psi(x xor a).
    n = 3
    for _ in range(10):
        psi = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
        psi /= np.linalg.norm(psi)
        for a in itertools.product((0, 1), repeat=n):
            for b in itertools.product((0, 1), repeat=n):
                p = kron(*(np.linalg.matrix_power(X, ai) @ np.linalg.matrix_power(Z, bi)
                           for ai, bi in zip(a, b)))
                pr = np.abs(p @ psi) ** 2
                aint = int("".join(map(str, a)), 2)
                assert close(pr, np.abs(psi[[x ^ aint for x in range(2 ** n)]]) ** 2)
    print("[ok] Lemma 3: measuring X^a Z^b |psi> gives outcome x with probability Pr_psi(x xor a)")

    # Theorem 3: table-only operations are product-preserving; CNOT is not.
    plus0 = kron(H, I2) @ np.array([1, 0, 0, 0], dtype=complex)
    bell = CNOT @ plus0
    schmidt = np.linalg.svd(bell.reshape(2, 2), compute_uv=False)
    assert np.sum(schmidt > 1e-10) == 2
    for p in itertools.product("IXYZ", repeat=2):
        for perm in (np.eye(4), SWAP):
            out = perm @ pauli("".join(p)) @ plus0
            assert np.sum(np.linalg.svd(out.reshape(2, 2), compute_uv=False) > 1e-10) == 1
    print("[ok] Theorem 3: every Pauli*permutation keeps |+>|0> a product state; CNOT makes a Bell state")


# ---------------------------------------------------------------- main theorem, end to end

GATE_1Q = {"H": H, "S": S, "Sdg": SDG, "T": T, "Tdg": TDG, "X": X, "Y": Y, "Z": Z}
GATE_2Q = {"CNOT": CNOT, "CZ": CZ, "SWAP": SWAP}
T_SWAP = {"T": "Tdg", "Tdg": "T"}


def apply(state, n, u, qubits):
    """Apply matrix u to the given qubits (qubit 0 = most significant)."""
    k = len(qubits)
    psi = state.reshape([2] * n)
    psi = np.moveaxis(psi, qubits, range(k))
    psi = (u @ psi.reshape(2 ** k, -1)).reshape([2] * n)
    return np.moveaxis(psi, range(k), qubits).reshape(-1)


def measure_branches(state, n, q):
    psi = state.reshape([2] * n)
    for b in (0, 1):
        branch = np.zeros_like(psi)
        idx = [slice(None)] * n
        idx[q] = b
        branch[tuple(idx)] = psi[tuple(idx)]
        pr = np.linalg.norm(branch) ** 2
        if pr > 1e-12:
            yield b, pr, branch.reshape(-1) / np.sqrt(pr)


def run_reference(ops, n):
    """Exact outcome distribution of the circuit as written (every gate applied)."""
    dist = {}

    def go(i, state, bits, pr):
        if i == len(ops):
            key = tuple(bits)
            dist[key] = dist.get(key, 0) + pr
            return
        op = ops[i]
        if op[0] == "measure":
            for b, p, s in measure_branches(state, n, op[1]):
                go(i + 1, s, bits + [b], pr * p)
            return
        if op[0] == "if":                      # ("if", bit_index, gate_op)
            if bits[op[1]] == 0:
                return go(i + 1, state, bits, pr)
            op = op[2]
        name, qs = op[0], list(op[1:])
        u = GATE_1Q.get(name) if len(qs) == 1 else GATE_2Q[name]
        go(i + 1, apply(state, n, u, qs), bits, pr)

    go(0, np.eye(2 ** n, dtype=complex)[0], [], 1.0)
    return dist


def run_compiled(ops, n):
    """Run with the address table + Pauli flags: no physical X, Y, Z or SWAP is ever applied.

    Returns (distribution of corrected outcomes, number of physical gates applied per branch).
    """
    dist = {}
    phys_counts = set()
    tables = {g: {r["in"]: r["out"] for r in e["rows"]}
              for g, e in json.loads(TABLES.read_text())["gates"].items()}

    def letter(x, z):
        return "IXZY"[x + 2 * z]

    def bits_of(c):
        return {"I": (0, 0), "X": (1, 0), "Z": (0, 1), "Y": (1, 1)}[c]

    def go(i, state, pi, fx, fz, bits, pr, nphys):
        if i == len(ops):
            key = tuple(bits)
            dist[key] = dist.get(key, 0) + pr
            phys_counts.add(nphys)
            return
        op = ops[i]
        if op[0] == "measure":
            q = op[1]
            for b, p, s in measure_branches(state, n, pi[q]):
                fz2 = list(fz)
                fz2[q] = 0                     # a Z flag on a measured qubit is a global phase
                go(i + 1, s, pi, fx, fz2, bits + [b ^ fx[q]], pr * p, nphys)
            return
        if op[0] == "if":
            if bits[op[1]] == 0:
                return go(i + 1, state, pi, fx, fz, bits, pr, nphys)
            op = op[2]
        name, qs = op[0], list(op[1:])
        fx, fz, pi = list(fx), list(fz), list(pi)
        if name in ("X", "Y", "Z"):            # rule P: flip flags only
            dx, dz = bits_of(name)
            fx[qs[0]] ^= dx
            fz[qs[0]] ^= dz
            return go(i + 1, state, pi, fx, fz, bits, pr, nphys)
        if name == "SWAP":                     # rule SWAP: edit the address table only
            a, b = qs
            pi[a], pi[b] = pi[b], pi[a]
            fx[a], fx[b] = fx[b], fx[a]
            fz[a], fz[b] = fz[b], fz[a]
            return go(i + 1, state, pi, fx, fz, bits, pr, nphys)
        phys = [pi[q] for q in qs]
        if name in T_SWAP:                     # rule T: substitute T <-> T^dag under an X or Y flag
            g = T_SWAP[name] if fx[qs[0]] else name
            return go(i + 1, apply(state, n, GATE_1Q[g], phys), pi, fx, fz, bits, pr, nphys + 1)
        # rule C: apply the Clifford on the mapped lines, push the flags through it
        u = GATE_1Q.get(name) if len(qs) == 1 else GATE_2Q[name]
        out = tables[name]["".join(letter(fx[q], fz[q]) for q in qs)] \
            if any(fx[q] or fz[q] for q in qs) else "I" * len(qs)
        for q, c in zip(qs, out):
            fx[q], fz[q] = bits_of(c)
        go(i + 1, apply(state, n, u, phys), pi, fx, fz, bits, pr, nphys + 1)

    go(0, np.eye(2 ** n, dtype=complex)[0], list(range(n)), [0] * n, [0] * n, [], 1.0, 0)
    return dist, phys_counts


def random_circuit(rng, n, depth, use_t, mid_measure):
    ones = ["H", "S", "Sdg", "X", "Y", "Z"] + (["T", "Tdg"] if use_t else [])
    ops, n_meas = [], 0
    for _ in range(depth):
        r = rng.random()
        if r < 0.5 or (n < 2 and r < 0.85):
            ops.append((rng.choice(ones), int(rng.integers(n))))
        elif r < 0.85:
            a, b = rng.choice(n, 2, replace=False)
            ops.append((rng.choice(["CNOT", "CZ", "SWAP"]), int(a), int(b)))
        elif mid_measure:
            ops.append(("measure", int(rng.integers(n))))
            n_meas += 1
            if rng.random() < 0.7:             # classically controlled Pauli or Clifford
                g = rng.choice(["X", "Z", "H", "S"])
                ops.append(("if", n_meas - 1, (g, int(rng.integers(n)))))
    for q in range(n):
        ops.append(("measure", q))
    return ops


def n_nonfree(ops):
    """Gates the compiled circuit must still apply physically, in the branch where every `if` fires."""
    count = 0
    for op in ops:
        if op[0] == "if":
            op = op[2]
        if op[0] not in ("measure", "X", "Y", "Z", "SWAP"):
            count += 1
    return count


def check_end_to_end():
    rng = np.random.default_rng(7)
    cases = [("Clifford", False, False), ("Clifford + mid-circuit measurement/feed-forward", False, True),
             ("Clifford+T", True, False), ("Clifford+T + mid-circuit measurement/feed-forward", True, True)]
    for label, use_t, mid in cases:
        for trial in range(150):
            n = int(rng.integers(1, 5))
            ops = random_circuit(rng, n, int(rng.integers(5, 30)), use_t, mid)
            ref = run_reference(ops, n)
            got, counts = run_compiled(ops, n)
            keys = set(ref) | set(got)
            assert all(abs(ref.get(k, 0) - got.get(k, 0)) < 1e-9 for k in keys), (label, ops)
            if not mid:
                assert counts == {n_nonfree(ops)}, (counts, n_nonfree(ops))
        print(f"[ok] Theorem 1/2 end to end, {label}: 150 random circuits, identical outcome "
              f"distributions with zero physical X/Y/Z/SWAP")


if __name__ == "__main__":
    if "--write" in sys.argv:
        TABLES.write_text(json.dumps(build_tables(), indent=1) + "\n")
        print(f"wrote {TABLES}")
    data = json.loads(TABLES.read_text())
    assert data == json.loads(json.dumps(build_tables())), "tables on disk differ from matrices"
    check_tables(data)
    check_lemmas()
    check_end_to_end()
    print("all checks passed")
