"""Benchmark circuit generators. All randomness goes through an explicit seed.

Gate-set choice matters for the savings numbers: the pass can only remove X/Y/Z/SWAP gates that
are present. ``random_layered`` draws from a stated distribution (``DEFAULT_LAYER_PROBS``) so the
Pauli and SWAP density is a documented input; ``random_clifford_synth`` uses Qiskit's synthesis of
a uniformly random Clifford, which emits whatever gates the synthesis chooses.
"""

from __future__ import annotations

import math

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import MCMTGate, ZGate
from qiskit.quantum_info import random_clifford

# Per-slot probabilities for random_layered. One-qubit slots pick from the 1q gates; with
# probability p2 a pair of free qubits gets a two-qubit gate instead.
DEFAULT_LAYER_PROBS = {
    "1q": {"h": 0.25, "s": 0.15, "sdg": 0.1, "x": 0.15, "y": 0.1, "z": 0.15, "sx": 0.1},
    "2q": {"cx": 0.6, "cz": 0.2, "swap": 0.2},
    "p2": 0.3,
}


def random_layered(n: int, depth: int, seed: int, probs: dict | None = None,
                   t_prob: float = 0.0, measure: bool = False) -> QuantumCircuit:
    """Layered random circuit. ``t_prob`` mixes in T/Tdg on one-qubit slots (non-Clifford)."""
    probs = probs or DEFAULT_LAYER_PROBS
    rng = np.random.default_rng(seed)
    g1, w1 = zip(*probs["1q"].items())
    g2, w2 = zip(*probs["2q"].items())
    w1 = np.array(w1) / sum(w1)
    w2 = np.array(w2) / sum(w2)
    qc = QuantumCircuit(n, n if measure else 0)
    for _ in range(depth):
        order = rng.permutation(n)
        i = 0
        while i < n:
            if i + 1 < n and rng.random() < probs["p2"]:
                a, b = int(order[i]), int(order[i + 1])
                getattr(qc, rng.choice(g2, p=w2))(a, b)
                i += 2
            else:
                q = int(order[i])
                if t_prob and rng.random() < t_prob:
                    getattr(qc, rng.choice(["t", "tdg"]))(q)
                else:
                    getattr(qc, rng.choice(g1, p=w1))(q)
                i += 1
    if measure:
        qc.measure(range(n), range(n))
    return qc


def random_clifford_synth(n: int, seed: int) -> QuantumCircuit:
    return random_clifford(n, seed=seed).to_circuit()


def ghz(n: int, measure: bool = True) -> QuantumCircuit:
    qc = QuantumCircuit(n, n if measure else 0)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    if measure:
        qc.measure(range(n), range(n))
    return qc


def bernstein_vazirani(n: int, seed: int) -> QuantumCircuit:
    """n input qubits plus one ancilla prepared in |->  with X then H (the X is absorbable)."""
    secret = np.random.default_rng(seed).integers(0, 2, n)
    qc = QuantumCircuit(n + 1, n)
    qc.x(n)
    qc.h(range(n + 1))
    for i, bit in enumerate(secret):
        if bit:
            qc.cx(i, n)
    qc.h(range(n))
    qc.measure(range(n), range(n))
    qc.metadata = {"secret": "".join(map(str, secret[::-1]))}
    return qc


def repetition_syndrome(d: int, rounds: int) -> QuantumCircuit:
    """Distance-d bit-flip repetition code: d data qubits, d-1 ancillas, measure+reset rounds."""
    n = 2 * d - 1
    qc = QuantumCircuit(n, rounds * (d - 1) + d)
    data = list(range(d))
    anc = list(range(d, n))
    k = 0
    for _ in range(rounds):
        for i, a in enumerate(anc):
            qc.cx(data[i], a)
            qc.cx(data[i + 1], a)
        for a in anc:
            qc.measure(a, k)
            qc.reset(a)
            k += 1
    qc.measure(data, range(k, k + d))
    return qc


def qft(n: int, swaps: bool = True, measure: bool = False) -> QuantumCircuit:
    """Textbook QFT: H and controlled phases, then the bit-reversal SWAP network (floor(n/2) SWAPs)."""
    qc = QuantumCircuit(n, n if measure else 0)
    for j in reversed(range(n)):
        qc.h(j)
        for k in reversed(range(j)):
            qc.cp(math.pi / 2 ** (j - k), k, j)
    if swaps:
        for i in range(n // 2):
            qc.swap(i, n - 1 - i)
    if measure:
        qc.measure(range(n), range(n))
    return qc


def _mcz(qc: QuantumCircuit, qubits: list[int], decompose: bool) -> None:
    if len(qubits) == 1:
        qc.z(qubits[0])
        return
    gate = MCMTGate(ZGate(), len(qubits) - 1, 1)
    if decompose:
        sub = QuantumCircuit(len(qubits))
        sub.append(gate, range(len(qubits)))
        from qiskit import transpile
        sub = transpile(sub, basis_gates=["h", "s", "sdg", "t", "tdg", "cx", "x"],
                        optimization_level=0)
        qc.compose(sub, qubits, inplace=True)
    else:
        qc.append(gate, qubits)


def grover(n: int, marked: int, iterations: int | None = None,
           decompose_mcz: bool = False) -> QuantumCircuit:
    """Grover search for one marked n-bit item. The oracle and diffusion are written with
    X conjugations around a multi-controlled Z, the usual textbook form. With
    ``decompose_mcz`` the MCZ is lowered to Clifford+T; otherwise it is kept as one gate."""
    if iterations is None:
        iterations = max(1, int(math.floor(math.pi / 4 * math.sqrt(2 ** n))))
    qs = list(range(n))
    qc = QuantumCircuit(n, n)
    qc.h(qs)
    for _ in range(iterations):
        zeros = [q for q in qs if not (marked >> q) & 1]
        for q in zeros:
            qc.x(q)
        _mcz(qc, qs, decompose_mcz)
        for q in zeros:
            qc.x(q)
        qc.h(qs)
        qc.x(qs)
        _mcz(qc, qs, decompose_mcz)
        qc.x(qs)
        qc.h(qs)
    qc.measure(qs, qs)
    return qc


def shor15(a: int = 7, n_count: int = 8) -> QuantumCircuit:
    """Textbook order-finding circuit for N = 15 (Qiskit textbook construction).
    Controlled modular multiplication by a^(2^j) mod 15 is a controlled permutation, built
    from controlled SWAPs and CNOTs; the inverse QFT includes its SWAP network."""
    if a not in (2, 4, 7, 8, 11, 13):
        raise ValueError("a must be coprime to 15 with a non-trivial order")

    def c_amod15(power: int) -> QuantumCircuit:
        u = QuantumCircuit(5)  # qubit 0 = control, 1..4 = work register
        w = [1, 2, 3, 4]
        for _ in range(power):
            if a in (2, 13):
                u.cswap(0, w[2], w[3]); u.cswap(0, w[1], w[2]); u.cswap(0, w[0], w[1])
            if a in (7, 8):
                u.cswap(0, w[0], w[1]); u.cswap(0, w[1], w[2]); u.cswap(0, w[2], w[3])
            if a in (4, 11):
                u.cswap(0, w[1], w[3]); u.cswap(0, w[0], w[2])
            if a in (7, 11, 13):
                for q in w:
                    u.cx(0, q)
        return u

    qc = QuantumCircuit(n_count + 4, n_count)
    qc.h(range(n_count))
    qc.x(n_count)  # work register starts in |1>
    for j in range(n_count):
        qc.compose(c_amod15(2 ** j), [j] + list(range(n_count, n_count + 4)), inplace=True)
    qc.compose(qft(n_count).inverse(), range(n_count), inplace=True)
    qc.measure(range(n_count), range(n_count))
    return qc
