"""Memory-lookup circuits for step 4 of the research plan (the oracle is not free).

Two constructions, both reading a classical table ``data`` of N = 2^k entries of ``d`` bits
into an output register, controlled on a k-qubit address in superposition:

``bucket_brigade``  Giovannetti–Lloyd–Maccone tree, in the qubit-router circuit model used by
    Arunachalam et al. (2015) and Di Matteo et al. (2019). Each internal node has a router
    qubit and an input qubit. Address bits are sent down the tree one at a time with
    controlled-SWAPs (2 per node on the way, with an X pair to select the other child),
    each bit becoming the router of its level. A bus then descends, copies the addressed
    cell, and ascends; finally the address routing is uncomputed.
        CSWAP = 4(N-1-k) + 8d(N-1)  (= 12N - 4k - 12 for d = 1),  CNOT(data copy) = 2N·d,
    X = CSWAP,  qubits = k + d + 2(N-1) + N + N·d + d.

``table_lookup``  Babbush et al. (2018) unary-iteration QROM (Fig. 7), written out with the
    simple compute/uncompute Toffoli pair per node rather than the measurement-based
    uncomputation, so it costs 2(N-2) Toffolis instead of N-1; the N-1 figure is given as the
    closed form.  Ancillas: k-1.  Data copy: one CNOT per set bit.

``select_swap_counts`` closed form only, written from memory of Low, Kliuchnikov,
    Schaeffer (2018) and not yet checked against the paper: Toffoli
    ≈ ceil(N/λ) + d(λ-1) with d·λ extra qubits, λ a power of two.

All three are permutation circuits (X, CX, SWAP, CSWAP, CCX), so ``classical_run`` checks
them exactly on every address without a quantum simulator. The X gates that select a child
are exactly the gates the address-table pass absorbs, which links step 4 back to step 2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from qiskit.circuit import AncillaRegister, QuantumCircuit, QuantumRegister

from .circuits import toffoli


@dataclass
class Lookup:
    circuit: QuantumCircuit
    address: list[int]
    output: list[int]
    N: int
    d: int


def _bits(k: int, v: int) -> list[int]:
    return [(v >> i) & 1 for i in range(k)]


# --------------------------------------------------------------------------- bucket brigade

def bucket_brigade(data: list[int], d: int = 1, exact: bool = False) -> Lookup:
    """``exact`` lowers every CSWAP to CNOT + exact Clifford+T Toffoli + CNOT."""
    N = len(data)
    k = int(math.log2(N))
    assert 2 ** k == N and k >= 1
    addr = QuantumRegister(k, "a")
    out = QuantumRegister(d, "out")
    # node i (heap order, root 0) has router r[i] and input in_[i]; leaves hold the cells
    n_internal = N - 1
    r = QuantumRegister(n_internal, "r")
    inp = QuantumRegister(n_internal, "in")
    leaf = QuantumRegister(N, "leaf")          # bus arrives here at the leaves
    cell = QuantumRegister(N * d, "cell")      # memory contents, prepared classically
    bus = QuantumRegister(d, "bus")
    qc = QuantumCircuit(addr, out, r, inp, leaf, cell, bus)

    for i, v in enumerate(data):
        for b in range(d):
            if (v >> b) & 1:
                qc.x(cell[i * d + b])

    def children(i):
        return 2 * i + 1, 2 * i + 2

    def child_slot(c):  # where the bus lands when sent to child c
        return inp[c] if c < n_internal else leaf[c - n_internal]

    def cswap(ctrl, a, b):
        if exact:
            qc.cx(b, a); toffoli(qc, ctrl, a, b); qc.cx(b, a)
        else:
            qc.cswap(ctrl, a, b)

    def route_down(level_from: int, level_to: int, payload_at_root):
        """Send what sits in the root's input down to depth ``level_to`` (all nodes in
        superposition), through routers at depths < level_to."""
        for depth in range(level_from, level_to):
            for i in range(2 ** depth - 1, 2 ** (depth + 1) - 1):
                left, right = children(i)
                cswap(r[i], inp[i], child_slot(right))
                qc.x(r[i])
                cswap(r[i], inp[i], child_slot(left))
                qc.x(r[i])

    def route_up(level_to: int):
        for depth in reversed(range(level_to)):
            for i in reversed(range(2 ** depth - 1, 2 ** (depth + 1) - 1)):
                left, right = children(i)
                qc.x(r[i])
                cswap(r[i], inp[i], child_slot(left))
                qc.x(r[i])
                cswap(r[i], inp[i], child_slot(right))

    def set_routers(depth: int, forward: bool):
        # the address bit that arrived in in_[i] becomes router r[i]
        for i in range(2 ** depth - 1, 2 ** (depth + 1) - 1):
            qc.swap(inp[i], r[i])

    # phase 1: load address bits top-down (bit j = most significant selects the root)
    for j in range(k):
        bit = addr[k - 1 - j]
        qc.swap(bit, inp[0])
        route_down(0, j, inp[0])
        set_routers(j, True)

    # phase 2: data bus, one bit at a time (d passes), descend, copy, ascend
    for b in range(d):
        qc.swap(bus[b], inp[0])
        route_down(0, k, inp[0])
        for i in range(N):
            qc.cx(cell[i * d + b], leaf[i])
        route_up(k)
        qc.swap(bus[b], inp[0])
        qc.cx(bus[b], out[b])  # bus now holds data[address]; copy it out, then uncompute
        # uncompute the bus value by reading again (the lookup is its own inverse on the bus)
        qc.swap(bus[b], inp[0])
        route_down(0, k, inp[0])
        for i in range(N):
            qc.cx(cell[i * d + b], leaf[i])
        route_up(k)
        qc.swap(bus[b], inp[0])

    # phase 3: return the address bits
    for j in reversed(range(k)):
        bit = addr[k - 1 - j]
        set_routers(j, False)
        route_up(j)
        qc.swap(bit, inp[0])

    idx = {q: i for i, q in enumerate(qc.qubits)}
    return Lookup(qc, [idx[q] for q in addr], [idx[q] for q in out], N, d)


def bucket_brigade_counts(N: int, d: int = 1) -> dict:
    k = int(math.log2(N))
    cswap_addr = 2 * sum(2 ** j - 1 for j in range(k))      # routing address bits down
    cswap_bus = 2 * (N - 1)                                 # one descent (or ascent)
    cswap = 2 * cswap_addr + d * 4 * cswap_bus              # addr down+up, bus 2 descents + 2 ascents per bit
    return {"N": N, "d": d, "qubits": k + d + 2 * (N - 1) + N + N * d + d,
            "cswap": cswap, "toffoli": cswap, "t_count": 7 * cswap,
            "cnot_data": 2 * N * d, "x_absorbable": cswap}


# ---------------------------------------------------------------------------- table lookup

def table_lookup(data: list[int], d: int = 1, exact: bool = False) -> Lookup:
    """``exact`` lowers every Toffoli to the 7-T Clifford+T circuit."""
    N = len(data)
    k = int(math.log2(N))
    assert 2 ** k == N and k >= 1
    addr = QuantumRegister(k, "a")
    out = QuantumRegister(d, "out")
    anc = AncillaRegister(max(k - 1, 0), "anc")
    qc = QuantumCircuit(addr, out, anc)

    def tof(a, b, c):
        if exact:
            toffoli(qc, a, b, c)
        else:
            qc.ccx(a, b, c)

    def write(idx: int, ctrl):
        for b in range(d):
            if (data[idx] >> b) & 1:
                qc.cx(ctrl, out[b])

    def iterate(ctrl, level: int, idx: int):
        # ctrl is set iff the address bits above ``level`` match idx's prefix
        if level < 0:
            write(idx, ctrl)
            return
        bit = addr[level]
        a = anc[level]  # ancilla for this level (levels k-2 .. 0)
        qc.x(bit)
        tof(ctrl, bit, a)                  # a = ctrl & ~bit
        qc.x(bit)
        iterate(a, level - 1, 2 * idx)
        qc.cx(ctrl, a)                     # a = ctrl & bit
        iterate(a, level - 1, 2 * idx + 1)
        tof(ctrl, bit, a)                  # uncompute

    top = addr[k - 1]
    if k == 1:
        qc.x(top); write(0, top); qc.x(top); write(1, top)
    else:
        qc.x(top)
        iterate(top, k - 2, 0)
        qc.x(top)
        iterate(top, k - 2, 1)
    idx = {q: i for i, q in enumerate(qc.qubits)}
    return Lookup(qc, [idx[q] for q in addr], [idx[q] for q in out], N, d)


def table_lookup_counts(N: int, d: int = 1) -> dict:
    k = int(math.log2(N))
    toff = 2 * (N - 2) if k >= 2 else 0
    return {"N": N, "d": d, "qubits": k + d + max(k - 1, 0),
            "toffoli": toff, "toffoli_babbush": N - 1, "t_count": 7 * toff,
            "t_count_babbush": 4 * (N - 1), "cnot_data_max": N * d,
            "x_absorbable": 2 * (N - 1)}


def select_swap_counts(N: int, d: int = 1, lam: int | None = None) -> dict:
    """Closed form for the Low–Kliuchnikov–Schaeffer select-swap QROM. λ defaults to the
    power of two nearest sqrt(N/d)."""
    if lam is None:
        lam = 2 ** max(0, round(math.log2(max(1, math.sqrt(N / d)))))
    k = int(math.log2(N))
    toff = math.ceil(N / lam) + d * (lam - 1)
    return {"N": N, "d": d, "lambda": lam, "qubits": k + d * lam + max(k - 1, 0),
            "toffoli": toff, "t_count": 4 * toff}


# ------------------------------------------------------------------- exact classical check

def classical_run(qc: QuantumCircuit, init: dict[int, int]) -> np.ndarray:
    """Run a permutation circuit (x, cx, ccx, swap, cswap) on a basis state."""
    idx = {q: i for i, q in enumerate(qc.qubits)}
    s = np.zeros(qc.num_qubits, dtype=np.uint8)
    for q, v in init.items():
        s[q] = v
    for inst in qc.data:
        n = inst.operation.name
        qs = [idx[q] for q in inst.qubits]
        if n == "x":
            s[qs[0]] ^= 1
        elif n == "cx":
            s[qs[1]] ^= s[qs[0]]
        elif n == "ccx":
            s[qs[2]] ^= s[qs[0]] & s[qs[1]]
        elif n == "swap":
            s[qs[0]], s[qs[1]] = s[qs[1]], s[qs[0]]
        elif n == "cswap":
            if s[qs[0]]:
                s[qs[1]], s[qs[2]] = s[qs[2]], s[qs[1]]
        elif n in ("h", "t", "tdg"):
            raise ValueError("not a permutation circuit; build with exact=False")
        elif n in ("barrier",):
            continue
        else:
            raise ValueError(f"unsupported gate {n}")
    return s


def check_lookup(lk: Lookup, data: list[int]) -> bool:
    """Every address reads its entry, the address comes back unchanged, and every other
    qubit is clean except the memory cells, which hold the data they were prepared with."""
    k = len(lk.address)
    qc = lk.circuit
    cell = next((r for r in qc.qregs if r.name == "cell"), None)
    cell_idx = {qc.find_bit(q).index for q in cell} if cell else set()
    for a in range(lk.N):
        init = {q: b for q, b in zip(lk.address, _bits(k, a))}
        s = classical_run(qc, init)
        got = sum(int(s[q]) << b for b, q in enumerate(lk.output))
        if got != data[a]:
            return False
        for i in range(qc.num_qubits):
            if i in lk.output or i in cell_idx:
                continue
            if s[i] != init.get(i, 0):
                return False
    return True
