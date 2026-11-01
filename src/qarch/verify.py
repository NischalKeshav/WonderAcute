"""Equivalence checks between a source circuit and its compressed form.

* ``unitary_equivalent``: exact operator comparison (up to global phase), small n.
* ``tableau_equivalent``: Stim tableau comparison for Clifford circuits, any n.
* ``measurements_equivalent``: Stim lock-step check for Clifford circuits with measurements.
  Both circuits are simulated side by side; at every measurement the two outcomes must have the
  same distribution (both random, or both deterministic with values differing by the recorded
  flip), and the compressed run is postselected on (source outcome XOR flip) so the two runs stay
  correlated. This is exact, not a sampling test.
"""

from __future__ import annotations

import numpy as np
import stim
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Operator

from .passes import CompressionResult, materialize

STIM_GATE = {
    "h": "H", "s": "S", "sdg": "S_DAG", "sx": "SQRT_X", "sxdg": "SQRT_X_DAG",
    "x": "X", "y": "Y", "z": "Z", "id": "I",
    "cx": "CX", "cz": "CZ", "cy": "CY", "swap": "SWAP", "iswap": "ISWAP", "ecr": None,
}


def unitary_equivalent(source: QuantumCircuit, result: CompressionResult) -> bool:
    return Operator(source.remove_final_measurements(inplace=False)).equiv(
        Operator(materialize(result).remove_final_measurements(inplace=False)))


def _stim_ops(qc: QuantumCircuit):
    """Yield ('gate', stim_name, targets) | ('measure', q, clbit) | ('reset', q, None)."""
    qi = {q: i for i, q in enumerate(qc.qubits)}
    ci = {c: i for i, c in enumerate(qc.clbits)}
    for inst in qc.data:
        name = inst.operation.name
        qs = [qi[q] for q in inst.qubits]
        if name == "measure":
            yield "measure", qs[0], ci[inst.clbits[0]]
        elif name == "reset":
            yield "reset", qs[0], None
        elif name in ("barrier", "delay"):
            continue
        else:
            sname = STIM_GATE.get(name)
            if sname is None:
                raise ValueError(f"gate {name} has no Stim equivalent")
            yield "gate", sname, qs


def to_stim(qc: QuantumCircuit) -> stim.Circuit:
    c = stim.Circuit()
    for kind, a, b in _stim_ops(qc):
        if kind == "gate":
            c.append(a, b)
        elif kind == "measure":
            c.append("M", [a])
        else:
            c.append("R", [a])
    return c


def tableau_equivalent(source: QuantumCircuit, result: CompressionResult) -> bool:
    src = source.remove_final_measurements(inplace=False)
    out = materialize(result).remove_final_measurements(inplace=False)
    n = source.num_qubits
    ta = stim.Tableau.from_circuit(to_stim(src), ignore_measurement=True)
    tb = stim.Tableau.from_circuit(to_stim(out), ignore_measurement=True)
    # from_circuit sizes the tableau by the largest qubit touched
    if len(ta) < n:
        ta = ta + stim.Tableau(n - len(ta))
    if len(tb) < n:
        tb = tb + stim.Tableau(n - len(tb))
    return ta == tb


def measurements_equivalent(source: QuantumCircuit, result: CompressionResult,
                            seed: int = 0) -> bool:
    n = source.num_qubits
    src_ops = list(_stim_ops(source))
    out_ops = list(_stim_ops(result.circuit))
    sa = stim.TableauSimulator(seed=seed)
    sb = stim.TableauSimulator(seed=seed)
    sa.set_num_qubits(n)
    sb.set_num_qubits(n)
    flips = iter(result.measurements)

    def run_until_measure(sim, ops, i):
        while i < len(ops):
            kind, a, b = ops[i]
            if kind == "measure":
                return i
            if kind == "reset":
                sim.reset(a)
            else:
                sim.do(stim.CircuitInstruction(a, b))
            i += 1
        return i

    i = j = 0
    while True:
        i = run_until_measure(sa, src_ops, i)
        j = run_until_measure(sb, out_ops, j)
        if i == len(src_ops) or j == len(out_ops):
            return i == len(src_ops) and j == len(out_ops) and next(flips, None) is None
        _, qa, ca = src_ops[i]
        _, qb, cb = out_ops[j]
        p, c, flip = next(flips)
        if (p, c) != (qb, cb) or ca != cb:
            return False
        pa, pb = sa.peek_z(qa), sb.peek_z(qb)
        if (pa == 0) != (pb == 0):
            return False
        if pa != 0 and (pa == pb) == flip:
            return False
        m = sa.measure(qa)
        try:
            sb.postselect_z(qb, desired_value=m ^ flip)
        except ValueError:
            return False
        i += 1
        j += 1


def counts_with_flips(counts: dict[str, int], result: CompressionResult,
                      num_clbits: int) -> dict[str, int]:
    """Apply the classical corrections to Qiskit counts (bitstrings are little-endian)."""
    mask = np.zeros(num_clbits, dtype=bool)
    for c, f in result.clbit_flips.items():
        mask[c] = f
    fixed: dict[str, int] = {}
    for key, v in counts.items():
        bits = [b == "1" for b in key.replace(" ", "")[::-1]]
        bits = [b ^ mask[i] for i, b in enumerate(bits)]
        k = "".join("1" if b else "0" for b in bits[::-1])
        fixed[k] = fixed.get(k, 0) + v
    return fixed
