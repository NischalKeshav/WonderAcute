"""Compiler pass that absorbs X, Y, Z and SWAP gates into the address table.

``compress(circuit)`` walks the circuit once, in order, and keeps an ``AddressTable``:

* X, Y, Z, id     -> toggle frame bits, emit nothing
* SWAP            -> exchange two entries of the permutation, emit nothing
* Clifford G      -> emit G on the mapped lines; frame F becomes G F G^dagger
* diagonal non-Clifford (T, Tdg, RZ, P) on a line with pending X/Y:
    policy "rewrite": emit the angle-negated gate (X T X ~ Tdg), frame unchanged
    policy "flush":   emit the pending Pauli as a gate, clear it, then emit the gate
  pending Z commutes with diagonal gates, so it is kept either way
* any other gate  -> flush the pending Paulis on its lines, then emit it
* measure         -> emit on the mapped line; record a classical flip if X/Y is pending
* reset           -> emit on the mapped line; the frame on that line is cleared

The result carries the final table so that the frame can be applied to unmeasured outputs or
handed on to whatever runs next (``materialize`` appends it as gates for checking).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qiskit.circuit import Instruction, QuantumCircuit
from qiskit.quantum_info import Operator
from qiskit.circuit.library import (CPhaseGate, PhaseGate, RZGate, RZZGate, TdgGate, TGate,
                                    XGate, YGate, ZGate)
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler.basepasses import TransformationPass

from .frame import (ANGLE_FLIP, ANGLE_NEGATE, PAULI_GATES, PERMUTATION_GATES, AddressTable,
                    is_clifford)

PAULI_GATE = {"X": XGate, "Y": YGate, "Z": ZGate}
POLICIES = ("rewrite", "flush")


@dataclass
class CompressionResult:
    circuit: QuantumCircuit
    table: AddressTable
    # One entry per measurement, in emission order: (physical qubit, clbit index, flipped?).
    measurements: list[tuple[int, int, bool]] = field(default_factory=list)
    removed: dict[str, int] = field(default_factory=dict)
    flushed: int = 0  # Pauli gates that had to be emitted before non-Clifford gates
    rewritten_2q: int = 0  # CP/RZZ gates rewritten under a pending X (adds a P gate for CP)

    @property
    def clbit_flips(self) -> dict[int, bool]:
        """Final classical correction per clbit (the last measurement into a clbit wins)."""
        flips: dict[int, bool] = {}
        for _, c, f in self.measurements:
            flips[c] = f
        return flips


def _flush(out: QuantumCircuit, table: AddressTable, physical: list[int]) -> int:
    n = 0
    for p in physical:
        label = table.pending(p)
        if label != "I":
            out.append(PAULI_GATE[label](), [p])
            table.clear(p)
            n += 1
    return n


def compress(circuit: QuantumCircuit, policy: str = "rewrite",
             rewrite_cp: bool = False) -> CompressionResult:
    """``rewrite_cp`` applies the CP/RZZ rewrite under policy "rewrite". It trades one flushed X
    for one extra P gate, which is free on virtual-Z hardware and a cost otherwise, so it is off
    by default (the plan's metric is total gate count); compare ``gates`` and ``gates_vz``."""
    if policy not in POLICIES:
        raise ValueError(f"policy must be one of {POLICIES}")
    qindex = {q: i for i, q in enumerate(circuit.qubits)}
    cindex = {c: i for i, c in enumerate(circuit.clbits)}
    table = AddressTable(circuit.num_qubits)
    out = circuit.copy_empty_like()
    result = CompressionResult(out, table)
    removed = result.removed

    for inst in circuit.data:
        op = inst.operation
        name = op.name
        logical = [qindex[q] for q in inst.qubits]
        physical = [table.phys[q] for q in logical]

        if getattr(op, "blocks", None):
            raise NotImplementedError(f"control flow ({name}) is not supported")

        if name in PAULI_GATES:
            table.apply_pauli(name, logical[0])
            removed[name] = removed.get(name, 0) + 1
            continue
        if name in PERMUTATION_GATES:
            table.swap(*logical)
            removed[name] = removed.get(name, 0) + 1
            continue

        qargs = [out.qubits[p] for p in physical]
        if name == "measure":
            (p,) = physical
            result.measurements.append((p, cindex[inst.clbits[0]], table.measurement_flip(p)))
            table.z[p] = False  # Z on a computational-basis state is a global phase
            out.append(op, qargs, inst.clbits)
            continue
        if name == "reset":
            table.clear(physical[0])
            out.append(op, qargs, inst.clbits)
            continue
        if name in ("barrier", "delay"):
            out.append(op, qargs, inst.clbits)
            continue

        if is_clifford(op):
            out.append(op, qargs, inst.clbits)
            table.conjugate_through(op, physical)
            continue

        if len(physical) == 1 and (name in ANGLE_FLIP or name in ANGLE_NEGATE):
            (p,) = physical
            if table.x[p]:
                if policy == "rewrite":
                    if name in ANGLE_FLIP:
                        new = TGate() if ANGLE_FLIP[name] == "t" else TdgGate()
                    else:
                        new = (RZGate if name == "rz" else PhaseGate)(-op.params[0])
                    out.append(new, qargs)
                    continue
                # flush only the X part; a pending Z commutes with the diagonal gate
                out.append(XGate(), [p])
                table.x[p] = False
                result.flushed += 1
            out.append(op, qargs, inst.clbits)
            continue

        if policy == "rewrite" and rewrite_cp and name in ("cp", "cu1", "rzz") and len(physical) == 2:
            a, b = physical
            xa, xb = bool(table.x[a]), bool(table.x[b])
            theta = op.params[0]
            if xa or xb:
                if name == "rzz":
                    # X_a RZZ(t) X_a = RZZ(-t); X on both leaves it unchanged
                    out.append(RZZGate(theta if (xa and xb) else -theta), qargs)
                elif xa and xb:
                    # X_a X_b CP(t) X_a X_b = e^{it} P_a(-t) P_b(-t) CP(t)
                    out.append(PhaseGate(-theta), [out.qubits[a]])
                    out.append(PhaseGate(-theta), [out.qubits[b]])
                    out.append(CPhaseGate(theta), qargs)
                else:
                    # X_a CP(t) X_a = CP(-t) P_b(t)   (proofs/formal-model.md section 7)
                    other = b if xa else a
                    out.append(CPhaseGate(-theta), qargs)
                    out.append(PhaseGate(theta), [out.qubits[other]])
                result.rewritten_2q += 1
                continue

        if is_diagonal(op):
            # pending Z commutes with any diagonal gate; only the X parts must be flushed
            for p in physical:
                if table.x[p]:
                    out.append(XGate(), [p])
                    table.x[p] = False
                    result.flushed += 1
        else:
            result.flushed += _flush(out, table, physical)
        out.append(op, qargs, inst.clbits)

    return result


DIAGONAL_NAMES = {"p", "rz", "t", "tdg", "z", "s", "sdg", "u1", "cp", "cz", "ccz", "cu1",
                  "crz", "rzz", "mcphase", "mcp", "mcu1"}


_diag_cache: dict[tuple, bool] = {}


def is_diagonal(op: Instruction) -> bool:
    """Diagonal in the computational basis: checked by name, else by matrix for <= 8 qubits."""
    if op.name in DIAGONAL_NAMES:
        return True
    if op.num_qubits > 8:
        return False
    key = (op.name, repr(op.params), op.num_qubits)
    if key not in _diag_cache:
        try:
            m = Operator(op).data
            _diag_cache[key] = bool(np.allclose(m - np.diag(np.diag(m)), 0))
        except Exception:
            _diag_cache[key] = False
    return _diag_cache[key]


def materialize(result: CompressionResult, include_permutation: bool = True) -> QuantumCircuit:
    """Compressed circuit followed by the final frame as gates (and swaps undoing the
    permutation), so it equals the source circuit up to global phase on unmeasured lines."""
    qc = result.circuit.copy()
    table = result.table
    for p in range(table.num_qubits):
        label = table.pending(p)
        if label != "I":
            qc.append(PAULI_GATE[label](), [p])
    if include_permutation:
        # Move the content of physical phys[q] to wire q with swaps.
        where = list(table.phys)              # where[q] = wire currently holding logical q
        holder = {p: q for q, p in enumerate(where)}
        for q in range(table.num_qubits):
            p = where[q]
            if p != q:
                other = holder[q]             # logical line currently on wire q
                qc.swap(q, p)
                where[q], where[other] = q, p
                holder[q], holder[p] = q, other
    return qc


class AddressTablePass(TransformationPass):
    """Qiskit pass wrapper. Stores the final table in ``property_set['address_table']`` and the
    per-measurement flips in ``property_set['measurement_flips']``."""

    def __init__(self, policy: str = "rewrite", rewrite_cp: bool = False):
        super().__init__()
        self.policy = policy
        self.rewrite_cp = rewrite_cp

    def run(self, dag):
        res = compress(dag_to_circuit(dag), self.policy, self.rewrite_cp)
        self.property_set["address_table"] = res.table
        self.property_set["measurement_flips"] = res.measurements
        self.property_set["absorbed_gates"] = res.removed
        return circuit_to_dag(res.circuit)
