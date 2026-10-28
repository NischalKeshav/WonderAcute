"""Clock-tick scheduling and a noise model with idle error (step 6 of the research plan).

``schedule_ticks(circuit, ops_per_tick)`` turns a circuit into a list of ticks. Each tick holds
at most ``ops_per_tick`` gates from one dependency layer (a layer that does not fit is split
across consecutive ticks). ``with_idle(circuit, ops_per_tick)`` rebuilds the circuit with an
``id`` gate on every qubit that is idle during a tick, so an Aer noise model can attach an
idle channel to it. Measurements are left at the end and are not scheduled.

``noise_model(p1, p2, p_idle)``: amplitude damping and dephasing with probability ``p1`` after
every one-qubit gate, the same channel on both qubits with probability ``p2`` after a
two-qubit gate, and with probability ``p_idle`` on every idle qubit per tick.

There is no crosstalk term, so packing more operations into a tick is always better here;
an optimum number of operations per tick only appears once simultaneous gates are penalised.
"""

from __future__ import annotations

from qiskit.circuit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from qiskit_aer.noise import NoiseModel, amplitude_damping_error, phase_damping_error

ONE_Q = ["h", "s", "sdg", "sx", "sxdg", "x", "y", "z", "t", "tdg", "rz", "p", "u", "u3", "rx", "ry"]
TWO_Q = ["cx", "cz", "cy", "swap", "cp", "rzz", "ecr"]


def schedule_ticks(qc: QuantumCircuit, ops_per_tick: int | None):
    body = qc.remove_final_measurements(inplace=False)
    dag = circuit_to_dag(body)
    ticks = []
    for layer in dag.layers():
        ops = [n for n in layer["graph"].op_nodes() if n.op.name not in ("barrier", "measure")]
        if not ops:
            continue
        if ops_per_tick is None:
            ticks.append(ops)
        else:
            for i in range(0, len(ops), ops_per_tick):
                ticks.append(ops[i:i + ops_per_tick])
    return ticks


def with_idle(qc: QuantumCircuit, ops_per_tick: int | None) -> QuantumCircuit:
    out = qc.copy_empty_like()
    qi = {q: i for i, q in enumerate(qc.qubits)}
    for ops in schedule_ticks(qc, ops_per_tick):
        busy = set()
        for n in ops:
            out.append(n.op, [out.qubits[qi[q]] for q in n.qargs], n.cargs)
            busy.update(qi[q] for q in n.qargs)
        for i in range(qc.num_qubits):
            if i not in busy:
                out.id(out.qubits[i])
    # measurements, in their original order
    for inst in qc.data:
        if inst.operation.name == "measure":
            out.append(inst.operation, [out.qubits[qi[q]] for q in inst.qubits], inst.clbits)
    return out


def noise_model(p1: float, p2: float, p_idle: float, channel: str = "both") -> NoiseModel:
    """``channel``: "both" (amplitude damping then dephasing), "dephasing" or "damping" only.
    Separating them matters: a Pauli frame keeps the physical register in |0> while the
    answer lives in classical bits, and amplitude damping cannot touch |0>."""
    def chan(p):
        if channel == "dephasing":
            return phase_damping_error(p)
        if channel == "damping":
            return amplitude_damping_error(p)
        return amplitude_damping_error(p).compose(phase_damping_error(p))

    nm = NoiseModel(basis_gates=ONE_Q + TWO_Q + ["id"])
    nm.add_all_qubit_quantum_error(chan(p1), ONE_Q)
    nm.add_all_qubit_quantum_error(chan(p2).tensor(chan(p2)), TWO_Q)
    nm.add_all_qubit_quantum_error(chan(p_idle), ["id"])
    return nm
