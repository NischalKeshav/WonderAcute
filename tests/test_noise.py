from qiskit import QuantumCircuit

from qarch import circuits
from qarch.noise import noise_model, schedule_ticks, with_idle


def test_ticks_respect_ops_per_tick_and_idle_ids():
    qc = circuits.ghz(4)                      # h, cx, cx, cx (a chain: 4 layers)
    assert len(schedule_ticks(qc, None)) == 4
    qc2 = QuantumCircuit(4, 4)
    qc2.h(range(4))                          # one layer of 4 ops
    qc2.measure(range(4), range(4))
    assert len(schedule_ticks(qc2, 2)) == 2
    sched = with_idle(qc2, 2)
    names = [i.operation.name for i in sched.data]
    assert names.count("id") == 4            # two ticks, two idle qubits each
    assert names.count("measure") == 4 and names[-1] == "measure"


def test_noise_model_covers_gates():
    nm = noise_model(1e-3, 1e-2, 1e-3)
    assert {"id", "cx", "h"} <= set(nm.noise_instructions)
