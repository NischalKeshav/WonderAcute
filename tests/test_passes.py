"""The compiler pass preserves behaviour: exact operators (small n), Stim tableaux (large n),
and Stim lock-step measurement checks. All seeds are fixed."""

import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import PassManager

from qarch import circuits
from qarch.passes import AddressTablePass, compress
from qarch.verify import measurements_equivalent, tableau_equivalent, unitary_equivalent

SEEDS = range(40)


@pytest.mark.parametrize("policy", ["rewrite", "flush"])
@pytest.mark.parametrize("seed", SEEDS)
def test_random_clifford_t_unitary(seed, policy):
    n = 2 + seed % 6
    qc = circuits.random_layered(n, depth=12, seed=seed, t_prob=0.2)
    assert unitary_equivalent(qc, compress(qc, policy))


@pytest.mark.parametrize("seed", SEEDS)
def test_random_rz_unitary(seed):
    qc = circuits.random_layered(4, depth=8, seed=1000 + seed)
    qc.rz(0.37, 1)
    qc.p(1.1, 2)
    qc.x(3)
    qc.rx(0.5, 3)  # not diagonal: forces a flush
    qc.cp(0.4, 0, 3)  # two-qubit non-Clifford: flush
    qc.h(0)
    assert unitary_equivalent(qc, compress(qc))


@pytest.mark.parametrize("seed", range(10))
def test_large_clifford_tableau(seed):
    qc = circuits.random_layered(64, depth=40, seed=seed)
    res = compress(qc)
    assert res.removed  # something was absorbed
    assert tableau_equivalent(qc, res)


@pytest.mark.parametrize("seed", range(10))
def test_synthesized_clifford_tableau(seed):
    qc = circuits.random_clifford_synth(20, seed)
    assert tableau_equivalent(qc, compress(qc))


@pytest.mark.parametrize("seed", range(20))
def test_measurements_lockstep(seed):
    n = 3 + seed % 30
    qc = circuits.random_layered(n, depth=20, seed=seed, measure=True)
    assert measurements_equivalent(qc, compress(qc), seed=seed)


def test_mid_circuit_measure_and_reset():
    qc = circuits.repetition_syndrome(d=5, rounds=3)
    # inject absorbable gates between rounds
    qc2 = QuantumCircuit(qc.num_qubits, qc.num_clbits)
    for inst in qc.data:
        qc2.append(inst)
        if inst.operation.name == "reset":
            qc2.x(0)
            qc2.swap(1, 2)
    res = compress(qc2)
    assert res.removed["x"] and res.removed["swap"]
    assert measurements_equivalent(qc2, res)


def test_lockstep_detects_wrong_flip():
    qc = QuantumCircuit(1, 1)
    qc.x(0)
    qc.measure(0, 0)
    res = compress(qc)
    assert res.measurements == [(0, 0, True)]
    assert measurements_equivalent(qc, res)
    res.measurements[0] = (0, 0, False)
    assert not measurements_equivalent(qc, res)


def test_t_rewrite_emits_tdg():
    qc = QuantumCircuit(1)
    qc.x(0)
    qc.t(0)
    res = compress(qc, "rewrite")
    assert [i.operation.name for i in res.circuit.data] == ["tdg"]
    assert res.table.pending(0) == "X"
    assert unitary_equivalent(qc, res)
    res = compress(qc, "flush")
    assert [i.operation.name for i in res.circuit.data] == ["x", "t"]
    assert res.flushed == 1


def test_x_through_h_becomes_z():
    qc = QuantumCircuit(1)
    qc.x(0)
    qc.h(0)
    res = compress(qc)
    assert res.table.pending(0) == "Z"


@pytest.mark.parametrize("n", [2, 3, 4, 5, 8, 9])
def test_qft_removes_floor_half_swaps(n):
    qc = circuits.qft(n)
    res = compress(qc)
    assert res.removed == {"swap": n // 2}
    if n <= 8:
        assert unitary_equivalent(qc, res)


def test_known_algorithms_preserve_counts():
    from qiskit_aer import AerSimulator

    from qarch.verify import counts_with_flips

    sim = AerSimulator(seed_simulator=7)
    for qc in [circuits.bernstein_vazirani(6, seed=3), circuits.grover(4, marked=0b1010),
               circuits.grover(5, marked=0b10110, decompose_mcz=True)]:
        res = compress(qc)
        a = sim.run(transpile(qc, sim), shots=2000).result().get_counts()
        b = sim.run(transpile(res.circuit, sim), shots=2000).result().get_counts()
        b = counts_with_flips(b, res, qc.num_clbits)
        top_a = max(a, key=a.get)
        assert top_a == max(b, key=b.get)
        assert b[top_a] > 1500


def test_pass_manager_wrapper():
    qc = circuits.random_layered(5, depth=10, seed=1)
    pm = PassManager([AddressTablePass()])
    out = pm.run(qc)
    assert out.size() <= qc.size()
    assert pm.property_set["address_table"] is not None


def test_diagonal_gate_keeps_pending_z():
    from qiskit.circuit.library import MCMTGate, ZGate

    qc = QuantumCircuit(3)
    qc.z(0)
    qc.x(1)
    qc.cp(0.3, 0, 1)
    qc.append(MCMTGate(ZGate(), 2, 1), [0, 1, 2])
    qc.rzz(0.7, 1, 2)
    res = compress(qc)
    assert res.flushed == 1            # only the X on line 1; the Z stays pending
    assert res.table.pending(0) == "Z"
    assert unitary_equivalent(qc, res)


def test_exact_toffoli_chain():
    from qiskit.quantum_info import Statevector

    # The chain needs clean ancillas, so compare on every input with the ancillas in |0>.
    for k in (2, 3, 4, 5):
        nq = 2 * k - 1
        controls, target, anc = list(range(k)), k, list(range(k + 1, nq))
        qc = QuantumCircuit(nq)
        circuits.mcx_vchain(qc, controls, target, anc)
        ref = QuantumCircuit(nq)
        ref.mcx(controls, target)
        for b in range(2 ** (k + 1)):
            psi = Statevector.from_int(b, 2 ** nq)
            assert psi.evolve(qc).equiv(psi.evolve(ref)), (k, b)


@pytest.mark.parametrize("pattern", ["a", "b", "ab"])
def test_cp_and_rzz_rewrite_under_pending_x(pattern):
    qc = QuantumCircuit(2)
    if "a" in pattern:
        qc.x(0)
    if "b" in pattern:
        qc.x(1)
    qc.cp(0.7, 0, 1)
    qc.rzz(0.4, 0, 1)
    qc.h(0)
    qc.cp(1.3, 1, 0)
    res = compress(qc, "rewrite")
    assert res.rewritten_2q >= 2
    assert res.flushed == 0
    assert unitary_equivalent(qc, res)
    res_f = compress(qc, "flush")
    assert res_f.rewritten_2q == 0 and res_f.flushed >= 1
    assert unitary_equivalent(qc, res_f)
    res_n = compress(qc, "rewrite", rewrite_cp=False)
    assert res_n.rewritten_2q == 0 and res_n.flushed >= 1
    assert unitary_equivalent(qc, res_n)
