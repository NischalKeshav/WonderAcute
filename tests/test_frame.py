"""Frame propagation rules, checked against the proofs agent's tables and exact operators."""

import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from qiskit.circuit.library import get_standard_gate_name_mapping
from qiskit.quantum_info import Operator, Pauli

from qarch.frame import FAST_CLIFFORDS, AddressTable

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "proofs" / "propagation-tables.json"
TABLE_NAME = {"H": "h", "S": "s", "Sdg": "sdg", "CNOT": "cx", "CZ": "cz", "SWAP": "swap"}
GATES = get_standard_gate_name_mapping()


def set_frame(table: AddressTable, letters: str) -> None:
    for i, ch in enumerate(letters):
        table.x[i] = ch in "XY"
        table.z[i] = ch in "ZY"


def read_frame(table: AddressTable, k: int) -> str:
    return "".join(table.pending(i) for i in range(k))


def proof_rows():
    if not TABLES.exists():
        return []
    data = json.loads(TABLES.read_text())
    return [(g, row) for g, spec in data["gates"].items() for row in spec["rows"]]


@pytest.mark.skipif(not TABLES.exists(), reason="proofs/propagation-tables.json not present")
@pytest.mark.parametrize("gate,row", proof_rows(),
                         ids=lambda v: v if isinstance(v, str) else v["in"])
def test_matches_proof_tables(gate, row):
    op = GATES[TABLE_NAME[gate]]
    k = op.num_qubits
    t = AddressTable(k)
    set_frame(t, row["in"])
    t.conjugate_through(op, list(range(k)))
    assert read_frame(t, k) == row["out"]


def _label(letters: str) -> str:
    # qiskit Pauli labels are little-endian; letters are in argument order (qubit 0 first)
    return letters[::-1]


@pytest.mark.parametrize("name", sorted(FAST_CLIFFORDS | {"swap", "iswap", "ecr", "sxdg"}))
def test_every_pauli_against_operator(name):
    op = GATES[name]
    k = op.num_qubits
    U = Operator(op)
    for letters in itertools.product("IXYZ", repeat=k):
        letters = "".join(letters)
        t = AddressTable(k)
        set_frame(t, letters)
        t.conjugate_through(op, list(range(k)))
        expected = Operator(Pauli(_label(letters))).compose(U.adjoint(), front=True).compose(U)
        # expected = U P U^dagger; compare with the frame's Pauli up to sign
        got = Operator(Pauli(_label(read_frame(t, k))))
        assert expected.equiv(got), (name, letters, read_frame(t, k))


def test_generic_path_matches_fast_path():
    # cx on reversed arguments exercises the qubit-order convention of the generic path
    from qiskit.circuit.library import CXGate
    from qiskit.quantum_info import Clifford

    for letters in itertools.product("IXYZ", repeat=2):
        a = AddressTable(3)
        b = AddressTable(3)
        set_frame(a, "".join(letters) + "I")
        set_frame(b, "".join(letters) + "I")
        a.conjugate_through(CXGate(), [1, 0])
        local = Pauli((b.z[[1, 0]].copy(), b.x[[1, 0]].copy()))
        out = local.evolve(Clifford(CXGate()), frame="s")
        b.x[[1, 0]] = out.x
        b.z[[1, 0]] = out.z
        assert np.array_equal(a.x, b.x) and np.array_equal(a.z, b.z)
