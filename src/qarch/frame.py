"""The classical controller's address table.

The table holds two things:

* a permutation ``phys[q]``: logical line ``q`` of the source program currently lives on
  physical line ``phys[q]`` of the compressed circuit (this is how SWAP is absorbed);
* a Pauli frame ``x[p], z[p]`` on physical lines: a Pauli X^x Z^z that is pending on line
  ``p`` and has not been applied as a gate (this is how X, Y and Z are absorbed).

Invariant kept by the compiler pass: running the source prefix is equivalent, up to global
phase, to running the compressed prefix, then applying the pending frame, then routing
physical line ``phys[q]`` back to logical line ``q``.

Moving a pending Pauli P past a Clifford gate G uses G P = (G P G^dagger) G, so the frame is
replaced by G P G^dagger. These are the rows of proofs/propagation-tables.json with the signs
dropped (a sign is a global phase).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from qiskit.circuit import Instruction
from qiskit.quantum_info import Clifford, Pauli

# Gates absorbed into the table without emitting anything.
PAULI_GATES = {"x", "y", "z", "id"}
PERMUTATION_GATES = {"swap"}

# Clifford gates with hand-written bit rules (checked against Qiskit and the proofs tables in tests).
FAST_CLIFFORDS = {"h", "s", "sdg", "sx", "sxdg", "cx", "cz", "cy"}

# Diagonal non-Clifford gates whose conjugation by X is the same gate with a negated angle,
# up to global phase: X T X ~ Tdg, X RZ(t) X = RZ(-t), X P(l) X ~ P(-l).
ANGLE_FLIP = {"t": "tdg", "tdg": "t"}
ANGLE_NEGATE = {"rz", "p", "u1"}
NON_ENTANGLING = {"measure", "reset", "barrier", "delay"}


@lru_cache(maxsize=512)
def _clifford_for(name: str, params: tuple, num_qubits: int) -> Clifford | None:
    from qiskit.circuit.library import get_standard_gate_name_mapping

    gate = get_standard_gate_name_mapping().get(name)
    if gate is None or gate.num_qubits != num_qubits:
        return None
    try:
        if params:
            gate = type(gate)(*params)
        return Clifford(gate)
    except Exception:
        return None


def clifford_of(op: Instruction) -> Clifford | None:
    """The Clifford for ``op`` if it is a parameter-bound Clifford gate, else None."""
    if op.name in NON_ENTANGLING:
        return None
    try:
        params = tuple(float(p) for p in op.params)
    except TypeError:
        return None  # unbound parameter
    cliff = _clifford_for(op.name, params, op.num_qubits)
    if cliff is None and op.num_qubits <= 4:
        try:
            cliff = Clifford(op)
        except Exception:
            return None
    return cliff


def is_clifford(op: Instruction) -> bool:
    return (op.name in FAST_CLIFFORDS or op.name in PAULI_GATES
            or op.name in PERMUTATION_GATES or clifford_of(op) is not None)


@dataclass
class AddressTable:
    num_qubits: int
    phys: list[int] = field(default_factory=list)
    x: np.ndarray | None = None
    z: np.ndarray | None = None

    def __post_init__(self):
        if not self.phys:
            self.phys = list(range(self.num_qubits))
        if self.x is None:
            self.x = np.zeros(self.num_qubits, dtype=bool)
        if self.z is None:
            self.z = np.zeros(self.num_qubits, dtype=bool)

    def copy(self) -> AddressTable:
        return AddressTable(self.num_qubits, list(self.phys), self.x.copy(), self.z.copy())

    # --- absorbed operations (no quantum gate) ---------------------------------------------

    def apply_pauli(self, name: str, logical: int) -> None:
        p = self.phys[logical]
        if name in ("x", "y"):
            self.x[p] ^= True
        if name in ("z", "y"):
            self.z[p] ^= True

    def swap(self, a: int, b: int) -> None:
        self.phys[a], self.phys[b] = self.phys[b], self.phys[a]

    # --- frame propagation ------------------------------------------------------------------

    def conjugate_through(self, op: Instruction, physical: list[int]) -> None:
        """Replace the frame F on ``physical`` by G F G^dagger, for a Clifford gate G = ``op``."""
        x, z = self.x, self.z
        name = op.name
        if name == "h":
            (p,) = physical
            x[p], z[p] = z[p], x[p]
        elif name in ("s", "sdg"):
            (p,) = physical
            z[p] ^= x[p]
        elif name in ("sx", "sxdg"):
            (p,) = physical
            x[p] ^= z[p]
        elif name == "cx":
            c, t = physical
            x[t] ^= x[c]
            z[c] ^= z[t]
        elif name == "cz":
            a, b = physical
            z[a], z[b] = z[a] ^ x[b], z[b] ^ x[a]
        elif name == "cy":
            c, t = physical
            z[c] ^= x[t] ^ z[t]
            x[t] ^= x[c]
            z[t] ^= x[c]
        elif name in PAULI_GATES:
            pass  # Paulis commute with Paulis up to sign
        else:
            cliff = clifford_of(op)
            if cliff is None:
                raise ValueError(f"{name} is not a Clifford gate")
            # Qiskit Pauli: qubit i of the local Pauli is physical[i].
            local = Pauli((z[physical].copy(), x[physical].copy()))
            out = local.evolve(cliff, frame="s")  # C P C^dagger
            x[physical] = out.x
            z[physical] = out.z

    # --- queries ----------------------------------------------------------------------------

    def pending(self, p: int) -> str:
        """Pending Pauli on physical line p: 'I', 'X', 'Y' or 'Z'."""
        return "IZXY"[2 * int(self.x[p]) + int(self.z[p])]

    def clear(self, p: int) -> None:
        self.x[p] = False
        self.z[p] = False

    def is_trivial(self) -> bool:
        return not self.x.any() and not self.z.any() and self.phys == list(range(self.num_qubits))

    def measurement_flip(self, p: int) -> bool:
        """A Z-basis measurement of physical line p reads the opposite bit iff X or Y is pending."""
        return bool(self.x[p])

    def __str__(self) -> str:
        rows = [f"  q{q} -> p{p}  pending {self.pending(p)}" for q, p in enumerate(self.phys)]
        return "AddressTable(\n" + "\n".join(rows) + "\n)"
