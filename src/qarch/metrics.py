"""Gate-count metrics for a circuit under five treatments.

    raw        the circuit as generated
    pass       raw -> AddressTablePass
    o3         Qiskit transpile(optimization_level=3)      (the fair software baseline)
    o3+pass    o3 -> AddressTablePass
    pass+o3    AddressTablePass -> o3                       (the fair hardware comparison)

Disclosures.
* At optimization levels 2 and 3 Qiskit already runs ``ElidePermutations``, which removes SWAP
  gates into a final layout. So on the ``o3`` column SWAP savings overlap with the baseline;
  what the pass adds on top of O3 is Pauli absorption, frame propagation through Cliffords and
  the T/RZ rewrite. The ``pass`` column shows the full saving against the raw circuit.
* With ``coupling="line"`` the baseline is routed on a nearest-neighbour line. ``o3+pass`` then
  removes routing SWAPs by relabeling, so its output no longer obeys the line: it assumes the
  controller can retarget lines for free (the architecture's premise). ``pass+o3`` compresses
  first and routes afterwards, so its output is a valid line circuit.
* ``gates_vz`` leaves out diagonal one-qubit gates (rz, p, z, s, t, ...), which virtual-Z
  hardware applies as frame updates; it is the count that matters for the CP rewrite.
* ``two_qubit`` counts every gate on two or more qubits as one. ``cx_equiv`` counts a SWAP as
  three CNOTs (the convention of the research plan and the proofs' caveat 2); other multi-qubit
  gates still count as one, so for circuits with a native MCZ compare the decomposed columns.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

from .passes import compress

# SWAP kept in the basis so a SWAP the baseline leaves in place is counted as one gate, not three.
BASIS = ["h", "s", "sdg", "sx", "sxdg", "x", "y", "z", "t", "tdg", "rz", "p",
         "cx", "cz", "cp", "swap"]
TREATMENTS = ("raw", "pass", "o3", "o3+pass", "pass+o3")
# Diagonal one-qubit gates are implemented as frame updates on virtual-Z hardware (McKay 2017).
VZ_FREE = {"rz", "p", "u1", "z", "s", "sdg", "t", "tdg"}


@dataclass
class Counts:
    gates: int
    two_qubit: int
    cx_equiv: int
    depth: int
    qubits: int
    gates_vz: int = 0  # gates that are not free on virtual-Z hardware (diagonal 1q gates are free)
    absorbed: int = 0
    flushed: int = 0


def count(qc: QuantumCircuit) -> Counts:
    body = qc.remove_final_measurements(inplace=False)
    ops = [i for i in body.data if i.operation.name not in ("barrier", "measure", "delay")]
    multi = [i for i in ops if len(i.qubits) >= 2]
    return Counts(
        gates=len(ops),
        two_qubit=len(multi),
        cx_equiv=sum(3 if i.operation.name == "swap" else 1 for i in multi),
        depth=body.depth(lambda i: i.operation.name not in ("barrier", "delay")),
        qubits=qc.num_qubits,
        gates_vz=sum(1 for i in ops if i.operation.name not in VZ_FREE),
    )


def baseline(qc: QuantumCircuit, coupling: str = "all", seed: int = 0) -> QuantumCircuit:
    cmap = CouplingMap.from_line(qc.num_qubits) if coupling == "line" else None
    return transpile(qc, basis_gates=BASIS, coupling_map=cmap, optimization_level=3,
                     seed_transpiler=seed)


def evaluate(qc: QuantumCircuit, coupling: str = "all", seed: int = 0,
             policy: str = "rewrite", rewrite_cp: bool = False) -> dict[str, Counts]:
    """Counts under every treatment."""
    out = {"raw": count(qc)}
    r = compress(qc, policy, rewrite_cp)
    out["pass"] = count(r.circuit)
    out["pass"].absorbed = sum(r.removed.values())
    out["pass"].flushed = r.flushed
    o3 = baseline(qc, coupling, seed)
    out["o3"] = count(o3)
    r3 = compress(o3, policy, rewrite_cp)
    out["o3+pass"] = count(r3.circuit)
    out["o3+pass"].absorbed = sum(r3.removed.values())
    out["o3+pass"].flushed = r3.flushed
    out["pass+o3"] = count(baseline(r.circuit, coupling, seed))
    out["pass+o3"].absorbed = out["pass"].absorbed
    out["pass+o3"].flushed = out["pass"].flushed
    return out


def rows(family: str, n: int, counts: dict[str, Counts], **extra) -> list[dict]:
    return [{"family": family, "n": n, "treatment": t, **asdict(c), **extra}
            for t, c in counts.items()]
