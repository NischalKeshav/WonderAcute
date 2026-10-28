"""Step 5 + 6 together: does a shared oracle still pay once the waiting instances decohere?

    uv run python benchmarks/shared_oracle_noise.py [--quick]

k Grover instances (k = 1..3) search the same N = 8 table through a unary-iteration table
lookup used as a phase oracle (output qubit in |->). ``duplicated``: every instance has its
own lookup ancillas and output qubit, so all k run in parallel. ``shared``: one output qubit
and one ancilla pair serve all instances; the data dependency serialises the lookups, so
instances wait in a queue while the others are served, and the idle channel acts on them
every tick they wait. Everything is lowered to Clifford+T so the step-6 noise model applies.

Qubits: duplicated 6k, shared 3k + 3 (N = 8: 3 address, 1 output, 2 ancillas). Success is the
probability that an instance reads the marked item (mean over instances), with 95% CIs across
seeds. The idle rate is swept because it is what decides the answer. Writes
results/shared_oracle_noise.csv and .png.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qarch import qram  # noqa: E402
from qarch.circuits import toffoli  # noqa: E402
from qarch.noise import ONE_Q, TWO_Q, noise_model, with_idle  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
N = 8
NBITS = 3


def multi_grover(k: int, marked: int, shared: bool) -> QuantumCircuit:
    data = [1 if i == marked else 0 for i in range(N)]
    lookup = qram.table_lookup(data, 1, exact=True).circuit     # qubits: a0 a1 a2 | out | anc0 anc1
    n_oracle = lookup.num_qubits - NBITS
    nq = k * NBITS + (n_oracle if shared else k * n_oracle)
    qc = QuantumCircuit(nq, k * NBITS)
    iters = max(1, math.floor(math.pi / 4 * math.sqrt(N)))

    def regs(i):
        s = list(range(i * NBITS, (i + 1) * NBITS))
        base = k * NBITS + (0 if shared else i * n_oracle)
        return s, list(range(base, base + n_oracle))

    for i in range(k):
        s, o = regs(i)
        qc.h(s)
        if not shared or i == 0:
            qc.x(o[0]); qc.h(o[0])          # output qubit in |->, phase kickback
    for _ in range(iters):
        for i in range(k):
            s, o = regs(i)
            qc.compose(lookup, s + o, inplace=True)
            qc.h(s); qc.x(s)
            qc.h(s[2]); toffoli(qc, s[0], s[1], s[2]); qc.h(s[2])
            qc.x(s); qc.h(s)
    for i in range(k):
        s, _ = regs(i)
        qc.measure(s, s)
    return qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--p1", type=float, default=1e-3)
    ap.add_argument("--p2", type=float, default=1e-2)
    args = ap.parse_args()
    shots = 1000
    seeds = range(3 if args.quick else 10)
    idle_rates = [0.0, 1e-3, 3e-3] if args.quick else [0.0, 5e-4, 1e-3, 2e-3, 4e-3, 8e-3]
    ks = [1, 2, 3]
    marked = 5
    basis = ONE_Q + TWO_Q + ["id", "measure"]
    RESULTS.mkdir(exist_ok=True)
    rows = []
    for k in ks:
        for design in ("duplicated", "shared"):
            qc = multi_grover(k, marked, design == "shared")
            circ = with_idle(transpile(qc, basis_gates=basis, optimization_level=0), None)
            ticks = circ.depth(lambda i: i.operation.name != "measure")
            for p_idle in idle_rates:
                nm = noise_model(args.p1, args.p2, p_idle)
                per_seed = []
                for seed in seeds:
                    sim = AerSimulator(noise_model=nm, seed_simulator=seed)
                    counts = sim.run(circ, shots=shots).result().get_counts()
                    hits = np.zeros(k)
                    for key, v in counts.items():
                        bits = key.replace(" ", "")[::-1]        # little-endian clbit order
                        for i in range(k):
                            word = bits[i * NBITS:(i + 1) * NBITS][::-1]
                            hits[i] += v if int(word, 2) == marked else 0
                    per_seed.append(hits.mean() / shots)
                ps = np.array(per_seed)
                # across-seed spread, floored at the binomial error of the pooled shots
                m = float(ps.mean())
                sem = ps.std(ddof=1) / np.sqrt(len(ps)) if len(ps) > 1 else 0.0
                half = 1.96 * max(sem, np.sqrt(m * (1 - m) / (shots * len(ps))))
                rows.append({"k": k, "design": design, "qubits": qc.num_qubits, "ticks": ticks,
                             "p_idle": p_idle, "p_success": round(float(ps.mean()), 4),
                             "ci95": round(float(half), 4), "shots": shots, "seeds": len(ps),
                             "p1": args.p1, "p2": args.p2})
                print(f"k={k} {design:>10} qubits={qc.num_qubits:<3} ticks={ticks:<4} "
                      f"p_idle={p_idle:<7} p={ps.mean():.3f} ±{half:.3f}")

    with open(RESULTS / "shared_oracle_noise.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(ks), figsize=(4.2 * len(ks), 4), sharey=True)
    for ax, k in zip(np.atleast_1d(axes), ks):
        for design, m in (("duplicated", "o"), ("shared", "s")):
            pts = [r for r in rows if r["k"] == k and r["design"] == design]
            ax.errorbar([r["p_idle"] for r in pts], [r["p_success"] for r in pts],
                        yerr=[r["ci95"] for r in pts], marker=m, capsize=3,
                        label=f"{design} ({pts[0]['qubits']} qubits, {pts[0]['ticks']} ticks)")
        ax.set_title(f"k = {k} Grover instances, N = 8")
        ax.set_xlabel("idle error per tick")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    np.atleast_1d(axes)[0].set_ylabel("per-instance success probability")
    fig.suptitle(f"shared vs duplicated table-lookup oracle under noise; p1={args.p1} p2={args.p2}, "
                 f"{shots} shots x {len(seeds)} seeds", fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS / "shared_oracle_noise.png", dpi=150)
    print(f"wrote {RESULTS / 'shared_oracle_noise.csv'} and .png")


if __name__ == "__main__":
    main()
