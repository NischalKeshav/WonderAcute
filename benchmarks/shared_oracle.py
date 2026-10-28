"""Step 5 of the research plan: k Grover instances against one shared oracle versus k copies.

    uv run python benchmarks/shared_oracle.py

Model (analytic; oracle qubits and lookup depth are read off the *built* step-4 circuits,
so the table lookup is the 2(N-2)-Toffoli compute/uncompute version, not Babbush's N-1):
  each Grover instance searches an N-entry table, running floor(pi/4 * sqrt(N)) iterations,
  and each iteration makes one lookup (the oracle) followed by a diffusion step on log2(N)
  qubits whose depth we take as 2*log2(N) + 1 (H, X, MCZ, X, H layers).

  duplicated:  k oracles, instances run in parallel
      qubits = k * (log2 N + oracle_qubits)
      depth  = iterations * (lookup_depth + diffusion_depth)
  shared:      one oracle, lookups queued one after another (a bucket-brigade or table-lookup
               memory answers one query at a time, Xu et al. 2025)
      qubits = k * log2 N + oracle_qubits
      depth  = iterations * (k * lookup_depth + diffusion_depth)

The whitepaper claims 66% fewer qubits for 3 instances against a 1,024-line oracle; this
script checks that number and shows what it costs in depth. Writes results/shared_oracle.csv
and results/shared_oracle.png.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qarch import qram  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
NS = [8, 64, 256, 1024]
KS = list(range(1, 9))


def lookup_cost(design: str, N: int) -> tuple[int, int]:
    """(oracle qubits beyond the address register, lookup depth) from built circuits."""
    data = [i % 2 for i in range(N)]
    lk = qram.bucket_brigade(data) if design == "bucket_brigade" else qram.table_lookup(data)
    k = int(math.log2(N))
    return lk.circuit.num_qubits - k, lk.circuit.depth()


def main():
    RESULTS.mkdir(exist_ok=True)
    rows = []
    for design in ("bucket_brigade", "table_lookup"):
        for N in NS:
            n = int(math.log2(N))
            oracle_q, lookup_d = lookup_cost(design, N)
            iters = max(1, math.floor(math.pi / 4 * math.sqrt(N)))
            diffusion_d = 2 * n + 1
            for k in KS:
                dup_q = k * (n + oracle_q)
                dup_d = iters * (lookup_d + diffusion_d)
                sh_q = k * n + oracle_q
                sh_d = iters * (k * lookup_d + diffusion_d)
                rows.append({"design": design, "N": N, "k": k, "iterations": iters,
                             "oracle_qubits": oracle_q, "lookup_depth": lookup_d,
                             "dup_qubits": dup_q, "dup_depth": dup_d,
                             "shared_qubits": sh_q, "shared_depth": sh_d,
                             "qubits_saved_frac": round(1 - sh_q / dup_q, 4),
                             "depth_ratio": round(sh_d / dup_d, 4),
                             "qubit_depth_dup": dup_q * dup_d, "qubit_depth_shared": sh_q * sh_d})
    with open(RESULTS / "shared_oracle.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        if r["N"] == 1024 and r["k"] == 3:
            print(f"{r['design']:>15} N=1024 k=3: qubits {r['dup_qubits']} -> {r['shared_qubits']} "
                  f"({100 * r['qubits_saved_frac']:.1f}% fewer), depth x{r['depth_ratio']:.2f}, "
                  f"qubit*depth x{r['qubit_depth_shared'] / r['qubit_depth_dup']:.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for design, ls in (("bucket_brigade", "-"), ("table_lookup", "--")):
        for N in NS:
            pts = [r for r in rows if r["design"] == design and r["N"] == N]
            axes[0].plot([r["k"] for r in pts], [100 * r["qubits_saved_frac"] for r in pts], ls,
                         marker="o", label=f"{design} N={N}")
            axes[1].plot([r["k"] for r in pts], [r["depth_ratio"] for r in pts], ls, marker="o")
    axes[0].set_ylabel("qubits saved by sharing (%)")
    axes[1].set_ylabel("depth, shared / duplicated")
    for ax in axes:
        ax.set_xlabel("Grover instances k")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle("One queued oracle vs k copies: qubits saved against depth paid", fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS / "shared_oracle.png", dpi=150)
    print(f"wrote {RESULTS / 'shared_oracle.csv'} and shared_oracle.png")


if __name__ == "__main__":
    main()
