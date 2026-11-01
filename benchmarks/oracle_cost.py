"""Step 4 of the research plan: what one memory lookup costs.

    uv run python benchmarks/oracle_cost.py

For N = 8..1024 entries (d = 1 data bit, and d = 8 for the closed forms) records, per lookup,
qubits, total gates, Toffoli/CSWAP count, T count and depth for

  bucket_brigade   built circuit (native CSWAP) and its closed form
  table_lookup     built circuit (native Toffoli) and its closed form, plus Babbush's N-1 figure
  select_swap      closed form only

and, for the built circuits, how many X gates the address-table pass absorbs and how many it
has to flush again before a CSWAP/Toffoli (which is not diagonal). Writes
results/oracle_cost.csv and results/oracle_cost.png. The whitepaper treats this cost as zero.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qarch import qram  # noqa: E402
from qarch.passes import compress  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
NS = [8, 16, 32, 64, 128, 256, 512, 1024]
BUILD_UP_TO = 1024      # bucket brigade at N=1024 is ~4k qubits / ~8k CSWAPs: fine to build


def measured(lk: qram.Lookup) -> dict:
    qc = lk.circuit
    names = [i.operation.name for i in qc.data]
    toff = names.count("ccx") + names.count("cswap")
    res = compress(qc)
    return {"qubits": qc.num_qubits, "gates": len(names), "toffoli": toff, "t_count": 7 * toff,
            "cnot": names.count("cx"), "x": names.count("x"), "depth": qc.depth(),
            "x_absorbed": res.removed.get("x", 0), "x_flushed": res.flushed,
            "gates_after_pass": res.circuit.size()}


def main():
    RESULTS.mkdir(exist_ok=True)
    rows = []
    t0 = time.time()
    rng = np.random.default_rng(0)
    for N in NS:
        data1 = [int(v) for v in rng.integers(0, 2, N)]
        for design, build, closed in (("bucket_brigade", qram.bucket_brigade, qram.bucket_brigade_counts),
                                      ("table_lookup", qram.table_lookup, qram.table_lookup_counts)):
            row = {"design": design, "N": N, "d": 1, "source": "closed_form", **closed(N, 1)}
            rows.append(row)
            if N <= BUILD_UP_TO:
                m = measured(build(data1, 1))
                rows.append({"design": design, "N": N, "d": 1, "source": "built", **m})
                print(f"{design:>15} N={N:<5} qubits={m['qubits']:<6} toffoli={m['toffoli']:<6} "
                      f"depth={m['depth']:<6} x_absorbed={m['x_absorbed']:<5} "
                      f"x_flushed={m['x_flushed']:<5} ({time.time() - t0:5.1f}s)")
            rows.append({"design": design, "N": N, "d": 8, "source": "closed_form", **closed(N, 8)})
        for d in (1, 8):
            rows.append({"design": "select_swap", "N": N, "d": d, "source": "closed_form",
                         **qram.select_swap_counts(N, d)})

    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("design", "N", "d", "source"), k))
    with open(RESULTS / "oracle_cost.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, metric in zip(axes, ("qubits", "toffoli", "depth")):
        for design in ("bucket_brigade", "table_lookup", "select_swap"):
            for source, style in (("built", "-o"), ("closed_form", "--")):
                pts = [(r["N"], r[metric]) for r in rows
                       if r["design"] == design and r["source"] == source and r["d"] == 1
                       and metric in r]
                if pts:
                    ax.plot(*zip(*pts), style, label=f"{design} ({source})")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.set_xlabel("entries N")
        ax.set_ylabel(f"{metric} per lookup (d = 1)")
        ax.grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=7)
    fig.suptitle("Cost of one oracle lookup; the whitepaper counted this as zero", fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS / "oracle_cost.png", dpi=150)
    print(f"wrote {RESULTS / 'oracle_cost.csv'} and oracle_cost.png")


if __name__ == "__main__":
    main()
