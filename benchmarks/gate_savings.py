"""Step 3 of the research plan: gate savings from the address-table pass, per circuit family,
against Qiskit's optimization_level=3 baseline.

    uv run python benchmarks/gate_savings.py            # full sweep, n = 4..64
    uv run python benchmarks/gate_savings.py --quick    # n = 4..16, fewer seeds
    uv run python benchmarks/gate_savings.py --coupling line

Writes results/gate_savings.csv (one row per family, n, seed, treatment),
results/gate_savings_regression.csv (slope and R^2 of gates saved vs n) and
results/gate_savings.png.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qarch import circuits  # noqa: E402
from qarch.metrics import TREATMENTS, evaluate, rows  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"

# family -> (builder(n, seed), needs seeds?)
# Fixed depth so the gate count, and hence the saving, is linear in n (depth = 2n made it quadratic).
LAYER_DEPTH = 32

CLIFFORD_FAMILIES = {
    "ghz": (lambda n, s: circuits.ghz(n), False),
    "bernstein_vazirani": (lambda n, s: circuits.bernstein_vazirani(n - 1, s), True),
    "repetition_syndrome": (lambda n, s: circuits.repetition_syndrome((n + 1) // 2, rounds=3), False),
    "random_layered": (lambda n, s: circuits.random_layered(n, depth=LAYER_DEPTH, seed=s), True),
    "random_clifford_synth": (lambda n, s: circuits.random_clifford_synth(n, s), True),
}
BEYOND_CLIFFORD = {
    "qft": (lambda n, s: circuits.qft(n), False),
    "grover_native_mcz": (lambda n, s: circuits.grover(n, marked=s % 2 ** n, iterations=1), True),
    "grover_clifford_t": (lambda n, s: circuits.grover(n, marked=s % 2 ** n, iterations=1,
                                                       decompose_mcz=True), True),
}


def regress(ns, saved):
    ns = np.asarray(ns, float)
    saved = np.asarray(saved, float)
    if len(set(ns)) < 2:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(ns, saved, 1)
    pred = slope * ns + intercept
    ss_res = ((saved - pred) ** 2).sum()
    ss_tot = ((saved - saved.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot else 1.0
    return float(slope), float(r2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--coupling", choices=["all", "line"], default="all")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--policy", choices=["rewrite", "flush"], default="rewrite")
    args = ap.parse_args()

    ns_clifford = list(range(4, 17, 4)) if args.quick else list(range(4, 65, 4))
    # Native MCZ cannot be synthesised by O3 at large n without ancillas, so it stops at 10;
    # the Clifford+T Grover (Toffoli v-chain) runs the full range with a single iteration.
    ns_grover_native = [4, 6] if args.quick else [4, 5, 6, 7, 8, 9, 10]
    ns_grover_ct = ns_clifford
    ns_shor = [8] if args.quick else [8, 12, 16]
    seeds = list(range(1 if args.quick else args.seeds))

    RESULTS.mkdir(exist_ok=True)
    all_rows: list[dict] = []
    t0 = time.time()
    plan = [(f, b, s, ns_clifford) for f, (b, s) in CLIFFORD_FAMILIES.items()]
    plan += [("qft", *BEYOND_CLIFFORD["qft"], ns_clifford)]
    plan += [("grover_native_mcz", *BEYOND_CLIFFORD["grover_native_mcz"], ns_grover_native)]
    plan += [("grover_clifford_t", *BEYOND_CLIFFORD["grover_clifford_t"], ns_grover_ct)]
    plan += [("shor15", lambda n, s: circuits.shor15(a=7, n_count=n), False, ns_shor)]

    # Record the inputs that decide the numbers next to the numbers (item 5 of the review).
    import json
    import subprocess
    try:
        commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                         cwd=RESULTS.parent, text=True).strip()
    except Exception:
        commit = "unknown"
    params = {"commit": commit, "coupling": args.coupling, "policy": args.policy,
              "seeds": seeds, "quick": args.quick, "basis": __import__("qarch.metrics").metrics.BASIS,
              "random_layered": {"probs": circuits.DEFAULT_LAYER_PROBS, "depth": LAYER_DEPTH,
                                 "note": "fixed depth so gates grow linearly in n; saved/gates is the scale-free statistic"},
              "grover": {"iterations": 1, "clifford_t_ancillas": "max(0, n-3) clean",
                         "note": "qubit counts of grover_clifford_t are not comparable to grover_native_mcz"},
              "shor15": "compiled order-aware N=15 circuit (exponents reduced mod ord(a)); not general Shor, "
                        "cf. Smolin, Smith & Vargo 2013",
              "notes": ["o3 runs ElidePermutations, so SWAP savings overlap with the baseline",
                        "o3+pass on coupling=line assumes free relabeling; pass+o3 is routed",
                        "control flow (feed-forward) is not supported by the pass"]}
    (RESULTS / f"gate_savings_{args.coupling}_params.json").write_text(json.dumps(params, indent=1))

    for family, build, seeded, ns in plan:
        for n in ns:
            for seed in (seeds if seeded else [0]):
                qc = build(n, seed)
                counts = evaluate(qc, coupling=args.coupling, seed=seed, policy=args.policy)
                all_rows += rows(family, n, counts, seed=seed, coupling=args.coupling)
            print(f"{family:>22} n={n:<3} "
                  + " ".join(f"{t}={counts[t].gates:<5}" for t in TREATMENTS)
                  + f"  ({time.time() - t0:5.1f}s)")

    fields = list(all_rows[0].keys())
    with open(RESULTS / f"gate_savings_{args.coupling}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)

    # regression of gates saved vs n, per family, for pass-vs-raw and o3+pass-vs-o3
    reg_rows = []
    by = {}
    for r in all_rows:
        by.setdefault((r["family"], r["n"], r["seed"]), {})[r["treatment"]] = r
    summary = {}
    for (fam, n, seed), tr in by.items():
        for label, a, b in (("pass_vs_raw", "pass", "raw"), ("o3pass_vs_o3", "o3+pass", "o3"),
                            ("passo3_vs_o3", "pass+o3", "o3")):
            for metric in ("gates", "gates_vz", "two_qubit", "cx_equiv", "depth"):
                summary.setdefault((fam, label, metric), []).append((n, tr[b][metric] - tr[a][metric]))
    for (fam, label, metric), pts in sorted(summary.items()):
        ns, saved = zip(*pts)
        slope, r2 = regress(ns, saved)
        reg_rows.append({"family": fam, "comparison": label, "metric": metric,
                         "slope": round(slope, 4), "r2": round(r2, 4),
                         "mean_saved": round(float(np.mean(saved)), 2)})
    with open(RESULTS / f"gate_savings_{args.coupling}_regression.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(reg_rows[0].keys()))
        w.writeheader()
        w.writerows(reg_rows)

    print("\ngates saved vs n (slope, R^2):")
    for r in reg_rows:
        if r["metric"] == "gates":
            print(f"  {r['family']:>22} {r['comparison']:<14} slope={r['slope']:7.3f} "
                  f"R2={r['r2']:6.3f} mean={r['mean_saved']:7.1f}")

    plot(all_rows, args.coupling)


def plot(all_rows, coupling):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fams = sorted({r["family"] for r in all_rows})
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    for ax, (a, b, title) in zip(axes, (("pass", "raw", "pass vs raw"),
                                        ("o3+pass", "o3", "O3 + pass vs O3"),
                                        ("pass+o3", "o3", "pass + O3 vs O3"))):
        for fam in fams:
            pts = {}
            for r in all_rows:
                if r["family"] == fam:
                    pts.setdefault(r["n"], {}).setdefault(r["treatment"], []).append(r["gates"])
            ns = sorted(pts)
            saved = [np.mean(pts[n][b]) - np.mean(pts[n][a]) for n in ns]
            ax.plot(ns, saved, marker="o", label=fam)
        ax.set_xlabel("qubits n")
        ax.set_ylabel("gates saved")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.suptitle(f"coupling = {coupling}; shor15 is a compiled order-aware N=15 circuit, not general Shor; "
                 "grover_clifford_t uses n-3 ancillas", fontsize=8)
    fig.savefig(RESULTS / f"gate_savings_{coupling}.png", dpi=150)
    print(f"\nwrote {RESULTS}/gate_savings_{coupling}[.csv|_regression.csv|_params.json|.png]")


if __name__ == "__main__":
    main()
