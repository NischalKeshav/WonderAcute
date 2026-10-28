"""Step 6 of the research plan: success probability under noise, original vs compressed,
across the number of operations packed into one clock tick.

    uv run python benchmarks/noise.py [--quick]

Two fair pairs are simulated with the same seeds and shots: raw vs pass, and Qiskit O3 vs
pass followed by O3 (the step-3 baseline). ``--p-idle-sweep`` sweeps the idle rate instead of
the ops-per-tick bound.

Circuits: Grover (n = 3, 4, 5, one marked item, full iteration count, MCZ as Clifford+T) and
Bernstein–Vazirani (n = 6, 10). Noise: amplitude damping + dephasing after every gate and on
every idle qubit per tick (qarch.noise). Compressed circuits get their classical bit flips
applied before scoring. >= 1000 shots per configuration, several seeds, 95% CI = 1.96 x max(across-seed
standard error, binomial error of the pooled shots). Writes results/noise.csv and results/noise.png.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qarch import circuits  # noqa: E402
from qarch.metrics import baseline  # noqa: E402
from qarch.noise import ONE_Q, TWO_Q, noise_model, with_idle  # noqa: E402
from qarch.passes import compress  # noqa: E402
from qarch.verify import counts_with_flips  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"


def success(counts: dict, target: str, n_clbits: int) -> float:
    total = sum(counts.values())
    return sum(v for k, v in counts.items() if k.replace(" ", "")[-n_clbits:] == target) / total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--p1", type=float, default=1e-3)
    ap.add_argument("--p2", type=float, default=1e-2)
    ap.add_argument("--p-idle", type=float, default=2e-3)
    ap.add_argument("--p-idle-sweep", action="store_true",
                    help="sweep the idle rate at unbounded ops per tick instead of sweeping ticks")
    args = ap.parse_args()
    shots = 1000
    seeds = range(3 if args.quick else 10)
    ticks = [1, 2, 4, None] if args.quick else [1, 2, 3, 4, 6, 8, None]
    idle_rates = [args.p_idle]
    if args.p_idle_sweep:
        ticks = [None]
        idle_rates = [0.0, 5e-4, 1e-3, 2e-3, 4e-3, 8e-3]
    tag = "_idle_sweep" if args.p_idle_sweep else ""

    cases = []
    for n in ([3, 4] if args.quick else [3, 4, 5]):
        marked = (0b1011010 >> (7 - n)) & (2 ** n - 1)
        qc = circuits.grover(n, marked, decompose_mcz=True)
        cases.append((f"grover_{n}", qc, format(marked, f"0{n}b")))
    for n in ([6] if args.quick else [6, 10]):
        qc = circuits.bernstein_vazirani(n, seed=n)
        cases.append((f"bv_{n}", qc, qc.metadata["secret"]))

    basis = ONE_Q + TWO_Q + ["id", "measure"]
    rows = []
    RESULTS.mkdir(exist_ok=True)
    for name, qc, target in cases:
        res = compress(qc)
        # two fair pairs: raw vs pass, and Qiskit O3 vs pass followed by O3 (the step-3 baseline)
        variants = {"original": (transpile(qc, basis_gates=basis, optimization_level=0), None),
                    "compressed": (transpile(res.circuit, basis_gates=basis, optimization_level=0), res),
                    "o3": (baseline(qc, seed=0), None),
                    "pass+o3": (baseline(res.circuit, seed=0), res)}
        for variant, (circ, cres) in variants.items():
          for p_idle in idle_rates:
            nm = noise_model(args.p1, args.p2, p_idle)
            for opt in ticks:
                sched = with_idle(circ, opt)
                n_ticks = sched.depth(lambda i: i.operation.name != "measure")
                ps = []
                for seed in seeds:
                    sim = AerSimulator(noise_model=nm, seed_simulator=seed)
                    counts = sim.run(sched, shots=shots).result().get_counts()
                    if cres is not None:
                        counts = counts_with_flips(counts, cres, qc.num_clbits)
                    ps.append(success(counts, target, qc.num_clbits))
                ps = np.array(ps)
                # across-seed spread, floored at the binomial error of the pooled shots
                m = float(ps.mean())
                sem = ps.std(ddof=1) / np.sqrt(len(ps)) if len(ps) > 1 else 0.0
                half = 1.96 * max(sem, np.sqrt(m * (1 - m) / (shots * len(ps))))
                rows.append({"circuit": name, "variant": variant,
                             "ops_per_tick": opt if opt is not None else "unbounded",
                             "ticks": n_ticks, "gates": circ.size() - qc.num_clbits,
                             "p_success": round(float(ps.mean()), 4), "ci95": round(float(half), 4),
                             "shots": shots, "seeds": len(ps),
                             "p1": args.p1, "p2": args.p2, "p_idle": p_idle})
                print(f"{name:>10} {variant:>10} ops/tick={str(opt):>4} p_idle={p_idle:<6} "
                      f"ticks={n_ticks:<5} p={ps.mean():.3f} ±{half:.3f}")

    with open(RESULTS / f"noise{tag}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = sorted({r["circuit"] for r in rows})
    fig, axes = plt.subplots(1, len(names), figsize=(4 * len(names), 4), sharey=True)
    for ax, name in zip(np.atleast_1d(axes), names):
        for variant, m in (("original", "o"), ("compressed", "s"), ("o3", "^"), ("pass+o3", "D")):
            pts = [r for r in rows if r["circuit"] == name and r["variant"] == variant]
            if args.p_idle_sweep:
                xs = [r["p_idle"] for r in pts]
            else:
                xs = [len(ticks) - 1 if r["ops_per_tick"] == "unbounded" else ticks.index(r["ops_per_tick"])
                      for r in pts]
            ax.errorbar(xs, [r["p_success"] for r in pts], yerr=[r["ci95"] for r in pts],
                        marker=m, capsize=3, label=variant)
        if args.p_idle_sweep:
            ax.set_xlabel("idle error per tick")
        else:
            ax.set_xticks(range(len(ticks)))
            ax.set_xticklabels([str(t) if t else "all" for t in ticks])
            ax.set_xlabel("ops per tick")
        ax.set_title(name)
        ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("success probability")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle(f"p1={args.p1} p2={args.p2} p_idle={args.p_idle}, {shots} shots x {len(seeds)} seeds",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS / f"noise{tag}.png", dpi=150)
    print(f"wrote {RESULTS}/noise{tag}.csv and .png")


if __name__ == "__main__":
    main()
