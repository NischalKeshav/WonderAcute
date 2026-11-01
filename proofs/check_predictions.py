"""Checks the closed-form predictions in proofs/predictions.md against the Step 3-5 results.

Each prediction is computed here from scratch (no import from src/qarch) and compared with
the CSV files the code agent's benchmarks wrote. The circuits are the ones in
src/qarch/circuits.py and src/qarch/qram.py; the formulas are derived in predictions.md.

Run:  python3 proofs/check_predictions.py
The results directory is found automatically (results/ next to proofs/, or the code
worktree's results/). With no results on disk the self-consistency checks still run.
"""

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CANDIDATES = [HERE.parent / "results", HERE.parent / ".claude" / "worktrees" / "code-compiler" / "results"]


def results_dir():
    for c in CANDIDATES:
        if (c / "gate_savings_all.csv").exists():
            return c
    return None


def popcount(v):
    return bin(v).count("1")


# ------------------------------------------------------------------ Step 3: gates removed

def predicted_saved(family, n, seed):
    """Gates removed by the pass alone (treatment 'pass' vs 'raw'), all-to-all coupling.

    Returns (gates_saved, two_qubit_saved, absorbed, flushed) or None when no closed form
    is claimed for that family (the random families).
    """
    if family == "ghz" or family == "repetition_syndrome":
        return 0, 0, 0, 0
    if family == "bernstein_vazirani":
        return 1, 0, 1, 0                           # the ancilla's X
    if family == "qft":
        s = n // 2                                  # final SWAP network
        return s, s, s, 0
    if family == "shor15":
        s = n // 2                                  # inverse-QFT SWAPs; work-register X is flushed
        return s, s, s + 1, 1
    marked = seed % 2 ** n
    z = n - popcount(marked)                        # qubits the oracle wraps in X
    if family == "grover_clifford_t":
        return 2 * z + 2 * n, 0, 2 * z + 2 * n, 0   # every X absorbed (Theorem 2)
    if family == "grover_native_mcz":
        # oracle X-before: z flushed (no saving). Oracle X-after + diffusion X-before merge on
        # z qubits into one Y flush: n flushed instead of z+n, saving z. Diffusion X-after: n
        # absorbed into the measurement correction.
        return n + z, 0, 2 * z + 2 * n, z + n
    return None


def predicted_raw(family, n, seed):
    """Raw gate counts (without measurements) where the construction is simple enough to count."""
    if family == "qft":
        return n + n * (n - 1) // 2 + n // 2
    if family == "ghz":
        return n
    if family == "grover_native_mcz":
        z = n - popcount(seed % 2 ** n)
        return 3 * n + 2 * z + 2 * n + 2            # H, X^z, MCZ, X^z, H, X^n, MCZ, X^n, H
    if family == "grover_clifford_t":
        z = n - popcount(seed % 2 ** n)
        tof = 2 * (n - 1) - 3 if n >= 4 else (1 if n == 3 else None)   # v-chain Toffolis
        if tof is None:
            return None
        mcz = 2 + 15 * tof                          # H, chain, H; each Toffoli = 2H + 6CX + 7T
        return 3 * n + 2 * z + 2 * n + 2 * mcz
    return None


def check_gate_savings(rdir):
    rows = list(csv.DictReader(open(rdir / "gate_savings_all.csv")))
    by = defaultdict(dict)
    for r in rows:
        by[(r["family"], int(r["n"]), int(r["seed"]))][r["treatment"]] = r
    n_checked, families = 0, set()
    for (fam, n, seed), t in sorted(by.items()):
        pred = predicted_saved(fam, n, seed)
        if pred is None or "raw" not in t or "pass" not in t:
            continue
        raw, p = t["raw"], t["pass"]
        got = (int(raw["gates"]) - int(p["gates"]), int(raw["two_qubit"]) - int(p["two_qubit"]),
               int(p["absorbed"]), int(p["flushed"]))
        assert got == pred, (fam, n, seed, "measured", got, "predicted", pred)
        pr = predicted_raw(fam, n, seed)
        if pr is not None:
            assert int(raw["gates"]) == pr, (fam, n, seed, raw["gates"], pr)
        if fam == "grover_clifford_t":
            assert int(raw["qubits"]) == 2 * n - 3
        n_checked += 1
        families.add(fam)
    print(f"[ok] S3   gates saved, 2q saved, absorbed and flushed match the closed forms on "
          f"{n_checked} (family, n, seed) rows: {', '.join(sorted(families))}")

    # Regression slopes implied by the closed forms (what the plan's data analysis reports).
    reg = {(r["family"], r["comparison"], r["metric"]): r
           for r in csv.DictReader(open(rdir / "gate_savings_all_regression.csv"))}
    expect = {("qft", "gates"): 0.5, ("qft", "two_qubit"): 0.5, ("qft", "cx_equiv"): 1.5,
              ("bernstein_vazirani", "gates"): 0.0, ("ghz", "gates"): 0.0,
              ("repetition_syndrome", "gates"): 0.0, ("shor15", "gates"): 0.5,
              ("grover_clifford_t", "gates"): 4.0}
    for (fam, metric), slope in expect.items():
        r = reg[(fam, "pass_vs_raw", metric)]
        assert abs(float(r["slope"]) - slope) < 0.05, (fam, metric, r["slope"], slope)
    # Grover Clifford+T: seeds 0,1,2 give z = n, n-1, n-1, so mean saved = 4n - 4/3.
    ns = sorted({int(r["n"]) for r in rows if r["family"] == "grover_clifford_t"})
    mean_pred = sum(4 * n - 4 / 3 for n in ns) / len(ns)
    assert abs(float(reg[("grover_clifford_t", "pass_vs_raw", "gates")]["mean_saved"]) - mean_pred) < 0.02
    print("[ok] S3   reported regression slopes equal the predicted slopes "
          "(QFT 1/2 gate, 3/2 CNOT-equiv; Grover Clifford+T 4; BV/GHZ/syndrome 0; Shor-15 1/2)")


def check_depth(rdir):
    rows = list(csv.DictReader(open(rdir / "gate_savings_all.csv")))
    by = defaultdict(dict)
    for r in rows:
        by[(r["family"], int(r["n"]), int(r["seed"]))][r["treatment"]] = r
    n_checked = 0
    for (fam, n, seed), t in by.items():
        if "raw" not in t or "pass" not in t:
            continue
        saved = int(t["raw"]["depth"]) - int(t["pass"]["depth"])
        if fam in ("ghz", "repetition_syndrome"):
            pred = 0
        elif fam in ("bernstein_vazirani", "qft"):
            pred = 1
        elif fam == "grover_native_mcz":
            pred = 2
        elif fam == "grover_clifford_t":
            pred = 4 if not (seed % 2 ** n >> 1) & 1 else 2      # qubit 1 wrapped in X?
        else:
            continue
        assert saved == pred, (fam, n, seed, saved, pred)
        n_checked += 1
    print(f"[ok] S3   depth saved matches the layer argument on {n_checked} rows "
          f"(BV/QFT 1, GHZ/syndrome 0, Grover native 2, Grover Clifford+T 4 or 2 by qubit 1's X)")


def check_line(rdir):
    path = rdir / "gate_savings_line.csv"
    if not path.exists():
        print("[--] S3   no line-coupling results; connectivity check skipped")
        return
    rows = list(csv.DictReader(open(path)))
    by = defaultdict(dict)
    for r in rows:
        by[(r["family"], int(r["n"]), int(r["seed"]))][r["treatment"]] = r
    worse = defaultdict(int)
    total = defaultdict(int)
    gate_exceptions = []
    for (fam, n, seed), t in by.items():
        if "o3" not in t or "pass+o3" not in t:
            continue
        d = int(t["pass+o3"]["cx_equiv"]) - int(t["o3"]["cx_equiv"])
        total[fam] += 1
        if d > 0:
            worse[fam] += 1
        # total gate count: lower or equal in every row but one (random Clifford n=4 seed 1,
        # 31 -> 33 gates at equal CNOT count, from O3's one-qubit resynthesis)
        if int(t["pass+o3"]["gates"]) > int(t["o3"]["gates"]):
            gate_exceptions.append((fam, n, seed, int(t["o3"]["gates"]), int(t["pass+o3"]["gates"])))
        # and never changes the pass-vs-raw numbers (routing happens after the pass)
        assert predicted_saved(fam, n, seed) is None or \
            int(t["raw"]["gates"]) - int(t["pass"]["gates"]) == predicted_saved(fam, n, seed)[0]
    structured = [f for f in total if f != "random_layered"]
    assert all(worse[f] == 0 for f in structured), {f: worse[f] for f in structured if worse[f]}
    assert worse["random_layered"] > 0            # the heuristic-router effect is real
    assert gate_exceptions == [("random_clifford_synth", 4, 1, 31, 33)], gate_exceptions
    print(f"[ok] S3   line coupling: routed pass+O3 is never worse than O3 in CNOT-equivalents on "
          f"{sum(total[f] for f in structured)} structured rows; random layered is worse in "
          f"{worse['random_layered']} of {total['random_layered']} rows (heuristic routing, section 6); "
          f"total gates lower or equal in {sum(total.values()) - len(gate_exceptions)} of {sum(total.values())} rows")


# ------------------------------------------------------------------ Step 6: noise lemma

def check_noise_lemma():
    import numpy as np
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    Y = 1j * X @ Z

    def channel(kraus):
        def apply(rho):
            return sum(k @ rho @ k.conj().T for k in kraus)
        return apply

    def conj(kraus, P):
        return [P @ k @ P for k in kraus]

    def same(c1, c2):
        rng = np.random.default_rng(0)
        for _ in range(20):
            v = rng.normal(size=2) + 1j * rng.normal(size=2)
            rho = np.outer(v, v.conj()) / np.vdot(v, v)
            if not np.allclose(c1(rho), c2(rho)):
                return False
        return True

    p, g = 0.3, 0.4
    deph = [math.sqrt(1 - p) * I2, math.sqrt(p) * Z]
    depol = [math.sqrt(1 - p) * I2] + [math.sqrt(p / 3) * P for P in (X, Y, Z)]
    ad = [np.diag([1, math.sqrt(1 - g)]).astype(complex), math.sqrt(g) * np.array([[0, 1], [0, 0]], dtype=complex)]
    for P in (X, Y, Z):
        assert same(channel(deph), channel(conj(deph, P)))
        assert same(channel(depol), channel(conj(depol, P)))
    assert same(channel(ad), channel(conj(ad, Z)))
    assert not same(channel(ad), channel(conj(ad, X)))
    assert not same(channel(ad), channel(conj(ad, Y)))
    one = np.diag([0.0, 1.0]).astype(complex)
    zero = np.diag([1.0, 0.0]).astype(complex)
    assert np.allclose(channel(ad)(zero), zero)                        # |0> is the fixed point
    assert np.allclose(channel(conj(ad, X))(one), one)                 # with an X flag, |1> is
    assert np.allclose(channel(ad)(one), np.diag([g, 1 - g]))          # stored directly, |1> decays
    print("[ok] S6   Lemma 4: dephasing and depolarizing commute with every Pauli flag; amplitude "
          "damping commutes with Z but not X/Y, and an X flag makes logical |1> its fixed point")


# ------------------------------------------------------------------ Step 4: lookup cost

def bucket_brigade_counts(N, d):
    """Derived in predictions.md section 2 from the qubit-router model."""
    k = int(math.log2(N))
    assert 2 ** k == N
    addr_routing = 2 * (N - 1 - k)                  # 2 CSWAP per internal node below the root... see text
    bus_pass = 2 * (N - 1)                          # one descent or ascent of the bus
    cswap = 2 * addr_routing + 4 * d * bus_pass     # route down + up; 2 descents + 2 ascents per data bit
    qubits = k + d + 2 * (N - 1) + N + N * d + d
    return {"cswap": cswap, "toffoli": cswap, "t_count": 7 * cswap, "cnot_data": 2 * N * d,
            "qubits": qubits, "x_absorbable": cswap}


def count_table_lookup_by_recursion(N):
    """Count Toffolis and X gates by walking the unary-iteration recursion exactly as written."""
    k = int(math.log2(N))
    tof = x = 0

    def iterate(level):
        nonlocal tof, x
        if level < 0:
            return
        x += 2
        tof += 2
        iterate(level - 1)
        iterate(level - 1)

    if k == 1:
        x = 2
    else:
        x += 2
        iterate(k - 2)
        iterate(k - 2)
    return tof, x


def table_lookup_counts(N, d):
    k = int(math.log2(N))
    tof = 2 * (N - 2) if k >= 2 else 0
    return {"toffoli": tof, "toffoli_babbush": N - 1, "t_count": 7 * tof,
            "t_count_babbush": 4 * (N - 1), "qubits": k + d + max(k - 1, 0),
            "x_absorbable": 2 * (N - 1), "cnot_data_max": N * d}


def select_swap_counts(N, d, lam):
    k = int(math.log2(N))
    return {"toffoli": math.ceil(N / lam) + d * (lam - 1), "qubits": k + d * lam + max(k - 1, 0)}


def check_oracle_cost(rdir):
    rows = list(csv.DictReader(open(rdir / "oracle_cost.csv")))
    n_rows = 0
    for r in rows:
        N, d = int(r["N"]), int(r["d"])
        if r["source"] != "closed_form":
            continue
        if r["design"] == "bucket_brigade":
            c = bucket_brigade_counts(N, d)
            for key in ("cswap", "toffoli", "t_count", "cnot_data", "qubits", "x_absorbable"):
                assert int(r[key]) == c[key], (r["design"], N, d, key, r[key], c[key])
        elif r["design"] == "table_lookup":
            c = table_lookup_counts(N, d)
            for key in ("toffoli", "toffoli_babbush", "t_count", "t_count_babbush", "qubits",
                        "x_absorbable", "cnot_data_max"):
                assert int(r[key]) == c[key], (r["design"], N, d, key, r[key], c[key])
            tof, x = count_table_lookup_by_recursion(N)
            assert tof == c["toffoli"] and x == c["x_absorbable"]
        elif r["design"] == "select_swap":
            c = select_swap_counts(N, d, int(r["lambda"]))
            assert int(r["toffoli"]) == c["toffoli"] and int(r["qubits"]) == c["qubits"]
        n_rows += 1
    # Built circuits agree with the closed forms on the quantities both report.
    for r in rows:
        if r["source"] == "built":
            N, d = int(r["N"]), int(r["d"])
            c = bucket_brigade_counts(N, d) if r["design"] == "bucket_brigade" else table_lookup_counts(N, d)
            assert int(r["toffoli"]) == c["toffoli"] and int(r["qubits"]) == c["qubits"]
    print(f"[ok] S4   {n_rows} closed-form rows of oracle_cost.csv match the independently derived "
          f"counts; table-lookup Toffoli/X counts also match a direct walk of the recursion")

    # What the pass does to the lookups themselves (built rows).
    fr = []
    for r in rows:
        if r["source"] != "built":
            continue
        N, x, fl = int(r["N"]), int(r["x"]), int(r["x_flushed"])
        assert int(r["x_absorbed"]) == x                      # every X is absorbed at first...
        if r["design"] == "table_lookup":
            # ...and in the unary iteration every one is flushed again: each X sits next to a
            # Toffoli that uses the same qubit as a control (section 7: not absorbable).
            assert fl == x == 2 * (N - 1)
            assert int(r["gates_after_pass"]) == int(r["gates"])
        else:
            # bucket brigade (predictions.md 2a): per router with m address passes and d data
            # bits the pass flushes 4m + 4d of 4m + 8d X gates; the 4d that cancel are the pairs
            # where a descent's trailing X meets an ascent's leading X. Summed over the tree,
            # with w set bits in the table (each costs 2d flushes on its leaf):
            k = int(math.log2(N))
            struct_x = 4 * (N - 1 - k) + 8 * d * (N - 1)
            w = x - struct_x                                  # data-preparation X gates
            assert 0 <= w <= N * d
            pred_flushed = 4 * (N - 1 - k) + 4 * d * (N - 1) + 2 * d * w
            swaps = 2 * (N - 1 + k) + 4 * d
            assert fl == pred_flushed, (N, fl, pred_flushed)
            assert int(r["gates"]) - int(r["gates_after_pass"]) == x + swaps - fl, (N,)
            fr.append(fl / x)
    print("[ok] S4   table lookup: the pass removes nothing (all 2(N-1) X gates flushed), as predicted; "
          f"bucket brigade: flushed = 4(N-1-k) + 4d(N-1) + 2dw and removed = X + SWAP - flushed on "
          f"every built row ({min(fr):.0%}-{max(fr):.0%} of X flushed, tending to 2/3)")

    # Select-swap optimum: N/lam + d*lam is minimised at lam = sqrt(N/d), value 2 sqrt(N d).
    for N in (8, 64, 1024):
        for d in (1, 8):
            lam_star = math.sqrt(N / d)
            best = min(N / lam + d * lam for lam in [x / 100 for x in range(1, 100 * N)])
            assert abs(best - 2 * math.sqrt(N * d)) < 0.05
            assert abs((N / lam_star + d * lam_star) - 2 * math.sqrt(N * d)) < 1e-9
    print("[ok] S4   select-swap: N/lam + d*lam has its minimum 2*sqrt(N d) at lam = sqrt(N/d)")


# ------------------------------------------------------------------ Step 5: shared oracle

def check_shared_oracle(rdir):
    rows = list(csv.DictReader(open(rdir / "shared_oracle.csv")))
    n_rows = 0
    for r in rows:
        N, k = int(r["N"]), int(r["k"])
        a = int(math.log2(N))
        Q, L = int(r["oracle_qubits"]), int(r["lookup_depth"])
        I = int(r["iterations"])
        D = 2 * a + 1
        assert I == max(1, math.floor(math.pi / 4 * math.sqrt(N)))
        dup_q, dup_d = k * (a + Q), I * (L + D)
        sh_q, sh_d = k * a + Q, I * (k * L + D)
        assert (int(r["dup_qubits"]), int(r["dup_depth"])) == (dup_q, dup_d), (r, dup_q, dup_d)
        assert (int(r["shared_qubits"]), int(r["shared_depth"])) == (sh_q, sh_d)
        # Closed forms for the two ratios.
        frac = (1 - 1 / k) * Q / (a + Q)
        assert abs(float(r["qubits_saved_frac"]) - frac) < 1e-3
        assert abs(float(r["depth_ratio"]) - (k * L + D) / (L + D)) < 1e-3
        # Break-even rule: shared has smaller qubit*depth  iff  k * a * L < Q * D  (k >= 2).
        if k >= 2:
            shared_wins = sh_q * sh_d < dup_q * dup_d
            assert shared_wins == (k * a * L < Q * D), (N, k, shared_wins)
        n_rows += 1
    print(f"[ok] S5   {n_rows} rows of shared_oracle.csv match the closed forms; the break-even rule "
          f"k < Q*D / (log2(N) * L) predicts the qubit*depth winner on every row")

    # The whitepaper's "66% fewer qubits for 3 instances": the bound is 1 - 1/k = 2/3, approached
    # as the oracle dwarfs the address register.
    for r in rows:
        if int(r["k"]) == 3 and int(r["N"]) == 1024:
            f = float(r["qubits_saved_frac"])
            if r["design"] == "bucket_brigade":
                assert 0.66 < f < 2 / 3                 # Q = 4096 >> log2 N = 10
            else:
                assert abs(f - 1 / 3) < 1e-3            # Q = log2 N = 10: (2/3) * 10/20
    print("[ok] S5   k=3, N=1024: bucket brigade saves 66.5% (just under the 2/3 ceiling), "
          "table lookup saves exactly 1/3")


# ------------------------------------------------------------------ self-consistency

def check_self_consistency():
    # Bucket-brigade address routing: 2 * sum_{j<k} (2^j - 1) = 2 (N - 1 - k).
    for k in range(1, 12):
        N = 2 ** k
        assert 2 * sum(2 ** j - 1 for j in range(k)) == 2 * (N - 1 - k)
    # Table lookup: the recursion visits N - 2 internal nodes below the root.
    for k in range(2, 11):
        tof, x = count_table_lookup_by_recursion(2 ** k)
        assert tof == 2 * (2 ** k - 2) and x == 2 * (2 ** k - 1)
    # Grover v-chain: k controls -> 2k - 3 Toffolis (k >= 3).
    for n in range(4, 70):
        k = n - 1
        steps = 1 + (k - 3)
        assert 2 * steps + 1 == 2 * n - 5
    print("[ok] SELF sums behind the closed forms hold for N = 2..2048 and n = 4..69")


def check_whitepaper_table():
    """Whitepaper Figure 6 table: conventional 3N, oracle-based N + 3 ceil(log2 N), for k = 3."""
    table = {128: 0.61, 256: 0.64, 512: 0.65, 1024: 0.66}
    for N, frac in table.items():
        a = math.ceil(math.log2(N))
        assert abs((1 - (N + 3 * a) / (3 * N)) - frac) < 0.005          # reproduces the table
        corrected = (1 - 1 / 3) * N / (a + N)                            # both designs count addresses
        assert 1 - (N + 3 * a) / (3 * N) < corrected < 2 / 3
    print("[ok] WP   whitepaper's 61/64/65/66% reproduce 1 - (N + 3 log N)/(3N); counting the address "
          "registers on both sides gives (2/3) N/(N + log N) = 63.2/64.9/65.8/66.2%, still under 2/3")


if __name__ == "__main__":
    check_self_consistency()
    check_whitepaper_table()
    rdir = results_dir()
    if rdir is None:
        print("no results directory found; only self-consistency checks were run")
        sys.exit(0)
    print(f"using results in {rdir}")
    check_gate_savings(rdir)
    check_depth(rdir)
    check_line(rdir)
    check_noise_lemma()
    check_oracle_cost(rdir)
    check_shared_oracle(rdir)
    print("all prediction checks passed")
