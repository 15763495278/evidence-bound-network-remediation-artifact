#!/usr/bin/env python3
"""Independent one-click recomputation of the paper's reported effects.

The script reads only frozen per-run records (or, for F2-B, its frozen
per-run directory) and recomputes SER rates, paired discordances, risk
differences, and exact two-sided McNemar p-values.  It does not import any
experiment runner or analyzer.
"""
from __future__ import annotations

import csv
import json
import itertools
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent

ROOTS = {
    "F1": BASE / "data" / "f1",
    "F2B": BASE / "data" / "f2b",
    "F2C": BASE / "data" / "f2c",
    "F3": BASE / "data" / "f3" / "runs-summary-all-planned.jsonl",
    "F4": BASE / "data" / "f4" / "runs-summary-all-planned.jsonl",
}

OUT = Path(__file__).resolve().parent / "generated"


def load_json_dir(path: Path):
    if not path.is_dir():
        raise FileNotFoundError(path)
    return [json.loads(p.read_text(encoding="utf-8-sig")) for p in sorted(path.glob("*.json"))]


def load_jsonl(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def exact_mcnemar(b: int, c: int) -> float:
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    lower = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return min(1.0, 2.0 * lower)


def compare(rows, arm_key, hi, lo, pair_key, success_key, cluster_key):
    by_pair = defaultdict(dict)
    pair_clusters = {}
    for r in rows:
        key = pair_key(r)
        by_pair[key][r[arm_key]] = bool(success_key(r))
        cluster = cluster_key(r)
        if key in pair_clusters and pair_clusters[key] != cluster:
            raise ValueError(f"pair {key!r} maps to more than one task cluster")
        pair_clusters[key] = cluster
    complete = [(key, v) for key, v in by_pair.items() if hi in v and lo in v]
    units = [v for _, v in complete]
    b = sum(v[hi] and not v[lo] for v in units)
    c = sum(v[lo] and not v[hi] for v in units)
    hi_s = sum(v[hi] for v in units)
    lo_s = sum(v[lo] for v in units)
    n = len(units)
    result = {
        "n_pairs": n,
        "hi": hi,
        "lo": lo,
        "hi_success": hi_s,
        "lo_success": lo_s,
        "b_hi_wins": b,
        "c_lo_wins": c,
        "ties": n - b - c,
        "paired_risk_difference": (hi_s - lo_s) / n if n else None,
        "exact_mcnemar_p": exact_mcnemar(b, c),
    }
    paired_units = [
        {
            "pair": repr(key),
            "task_template_id": pair_clusters[key],
            "treatment_success": int(v[hi]),
            "comparator_success": int(v[lo]),
            "difference": int(v[hi]) - int(v[lo]),
        }
        for key, v in complete
    ]
    return result, paired_units


def task_cluster_sensitivity(paired_units, seed=20260930, draws=100_000):
    """Equal-weight task-template analysis frozen for the manuscript sensitivity check."""
    grouped = defaultdict(list)
    for unit in paired_units:
        grouped[unit["task_template_id"]].append(float(unit["difference"]))
    if not grouped:
        raise ValueError("no task clusters")
    cluster_ids = sorted(grouped)
    means = np.array([np.mean(grouped[c]) for c in cluster_ids], dtype=float)
    observed = float(means.mean())

    rng = np.random.Generator(np.random.PCG64(seed))
    sampled = rng.integers(0, len(means), size=(draws, len(means)))
    bootstrap = means[sampled].mean(axis=1)
    ci_low, ci_high = np.percentile(bootstrap, [2.5, 97.5])

    signflip = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(means)):
        signflip.append(float(np.mean(means * np.asarray(signs))))
    p_two_sided = sum(abs(v) >= abs(observed) - 1e-15 for v in signflip) / len(signflip)

    if len(means) > 1:
        loto = [float(np.delete(means, i).mean()) for i in range(len(means))]
    else:
        loto = [observed]
    return {
        "n_task_clusters": len(means),
        "cluster_rd": observed,
        "cluster_bootstrap_ci_low": float(ci_low),
        "cluster_bootstrap_ci_high": float(ci_high),
        "cluster_signflip_p_two_sided": float(p_two_sided),
        "n_positive_clusters": int(np.sum(means > 0)),
        "n_zero_clusters": int(np.sum(means == 0)),
        "n_negative_clusters": int(np.sum(means < 0)),
        "leave_one_cluster_rd_min": min(loto),
        "leave_one_cluster_rd_max": max(loto),
        "bootstrap_seed": seed,
        "bootstrap_draws": draws,
        "cluster_means": [
            {
                "task_template_id": cid,
                "n_pairs": len(grouped[cid]),
                "cluster_rd": float(mean),
            }
            for cid, mean in zip(cluster_ids, means)
        ],
    }


def f1_success(r):
    repairs = r.get("repairs") or []
    if isinstance(repairs, dict):
        repairs = [repairs]
    oracle_ok = any(bool(x.get("all_oracles_pass")) for x in repairs)
    return oracle_ok and bool((r.get("cleanup") or {}).get("clean")) and not r.get("agent_failure") and not r.get("infrastructure_error")


def f2b_success(r):
    return r.get("verdict") == "SER_PASS" or r.get("safe_effective_repair") is True


def f2c_success(r):
    return r.get("verdict") == "SER_PASS"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = {}
    paired = {}

    f1 = load_json_dir(ROOTS["F1"])
    key1 = lambda r: (r.get("model"), r.get("task_id"), r.get("env_seed"), r.get("inference_seed"))
    results["F1_C1_M2_vs_M0"], paired["F1_C1_M2_vs_M0"] = compare(f1, "method_id", "M2-static-grounded", "M0-free-agent", key1, f1_success, lambda r: r.get("task_id"))
    results["F1_C2_M3_vs_M1"], paired["F1_C2_M3_vs_M1"] = compare(f1, "method_id", "M3-ours", "M1-rc-active", key1, f1_success, lambda r: r.get("task_id"))

    f2b = load_json_dir(ROOTS["F2B"])
    key2b = lambda r: ((r.get("plan") or {}).get("model"), (r.get("plan") or {}).get("task_id"), (r.get("plan") or {}).get("env_seed"))
    method_key = "method"
    results["F2B_ReadyBind_vs_RCFree"], paired["F2B_ReadyBind_vs_RCFree"] = compare(f2b, method_key, "Ready-Bind", "RC-Free", key2b, f2b_success, lambda r: (r.get("plan") or {}).get("task_id"))

    f2c = load_json_dir(ROOTS["F2C"])
    key2c = lambda r: (r.get("plan") or {}).get("pair_id")
    results["F2C_ReadyBind_vs_DGFree"], paired["F2C_ReadyBind_vs_DGFree"] = compare(f2c, "method", "DG-ReadyBind", "DG-Free", key2c, f2c_success, lambda r: (r.get("plan") or {}).get("task_id"))

    f3 = load_jsonl(ROOTS["F3"])
    key3 = lambda r: (r.get("model"), r.get("task_id"), r.get("env_seed"), r.get("replication_seed"))
    results["F3_Full_vs_ReadinessOnly"], paired["F3_Full_vs_ReadinessOnly"] = compare(f3, "arm_id", "A3-Full", "A1-ReadinessOnly", key3, lambda r: r.get("safe_effective_repair"), lambda r: r.get("task_id"))
    results["F3_Full_vs_BinderOnly"], paired["F3_Full_vs_BinderOnly"] = compare(f3, "arm_id", "A3-Full", "A2-BinderOnly", key3, lambda r: r.get("safe_effective_repair"), lambda r: r.get("task_id"))

    f4 = load_jsonl(ROOTS["F4"])
    key4 = lambda r: (r.get("task_id"), r.get("env_seed"), r.get("replication_seed"))
    results["F4_ReadyBind_vs_TypedGuard"], paired["F4_ReadyBind_vs_TypedGuard"] = compare(f4, "arm_id", "B2-ReadyBind", "B1-TypedGuard", key4, lambda r: r.get("safe_effective_repair"), lambda r: r.get("task_id"))

    expected = {
        "F1_C1_M2_vs_M0": (72, 39, 1, 0.527778),
        "F1_C2_M3_vs_M1": (72, 53, 0, 0.736111),
        "F2B_ReadyBind_vs_RCFree": (54, 3, 4, -0.018519),
        "F2C_ReadyBind_vs_DGFree": (54, 16, 3, 0.240741),
        "F3_Full_vs_ReadinessOnly": (36, 18, 0, 0.5),
        "F3_Full_vs_BinderOnly": (36, 3, 0, 0.083333),
        "F4_ReadyBind_vs_TypedGuard": (36, 0, 1, -0.027778),
    }
    checks = {}
    for name, vals in expected.items():
        got = results[name]
        checks[name] = (
            got["n_pairs"] == vals[0]
            and got["b_hi_wins"] == vals[1]
            and got["c_lo_wins"] == vals[2]
            and abs(got["paired_risk_difference"] - vals[3]) < 1e-6
        )
    payload = {"source": "frozen per-run records", "results": results, "checks": checks, "all_checks_pass": all(checks.values())}
    (OUT / "recomputed_results.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    with (OUT / "effect_estimates.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["comparison", "n_pairs", "hi_success", "lo_success", "b_hi_wins", "c_lo_wins", "ties", "paired_risk_difference", "exact_mcnemar_p", "check_pass"])
        w.writeheader()
        for name, row in results.items():
            w.writerow({"comparison": name, **{k: row[k] for k in w.fieldnames if k in row}, "check_pass": checks[name]})
    lines = ["# Independent recomputation", "", f"Overall: **{'PASS' if all(checks.values()) else 'FAIL'}**", ""]
    for name, row in results.items():
        lines.append(f"- {name}: n={row['n_pairs']}, b/c={row['b_hi_wins']}/{row['c_lo_wins']}, RD={row['paired_risk_difference']:+.4f}, exact p={row['exact_mcnemar_p']:.6g}, check={'PASS' if checks[name] else 'FAIL'}")
    (OUT / "verification_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    cluster_results = {name: task_cluster_sensitivity(units) for name, units in paired.items()}
    cluster_expected = {
        "F1_C1_M2_vs_M0": (12, 0.5277777777777778, 0.3194444444444445, 0.7222222222222222, 0.001953125),
        "F1_C2_M3_vs_M1": (12, 0.7361111111111112, 0.6111111111111110, 0.8611111111111112, 0.00048828125),
        "F2B_ReadyBind_vs_RCFree": (9, -0.01851851851851851, -0.14814814814814814, 0.12962962962962965, 1.0),
        "F2C_ReadyBind_vs_DGFree": (9, 0.24074074074074078, -0.05555555555555555, 0.5555555555555556, 0.25),
        "F3_Full_vs_ReadinessOnly": (6, 0.5, 0.25, 0.75, 0.0625),
        "F3_Full_vs_BinderOnly": (6, 0.08333333333333333, 0.0, 0.25, 1.0),
        "F4_ReadyBind_vs_TypedGuard": (12, -0.027777777777777776, -0.08333333333333333, 0.0, 1.0),
    }
    cluster_checks = {}
    for name, exp in cluster_expected.items():
        got = cluster_results[name]
        cluster_checks[name] = (
            got["n_task_clusters"] == exp[0]
            and abs(got["cluster_rd"] - exp[1]) < 1e-12
            and abs(got["cluster_bootstrap_ci_low"] - exp[2]) < 1e-12
            and abs(got["cluster_bootstrap_ci_high"] - exp[3]) < 1e-12
            and abs(got["cluster_signflip_p_two_sided"] - exp[4]) < 1e-12
        )
    cluster_payload = {
        "method": "equal-weight task-template cluster means; 100,000 cluster bootstrap draws; exact sign-flip enumeration",
        "rng": "NumPy PCG64",
        "seed": 20260930,
        "results": cluster_results,
        "checks": cluster_checks,
        "all_checks_pass": all(cluster_checks.values()),
    }
    (OUT / "cluster_sensitivity_results.json").write_text(json.dumps(cluster_payload, indent=2, ensure_ascii=False), encoding="utf-8")
    with (OUT / "cluster_sensitivity_summary.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["comparison", "n_task_clusters", "cluster_rd", "cluster_bootstrap_ci_low", "cluster_bootstrap_ci_high", "cluster_signflip_p_two_sided", "n_positive_clusters", "n_zero_clusters", "n_negative_clusters", "leave_one_cluster_rd_min", "leave_one_cluster_rd_max", "bootstrap_seed", "bootstrap_draws", "check_pass"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, row in cluster_results.items():
            w.writerow({"comparison": name, **{k: row[k] for k in fields if k in row}, "check_pass": cluster_checks[name]})

    overall_pass = all(checks.values()) and all(cluster_checks.values())
    lines.extend(["", "## Task-template sensitivity", ""])
    for name, row in cluster_results.items():
        lines.append(f"- {name}: G={row['n_task_clusters']}, RD={row['cluster_rd']:+.4f}, 95% cluster bootstrap CI=[{row['cluster_bootstrap_ci_low']:+.4f}, {row['cluster_bootstrap_ci_high']:+.4f}], exact sign-flip p={row['cluster_signflip_p_two_sided']:.6g}, check={'PASS' if cluster_checks[name] else 'FAIL'}")
    lines[2] = f"Overall: **{'PASS' if overall_pass else 'FAIL'}**"
    (OUT / "verification_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    print(json.dumps(cluster_payload, indent=2, ensure_ascii=False))
    if not overall_pass:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
