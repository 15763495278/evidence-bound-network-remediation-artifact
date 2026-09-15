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
import math
from collections import defaultdict
from pathlib import Path

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


def compare(rows, arm_key, hi, lo, pair_key, success_key):
    by_pair = defaultdict(dict)
    for r in rows:
        by_pair[pair_key(r)][r[arm_key]] = bool(success_key(r))
    units = [v for v in by_pair.values() if hi in v and lo in v]
    b = sum(v[hi] and not v[lo] for v in units)
    c = sum(v[lo] and not v[hi] for v in units)
    hi_s = sum(v[hi] for v in units)
    lo_s = sum(v[lo] for v in units)
    n = len(units)
    return {
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

    f1 = load_json_dir(ROOTS["F1"])
    key1 = lambda r: (r.get("model"), r.get("task_id"), r.get("env_seed"), r.get("inference_seed"))
    results["F1_C1_M2_vs_M0"] = compare(f1, "method_id", "M2-static-grounded", "M0-free-agent", key1, f1_success)
    results["F1_C2_M3_vs_M1"] = compare(f1, "method_id", "M3-ours", "M1-rc-active", key1, f1_success)

    f2b = load_json_dir(ROOTS["F2B"])
    key2b = lambda r: ((r.get("plan") or {}).get("model"), (r.get("plan") or {}).get("task_id"), (r.get("plan") or {}).get("env_seed"))
    method_key = "method"
    results["F2B_ReadyBind_vs_RCFree"] = compare(f2b, method_key, "Ready-Bind", "RC-Free", key2b, f2b_success)

    f2c = load_json_dir(ROOTS["F2C"])
    key2c = lambda r: (r.get("plan") or {}).get("pair_id")
    results["F2C_ReadyBind_vs_DGFree"] = compare(f2c, "method", "DG-ReadyBind", "DG-Free", key2c, f2c_success)

    f3 = load_jsonl(ROOTS["F3"])
    key3 = lambda r: (r.get("model"), r.get("task_id"), r.get("env_seed"), r.get("replication_seed"))
    results["F3_Full_vs_ReadinessOnly"] = compare(f3, "arm_id", "A3-Full", "A1-ReadinessOnly", key3, lambda r: r.get("safe_effective_repair"))
    results["F3_Full_vs_BinderOnly"] = compare(f3, "arm_id", "A3-Full", "A2-BinderOnly", key3, lambda r: r.get("safe_effective_repair"))

    f4 = load_jsonl(ROOTS["F4"])
    key4 = lambda r: (r.get("task_id"), r.get("env_seed"), r.get("replication_seed"))
    results["F4_ReadyBind_vs_TypedGuard"] = compare(f4, "arm_id", "B2-ReadyBind", "B1-TypedGuard", key4, lambda r: r.get("safe_effective_repair"))

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
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if not all(checks.values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
