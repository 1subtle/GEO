#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
import re

from summarize_results_to_xlsx import summarize as summarize_xlsx
from summarize_results_to_xlsx import write_xlsx


def instance_key(name):
    m = re.search(r"B(\d+)_N(\d+)_R(\d+)", name)
    if not m:
        return (10**9, 10**9, 10**9, name)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)), name)


def parse_instance(name):
    m = re.search(r"B(\d+)_N(\d+)_R(\d+)", name)
    if not m:
        return "", "", ""
    return m.group(1), m.group(2), m.group(3)


def fmt_float(value, digits=2):
    return f"{value:.{digits}f}"


def fmt_time(value):
    return f"{value:.4f}".rstrip("0").rstrip(".")


def sample_std(values):
    n = len(values)
    if n <= 1:
        return 0.0
    avg = sum(values) / n
    return math.sqrt(sum((x - avg) ** 2 for x in values) / (n - 1))


def find_latest_descent_results_csv():
    candidates = sorted(
        glob.glob("descent_batch_results_*/results.csv"),
        key=os.path.getmtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError("No descent_batch_results_*/results.csv found in current directory")
    return candidates[0]


def summarize(src):
    groups = {}
    with open(src, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"instance", "seed", "best_profit", "best_time", "verified"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{src} missing required columns: {sorted(missing)}")
        for row in reader:
            groups.setdefault(row["instance"], []).append(row)

    rows = []
    for inst, items in groups.items():
        vals = []
        for item in items:
            vals.append(
                {
                    "seed": int(item["seed"]),
                    "best_profit": int(item["best_profit"]),
                    "best_time": float(item["best_time"]),
                    "verified": item["verified"],
                }
            )

        profits = [v["best_profit"] for v in vals]
        best = max(profits)
        avg = sum(profits) / len(profits)
        std = sample_std(profits)

        best_runs = [v for v in vals if v["best_profit"] == best]
        best_time = min(v["best_time"] for v in best_runs)
        best_seed = min(v["seed"] for v in best_runs if v["best_time"] == best_time)
        verified_runs = sum(1 for v in vals if v["verified"] == "YES")
        b, n, r = parse_instance(inst)

        rows.append(
            {
                "instance": inst,
                "B": b,
                "N": n,
                "R": r,
                "time_to_best_s": fmt_time(best_time),
                "best": str(best),
                "avg": fmt_float(avg),
                "std": fmt_float(std),
                "best_seed": str(best_seed),
                "runs": str(len(vals)),
                "verified_runs": str(verified_runs),
            }
        )

    rows.sort(key=lambda row: instance_key(row["instance"]))
    return rows


def write_outputs(rows, out_prefix):
    md_path = out_prefix + ".md"

    with open(md_path, "w", encoding="utf-8") as f:
        f.write("| instance | B | N | R | time_to_best_s | best | avg | std | best_seed | runs | verified_runs |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            f.write(
                f"| {row['instance']} | {row['B']} | {row['N']} | {row['R']} | "
                f"{row['time_to_best_s']} | {row['best']} | {row['avg']} | {row['std']} | "
                f"{row['best_seed']} | {row['runs']} | {row['verified_runs']} |\n"
            )

    return md_path


def main():
    parser = argparse.ArgumentParser(
        description="Summarize GEO descent batch results by instance with time-to-best, best, avg, and sample std."
    )
    parser.add_argument(
        "results_csv",
        nargs="?",
        default=None,
        help="Path to descent batch results.csv. Defaults to latest descent_batch_results_*/results.csv",
    )
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Output path prefix without extension. Defaults to <results_dir>/descent_summary_best_avg_std",
    )
    args = parser.parse_args()

    results_csv = args.results_csv if args.results_csv is not None else find_latest_descent_results_csv()
    src = os.path.abspath(results_csv)
    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = os.path.join(os.path.dirname(src), "descent_summary_best_avg_std")

    rows = summarize(src)
    md_path = write_outputs(rows, out_prefix)
    xlsx_path = out_prefix + ".xlsx"
    write_xlsx(summarize_xlsx(src), xlsx_path)
    print(md_path)
    print(xlsx_path)
    print(f"instances={len(rows)}")
    if rows:
        run_counts = sorted({int(row["runs"]) for row in rows})
        print("runs_per_instance=" + ",".join(str(x) for x in run_counts))


if __name__ == "__main__":
    main()
