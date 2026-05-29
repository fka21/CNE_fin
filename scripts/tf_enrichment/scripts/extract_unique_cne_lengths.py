#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Extract unique CNE lengths from a final_table.tsv file.")
    p.add_argument("--cne-table", required=True, help="Input final_table.tsv")
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.cne_table, sep="\t", dtype=str)
    required = {"seqnames", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["width"] = pd.to_numeric(df.get("width", pd.Series([None] * len(df))), errors="coerce")
    df = df.dropna(subset=["seqnames", "start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    u = df.drop_duplicates(["seqnames", "start", "end"]).copy()
    calc_len = (u["end"] - u["start"]).abs() + 1
    u["length"] = u["width"].fillna(calc_len)
    u.loc[u["length"] <= 0, "length"] = calc_len
    u["length"] = u["length"].astype(int)

    lengths = u["length"].sort_values().reset_index(drop=True)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    lengths.to_csv(str(out_prefix) + "_lengths.txt", index=False, header=False)
    u.to_csv(str(out_prefix) + "s.tsv", sep="\t", index=False)

    stats = {
        "Unique_CNE_count": len(lengths),
        "Min_length": int(lengths.min()),
        "Q1_length": float(lengths.quantile(0.25)),
        "Median_length": float(lengths.median()),
        "Mean_length": float(lengths.mean()),
        "Q3_length": float(lengths.quantile(0.75)),
        "Max_length": int(lengths.max()),
    }
    with open(str(out_prefix) + "_length_stats.txt", "w") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")

    print("Done.")
    print(f"Output lengths: {out_prefix}_lengths.txt")
    print(f"Output stats: {out_prefix}_length_stats.txt")
    for k, v in stats.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
