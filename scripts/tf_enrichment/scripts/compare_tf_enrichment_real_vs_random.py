#!/usr/bin/env python3

import argparse
import math
import numpy as np
import pandas as pd

try:
    from scipy.stats import fisher_exact
except ImportError:
    raise ImportError("Missing scipy. Install with: pip install scipy or conda install scipy")


def parse_args():
    p = argparse.ArgumentParser(description="Compare TF motif presence between real CNEs and random genomic background.")
    p.add_argument("--real-hits", required=True, help="Real CNE tf_hits.tsv")
    p.add_argument("--real-summary", required=True, help="Real CNE summary.tsv")
    p.add_argument("--bg-hits", required=True, help="Background/random tf_hits.tsv")
    p.add_argument("--bg-summary", required=True, help="Background/random summary.tsv")
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    p.add_argument("--min-score", type=float, default=0.85, help="Minimum relative_score")
    return p.parse_args()


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = np.empty(n, dtype=float)
    cumulative_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        cumulative_min = min(cumulative_min, val)
        adjusted[order[i]] = cumulative_min
    return np.minimum(adjusted, 1.0)


def load_hits(path, min_score):
    df = pd.read_csv(path, sep="\t")
    required = {"cne_id", "tf", "relative_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns from {path}: {missing}")
    df["relative_score"] = pd.to_numeric(df["relative_score"], errors="coerce")
    return df[df["relative_score"] >= min_score].copy()


def load_summary(path):
    df = pd.read_csv(path, sep="\t")
    if "cne_id" not in df.columns:
        raise ValueError(f"Missing cne_id column from {path}")
    return df


def tf_presence_counts(hits):
    return hits[["cne_id", "tf"]].drop_duplicates().groupby("tf")["cne_id"].nunique().to_dict()


def selected_hit_counts(hits):
    return hits.groupby("tf").size().to_dict()


def main():
    args = parse_args()
    real_hits = load_hits(args.real_hits, args.min_score)
    bg_hits = load_hits(args.bg_hits, args.min_score)
    real_summary = load_summary(args.real_summary)
    bg_summary = load_summary(args.bg_summary)
    n_real = real_summary["cne_id"].nunique()
    n_bg = bg_summary["cne_id"].nunique()

    real_presence = tf_presence_counts(real_hits)
    bg_presence = tf_presence_counts(bg_hits)
    real_selected = selected_hit_counts(real_hits)
    bg_selected = selected_hit_counts(bg_hits)
    all_tfs = sorted(set(real_presence) | set(bg_presence))

    rows = []
    for tf in all_tfs:
        a = int(real_presence.get(tf, 0))
        c = int(bg_presence.get(tf, 0))
        b = n_real - a
        d = n_bg - c
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        real_fraction = a / n_real if n_real > 0 else np.nan
        bg_fraction = c / n_bg if n_bg > 0 else np.nan
        real_fraction_pc = (a + 0.5) / (n_real + 1.0)
        bg_fraction_pc = (c + 0.5) / (n_bg + 1.0)
        fold_enrichment = real_fraction_pc / bg_fraction_pc
        log2_enrichment = math.log2(fold_enrichment)
        rows.append({
            "tf": tf,
            "real_cnes_with_tf": a,
            "real_cnes_without_tf": b,
            "bg_cnes_with_tf": c,
            "bg_cnes_without_tf": d,
            "real_fraction": real_fraction,
            "bg_fraction": bg_fraction,
            "fold_enrichment_pc": fold_enrichment,
            "log2_enrichment_pc": log2_enrichment,
            "odds_ratio": odds_ratio,
            "p_value_fisher_greater": p_value,
            "real_selected_hits": int(real_selected.get(tf, 0)),
            "bg_selected_hits": int(bg_selected.get(tf, 0)),
        })

    out = pd.DataFrame(rows)
    out["FDR_BH"] = bh_fdr(out["p_value_fisher_greater"].values)
    out = out.sort_values(["FDR_BH", "p_value_fisher_greater", "log2_enrichment_pc"], ascending=[True, True, False])

    out_tsv = f"{args.out_prefix}.tf_enrichment_real_vs_random_genomic.tsv"
    sig_tsv = f"{args.out_prefix}.tf_enrichment_real_vs_random_genomic.FDR005.enriched.tsv"
    import os
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    out.to_csv(out_tsv, sep="\t", index=False)
    sig = out[(out["FDR_BH"] < 0.05) & (out["log2_enrichment_pc"] > 0)].copy()
    sig.to_csv(sig_tsv, sep="\t", index=False)

    print("Done.")
    print(f"Real CNEs: {n_real}")
    print(f"Background/random CNEs: {n_bg}")
    print(f"TFs tested: {len(out)}")
    print(f"Significant enriched TFs, FDR < 0.05 and log2 enrichment > 0: {len(sig)}")
    print(f"Output: {out_tsv}")
    print(f"Output: {sig_tsv}")


if __name__ == "__main__":
    main()
