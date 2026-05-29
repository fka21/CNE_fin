#!/usr/bin/env python3

import argparse
import csv
import math
import os
from collections import defaultdict

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def revcomp(seq):
    return seq.translate(DNA_COMP)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description="Scan FASTA sequences for TF binding sites using SwissRegulon-style weight matrices.")
    p.add_argument("--fasta", required=True, help="Input FASTA")
    p.add_argument("--wm", required=True, help="TF weight matrix file")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--pseudocount", type=float, default=0.1)
    p.add_argument("--background", type=float, nargs=4, default=[0.25, 0.25, 0.25, 0.25], metavar=("A", "C", "G", "T"))
    p.add_argument("--min-relative-score", type=float, default=0.85)
    p.add_argument("--min-abs-score", type=float, default=None)
    p.add_argument("--overlap-mode", choices=["global", "per_tf"], default="global")
    return p.parse_args()


def read_fasta(path):
    records = []
    current_id = None
    current_desc = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append({"seq_id": current_id, "description": current_desc, "seq": "".join(seq_chunks).upper()})
                header = line[1:].strip()
                current_desc = header
                current_id = header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if current_id is not None:
        records.append({"seq_id": current_id, "description": current_desc, "seq": "".join(seq_chunks).upper()})
    if not records:
        raise ValueError("No FASTA records found")
    return records


def parse_weight_matrices(path):
    motifs = {}
    current_name = None
    current_rows = []
    with open(path, "r", encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            line = raw.strip()
            if not line:
                continue
            if line.startswith("//"):
                if current_name is not None and current_rows:
                    motifs[current_name] = current_rows
                current_name = None
                current_rows = []
                continue
            if line.startswith("NA"):
                if current_name is not None and current_rows:
                    motifs[current_name] = current_rows
                parts = line.split(maxsplit=1)
                if len(parts) < 2:
                    raise ValueError(f"Invalid NA row at line {lineno}: {raw}")
                current_name = parts[1].strip()
                current_rows = []
                continue
            if line.startswith("P0"):
                continue
            parts = line.split()
            if current_name is None:
                continue
            if len(parts) < 5:
                raise ValueError(f"Matrix row too short at line {lineno}: {raw}")
            if not parts[0].isdigit():
                continue
            try:
                current_rows.append([float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])])
            except ValueError:
                raise ValueError(f"Non-numeric matrix row at line {lineno}: {raw}")
    if current_name is not None and current_rows:
        motifs[current_name] = current_rows
    if not motifs:
        raise ValueError("No usable motifs found")
    return motifs


def counts_to_pwm(count_matrix, pseudocount, bg):
    pwm = []
    min_score = 0.0
    max_score = 0.0
    for row in count_matrix:
        total = sum(row) + 4.0 * pseudocount
        probs = [(x + pseudocount) / total for x in row]
        scores = [math.log2(probs[i] / bg[i]) for i in range(4)]
        pwm.append(scores)
        min_score += min(scores)
        max_score += max(scores)
    return pwm, min_score, max_score


def score_kmer(kmer, pwm):
    score = 0.0
    for i, base in enumerate(kmer):
        idx = BASE_INDEX.get(base)
        if idx is None:
            return None
        score += pwm[i][idx]
    return score


def relative_score(raw_score, min_score, max_score):
    return 0.0 if max_score == min_score else (raw_score - min_score) / (max_score - min_score)


def scan_sequence(seq, motif_name, pwm, min_score, max_score, min_relative_score, min_abs_score):
    hits = []
    L = len(pwm)
    n = len(seq)
    if n < L:
        return hits
    for i in range(0, n - L + 1):
        kmer = seq[i:i + L]
        raw = score_kmer(kmer, pwm)
        if raw is not None:
            rel = relative_score(raw, min_score, max_score)
            if rel >= min_relative_score and (min_abs_score is None or raw >= min_abs_score):
                hits.append({"motif": motif_name, "strand": "+", "start0": i, "end0": i + L, "matched_seq": kmer, "raw_score": raw, "rel_score": rel, "motif_len": L})
    rc_seq = revcomp(seq)
    for i in range(0, n - L + 1):
        rc_kmer = rc_seq[i:i + L]
        raw = score_kmer(rc_kmer, pwm)
        if raw is not None:
            rel = relative_score(raw, min_score, max_score)
            if rel >= min_relative_score and (min_abs_score is None or raw >= min_abs_score):
                orig_start = n - (i + L)
                orig_end = n - i
                hits.append({"motif": motif_name, "strand": "-", "start0": orig_start, "end0": orig_end, "matched_seq": revcomp(rc_kmer), "raw_score": raw, "rel_score": rel, "motif_len": L})
    return hits


def select_non_overlapping_hits(hits):
    hits_sorted = sorted(hits, key=lambda h: (-h["raw_score"], -h["rel_score"], -(h["end0"] - h["start0"]), h["start0"], h["motif"], h["strand"]))
    selected = []
    occupied = []
    for h in hits_sorted:
        if any(not (h["end0"] <= s or h["start0"] >= e) for s, e in occupied):
            continue
        selected.append(h)
        occupied.append((h["start0"], h["end0"]))
    selected.sort(key=lambda h: (h["start0"], h["end0"], h["motif"], h["strand"]))
    return selected


def select_hits_by_mode(hits, mode):
    if mode == "global":
        return select_non_overlapping_hits(hits)
    by_tf = defaultdict(list)
    for h in hits:
        by_tf[h["motif"]].append(h)
    out = []
    for sub in by_tf.values():
        out.extend(select_non_overlapping_hits(sub))
    out.sort(key=lambda h: (h["start0"], h["end0"], h["motif"], h["strand"]))
    return out


def main():
    args = parse_args()
    bg = args.background
    if len(bg) != 4 or any(x <= 0 for x in bg) or abs(sum(bg) - 1.0) > 1e-6:
        raise ValueError("Background frequencies must be four positive values summing to 1")
    fasta_records = read_fasta(args.fasta)
    motifs_raw = parse_weight_matrices(args.wm)
    motifs_pwm = {}
    for name, rows in motifs_raw.items():
        pwm, min_score, max_score = counts_to_pwm(rows, args.pseudocount, bg)
        motifs_pwm[name] = {"pwm": pwm, "min_score": min_score, "max_score": max_score, "len": len(pwm)}

    hits_out = args.out_prefix + ".tf_hits.tsv"
    summary_out = args.out_prefix + ".summary.tsv"
    motifs_out = args.out_prefix + ".motif_summary.tsv"
    os.makedirs(os.path.dirname(hits_out) or ".", exist_ok=True)

    motif_selected_counter = defaultdict(int)
    motif_candidate_counter = defaultdict(int)
    total_candidate_hits = 0
    total_selected_hits = 0
    with open(hits_out, "w", encoding="utf-8", newline="") as hout, open(summary_out, "w", encoding="utf-8", newline="") as sout:
        hw = csv.writer(hout, delimiter="\t")
        sw = csv.writer(sout, delimiter="\t")
        hw.writerow(["cne_id", "seq_id", "seq_length", "tf", "strand", "motif_len", "hit_start_in_cne_0b", "hit_end_in_cne_0b", "raw_score", "relative_score", "matched_seq"])
        sw.writerow(["cne_id", "seq_id", "seq_length", "candidate_hits", "selected_nonoverlapping_hits"])
        for rec in fasta_records:
            seq_id = rec["seq_id"]
            seq = rec["seq"]
            all_hits = []
            for motif_name, m in motifs_pwm.items():
                motif_hits = scan_sequence(seq, motif_name, m["pwm"], m["min_score"], m["max_score"], args.min_relative_score, args.min_abs_score)
                motif_candidate_counter[motif_name] += len(motif_hits)
                all_hits.extend(motif_hits)
            total_candidate_hits += len(all_hits)
            selected = select_hits_by_mode(all_hits, args.overlap_mode)
            total_selected_hits += len(selected)
            for h in selected:
                motif_selected_counter[h["motif"]] += 1
                hw.writerow([seq_id, seq_id, len(seq), h["motif"], h["strand"], h["motif_len"], h["start0"], h["end0"], f"{h['raw_score']:.6f}", f"{h['rel_score']:.6f}", h["matched_seq"]])
            sw.writerow([seq_id, seq_id, len(seq), len(all_hits), len(selected)])
    with open(motifs_out, "w", encoding="utf-8", newline="") as mout:
        mw = csv.writer(mout, delimiter="\t")
        mw.writerow(["tf", "motif_len", "candidate_hits", "selected_hits"])
        for motif_name in sorted(motifs_pwm):
            mw.writerow([motif_name, motifs_pwm[motif_name]["len"], motif_candidate_counter.get(motif_name, 0), motif_selected_counter.get(motif_name, 0)])
    print("Done.")
    print(f"FASTA records: {len(fasta_records)}")
    print(f"Motifs: {len(motifs_pwm)}")
    print(f"Candidate hits: {total_candidate_hits}")
    print(f"Selected non-overlapping hits: {total_selected_hits}")
    print(f"Output: {hits_out}")
    print(f"Output: {summary_out}")
    print(f"Output: {motifs_out}")


if __name__ == "__main__":
    main()
