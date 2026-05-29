#!/usr/bin/env python3

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.stderr.write("Missing pysam. Install with: pip install pysam or conda install -c bioconda pysam\n")
    sys.exit(1)

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def revcomp(seq):
    return seq.translate(DNA_COMP)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description="Direct TFBS scanning on CNE coordinates using indexed genome FASTA and weight matrices.")
    p.add_argument("--csv", required=True, help="Input CNE table, CSV or TSV")
    p.add_argument("--genome", required=True, help="Genome FASTA indexed with .fai")
    p.add_argument("--wm", required=True, help="TF weight matrix text file")
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    p.add_argument("--pseudocount", type=float, default=0.1)
    p.add_argument("--background", type=float, nargs=4, default=[0.25, 0.25, 0.25, 0.25], metavar=("A", "C", "G", "T"))
    p.add_argument("--min-relative-score", type=float, default=0.85)
    p.add_argument("--min-abs-score", type=float, default=None)
    p.add_argument("--overlap-mode", choices=["global", "per_tf"], default="global")
    p.add_argument("--keep-duplicates", action="store_true")
    p.add_argument("--delimiter", choices=["auto", "tab", "comma"], default="auto")
    return p.parse_args()


def detect_delimiter(path, mode="auto"):
    if mode == "tab":
        return "\t"
    if mode == "comma":
        return ","
    with open(path, "r", encoding="utf-8-sig") as f:
        line = f.readline()
    return "\t" if line.count("\t") >= line.count(",") else ","


def ensure_fai(genome_fasta):
    if not os.path.exists(genome_fasta + ".fai"):
        raise FileNotFoundError(f"Missing FASTA index: {genome_fasta}.fai. Run: samtools faidx {genome_fasta}")


def make_cne_id(chrom, start, end):
    return f"{chrom}:{int(start)}-{int(end)}"


def parse_cne_table(path, keep_duplicates=False, delimiter_mode="auto"):
    delimiter = detect_delimiter(path, delimiter_mode)
    required = {"seqnames", "start", "end"}
    records = []
    seen = set()
    duplicate_count = 0
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("Input table has no readable header")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        for i, row in enumerate(reader, start=1):
            chrom = row["seqnames"].strip()
            if not chrom:
                continue
            try:
                start_1b = int(float(row["start"]))
                end_1b = int(float(row["end"]))
            except ValueError:
                raise ValueError(f"Invalid coordinates at row {i}: start={row.get('start')}, end={row.get('end')}")
            if end_1b < start_1b:
                start_1b, end_1b = end_1b, start_1b
            key = (chrom, start_1b, end_1b)
            if not keep_duplicates and key in seen:
                duplicate_count += 1
                continue
            seen.add(key)
            records.append({
                "cne_id": make_cne_id(chrom, start_1b, end_1b),
                "chrom": chrom,
                "start_1b": start_1b,
                "end_1b": end_1b,
                "width_csv": row.get("width", ""),
                "strand_csv": row.get("strand", ""),
                "phastcons": row.get("phastcons", ""),
                "annotation": row.get("annotation", ""),
                "gene_name": row.get("gene_name", ""),
                "gene_id": row.get("geneId", ""),
                "transcript_id": row.get("transcriptId", ""),
                "distance_to_tss": row.get("distanceToTSS", ""),
                "in_atac_peak": row.get("in_atac_peak", ""),
                "nearby_gene_active": row.get("nearby_gene_active", ""),
                "nearby_gene_fin_dev": row.get("nearby_gene_fin_dev", ""),
                "source_row_number": i,
            })
    if not records:
        raise ValueError("No usable CNE records remained")
    return records, duplicate_count, delimiter


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
        raise ValueError("No usable motifs found in weight matrix file")
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
    ensure_fai(args.genome)
    cnes, duplicate_count, delimiter = parse_cne_table(args.csv, args.keep_duplicates, args.delimiter)
    motifs_raw = parse_weight_matrices(args.wm)
    motifs_pwm = {name: {"pwm": counts_to_pwm(rows, args.pseudocount, bg)[0], "min_score": counts_to_pwm(rows, args.pseudocount, bg)[1], "max_score": counts_to_pwm(rows, args.pseudocount, bg)[2], "len": len(rows)} for name, rows in motifs_raw.items()}

    hits_out = args.out_prefix + ".tf_hits.tsv"
    summary_out = args.out_prefix + ".summary.tsv"
    motifs_out = args.out_prefix + ".motif_summary.tsv"
    os.makedirs(os.path.dirname(hits_out) or ".", exist_ok=True)

    motif_selected_counter = defaultdict(int)
    motif_candidate_counter = defaultdict(int)
    total_candidate_hits = 0
    total_selected_hits = 0
    fasta = pysam.FastaFile(args.genome)
    with open(hits_out, "w", encoding="utf-8", newline="") as hout, open(summary_out, "w", encoding="utf-8", newline="") as sout:
        hw = csv.writer(hout, delimiter="\t")
        sw = csv.writer(sout, delimiter="\t")
        hw.writerow(["cne_id", "chrom", "cne_start_1b", "cne_end_1b", "cne_length", "tf", "strand", "motif_len", "hit_start_in_cne_0b", "hit_end_in_cne_0b", "genomic_start_1b", "genomic_end_1b", "raw_score", "relative_score", "matched_seq", "annotation", "gene_name", "gene_id", "transcript_id", "distance_to_tss", "in_atac_peak", "nearby_gene_active", "nearby_gene_fin_dev"])
        sw.writerow(["cne_id", "chrom", "cne_start_1b", "cne_end_1b", "cne_length", "candidate_hits", "selected_nonoverlapping_hits", "annotation", "gene_name", "gene_id", "in_atac_peak", "nearby_gene_active", "nearby_gene_fin_dev"])
        for rec in cnes:
            seq = fasta.fetch(rec["chrom"], rec["start_1b"] - 1, rec["end_1b"]).upper()
            if not seq:
                continue
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
                genomic_start_1b = rec["start_1b"] + h["start0"]
                genomic_end_1b = rec["start_1b"] + h["end0"] - 1
                hw.writerow([rec["cne_id"], rec["chrom"], rec["start_1b"], rec["end_1b"], len(seq), h["motif"], h["strand"], h["motif_len"], h["start0"], h["end0"], genomic_start_1b, genomic_end_1b, f"{h['raw_score']:.6f}", f"{h['rel_score']:.6f}", h["matched_seq"], rec["annotation"], rec["gene_name"], rec["gene_id"], rec["transcript_id"], rec["distance_to_tss"], rec["in_atac_peak"], rec["nearby_gene_active"], rec["nearby_gene_fin_dev"]])
            sw.writerow([rec["cne_id"], rec["chrom"], rec["start_1b"], rec["end_1b"], len(seq), len(all_hits), len(selected), rec["annotation"], rec["gene_name"], rec["gene_id"], rec["in_atac_peak"], rec["nearby_gene_active"], rec["nearby_gene_fin_dev"]])
    fasta.close()
    with open(motifs_out, "w", encoding="utf-8", newline="") as mout:
        mw = csv.writer(mout, delimiter="\t")
        mw.writerow(["tf", "motif_len", "candidate_hits", "selected_hits"])
        for motif_name in sorted(motifs_pwm):
            mw.writerow([motif_name, motifs_pwm[motif_name]["len"], motif_candidate_counter.get(motif_name, 0), motif_selected_counter.get(motif_name, 0)])
    print("Done.")
    print(f"Input delimiter: {'TAB' if delimiter == chr(9) else 'COMMA'}")
    print(f"CNE records: {len(cnes)}")
    print(f"Duplicate records skipped: {duplicate_count if not args.keep_duplicates else 0}")
    print(f"Motifs: {len(motifs_pwm)}")
    print(f"Candidate hits: {total_candidate_hits}")
    print(f"Selected non-overlapping hits: {total_selected_hits}")
    print(f"Output: {hits_out}")
    print(f"Output: {summary_out}")
    print(f"Output: {motifs_out}")


if __name__ == "__main__":
    main()
