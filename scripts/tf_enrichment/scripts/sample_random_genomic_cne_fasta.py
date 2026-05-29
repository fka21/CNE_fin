#!/usr/bin/env python3

import argparse
import random
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.stderr.write("Missing pysam. Install with: pip install pysam or conda install -c bioconda pysam\n")
    sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser(description="Sample length-matched random genomic intervals from an indexed genome FASTA.")
    p.add_argument("--genome", required=True, help="Indexed genome FASTA")
    p.add_argument("--lengths", required=True, help="One sequence length per line")
    p.add_argument("--out", required=True, help="Output FASTA")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--max-n-frac", type=float, default=0.0, help="Maximum allowed fraction of N bases")
    p.add_argument("--max-attempts-per-seq", type=int, default=1000)
    return p.parse_args()


def read_lengths(path):
    lengths = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                lengths.append(int(line))
    if not lengths:
        raise ValueError("Empty length list")
    return lengths


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def choose_contig(references, ref_lengths, requested_len):
    valid = [(r, l) for r, l in zip(references, ref_lengths) if l >= requested_len]
    if not valid:
        raise ValueError(f"No contig is long enough for requested length: {requested_len}")
    names = [x[0] for x in valid]
    weights = [x[1] - requested_len + 1 for x in valid]
    return random.choices(names, weights=weights, k=1)[0]


def n_fraction(seq):
    return 1.0 if not seq else seq.upper().count("N") / len(seq)


def main():
    args = parse_args()
    random.seed(args.seed)
    if not Path(args.genome + ".fai").exists():
        raise FileNotFoundError(f"Missing FASTA index: {args.genome}.fai. Run: samtools faidx {args.genome}")

    lengths = read_lengths(args.lengths)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fasta = pysam.FastaFile(args.genome)
    references = list(fasta.references)
    ref_lengths = list(fasta.lengths)

    written = 0
    skipped = 0
    with open(args.out, "w") as out:
        for i, length in enumerate(lengths, start=1):
            success = False
            for _ in range(args.max_attempts_per_seq):
                chrom = choose_contig(references, ref_lengths, length)
                chrom_len = fasta.get_reference_length(chrom)
                start0 = random.randint(0, chrom_len - length)
                end0 = start0 + length
                seq = fasta.fetch(chrom, start0, end0).upper()
                if len(seq) != length:
                    continue
                if n_fraction(seq) > args.max_n_frac:
                    continue
                start1 = start0 + 1
                end1 = end0
                out.write(f">random_genomic_CNE_{i:05d}|{chrom}:{start1}-{end1}|length={length}|seed={args.seed}\n")
                out.write(wrap(seq) + "\n")
                written += 1
                success = True
                break
            if not success:
                skipped += 1
                sys.stderr.write(f"Warning: no suitable interval found for index={i}, length={length}\n")

    fasta.close()
    print(f"Done: {args.out}")
    print(f"Requested sequences: {len(lengths)}")
    print(f"Written sequences: {written}")
    print(f"Skipped sequences: {skipped}")
    print(f"Min length: {min(lengths)}")
    print(f"Max length: {max(lengths)}")
    print(f"Mean length: {sum(lengths) / len(lengths):.3f}")
    print(f"Max N fraction: {args.max_n_frac}")


if __name__ == "__main__":
    main()
