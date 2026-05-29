# Workflow commands

Run separately for each dataset by changing the dataset name and input TSV.

Example for actinopterygii:

```bash
DATASET=actinopteriigy
INPUT_TSV=actinopteriigy_cne_final_table.tsv
GENOME=GCF_049306965.1_GRCz12tu_genomic.fna
WM=TF_weight_matrix.txt
OUTDIR=${DATASET}

mkdir -p ${OUTDIR}/{01_length_stats,02_real_CNE_scan,04_random_genomic_background,05_random_genomic_scan,06_enrichment_real_vs_random_genomic,logs}

python3 scripts/extract_unique_cne_lengths.py \
  --cne-table ${INPUT_TSV} \
  --out-prefix ${OUTDIR}/01_length_stats/${DATASET}_unique_CNE

python3 scripts/cne_tf_scan_direct.py \
  --csv ${INPUT_TSV} \
  --genome ${GENOME} \
  --wm ${WM} \
  --out-prefix ${OUTDIR}/02_real_CNE_scan/${DATASET}_cne_tf \
  --min-relative-score 0.85 \
  --overlap-mode global

python3 scripts/sample_random_genomic_cne_fasta.py \
  --genome ${GENOME} \
  --lengths ${OUTDIR}/01_length_stats/${DATASET}_unique_CNE_lengths.txt \
  --out ${OUTDIR}/04_random_genomic_background/${DATASET}_random_genomic_lengthmatched.fa \
  --seed 42

python3 scripts/fasta_tf_scan.py \
  --fasta ${OUTDIR}/04_random_genomic_background/${DATASET}_random_genomic_lengthmatched.fa \
  --wm ${WM} \
  --out-prefix ${OUTDIR}/05_random_genomic_scan/${DATASET}_random_genomic_cne_tf \
  --min-relative-score 0.85 \
  --overlap-mode global

python3 scripts/compare_tf_enrichment_real_vs_random.py \
  --real-hits ${OUTDIR}/02_real_CNE_scan/${DATASET}_cne_tf.tf_hits.tsv \
  --real-summary ${OUTDIR}/02_real_CNE_scan/${DATASET}_cne_tf.summary.tsv \
  --bg-hits ${OUTDIR}/05_random_genomic_scan/${DATASET}_random_genomic_cne_tf.tf_hits.tsv \
  --bg-summary ${OUTDIR}/05_random_genomic_scan/${DATASET}_random_genomic_cne_tf.summary.tsv \
  --out-prefix ${OUTDIR}/06_enrichment_real_vs_random_genomic/${DATASET}_vs_random_genomic \
  --min-score 0.85
```

For gnathostomata, set:

```bash
DATASET=gnathostomata
INPUT_TSV=gnathostomata_cne_final_table.tsv
```
