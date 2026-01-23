# CNE_fin

Summary
-------
This repository contains analyses and results for comparative non-coding element (CNE) discovery and motif enrichment used in the CNE_fin project. It includes input data (genome and annotation fragments, peak sets, PWM files), analysis scripts (R and shell), and generated outputs (tables, motif analysis results and a Shiny app for exploration).

Repository layout
-----------------
- [input](input): original and processed input files (BED/TSV/FASTA/PWM).
- [ancilliary_files](ancilliary_files): large reference genomes, chains and indexes (kept out of git by default).
- [output](output): generated results, tables and MEME/AME outputs.
- [scripts](scripts): R scripts and shell wrappers used to run analyses.
- [shiny_app](shiny_app): small Shiny app and Docker container files to explore results.

Key files
---------
- [scripts/cne_analysis.R](scripts/cne_analysis.R) — main analysis driver R script.
- [scripts/custom_functions.R](scripts/custom_functions.R) — helper functions used by analysis scripts.
- [scripts/meme_motif_enrichment.sh](scripts/meme_motif_enrichment.sh) — wrapper calling MEME/AME for motif enrichment.
- [shiny_app/app.R](shiny_app/app.R) — Shiny app for result exploration.

Structured list of outputs
-------------------------
```
output/
├─ cne_atac_overlaps_by_annotation.png
├─ cne_atac_overlaps_fin_specific_by_annotation.png
├─ cne_atac_overlaps.tsv
├─ cne_widths.png
├─ comb_heatmap.pdf
├─ comb_profile.pdf
├─ GO_enrichment_results.csv
├─ GO_enrichment_results.tsv
├─ GO-filtered_cne_annotations.pdf
├─ teleost_cne_atac_overlaps_annotated.tsv
├─ teleost_cne_enrichment_dotplot.pdf
├─ teleost_GO-fin_specific_cne.csv
├─ teleost_GO-fin_specific_cne.tsv
├─ teleost_GO-fin-skeletal_specific_cne.csv
├─ teleost_GO-fin-skeletal_specific_cne.tsv
├─ teleost_hits_with_atac_annotated.tsv
├─ teleost_specific_cne.csv
├─ teleost_specific_cne.tsv
├─ teleost_unique_cne_with_atac_annotated_active_genes.tsv
├─ teleost_unique_cne_with_atac_annotated.tsv
├─ upset_overlaps.pdf
└─ meme_analysis/
  ├─ motif_enrichment_summary.txt
  ├─ teleostei_motifs/
  │  ├─ teleostei_cne.fa
  │  └─ ame_output/
  │     ├─ ame.html
  │     ├─ ame.tsv
  │     └─ sequences.tsv          # derived sequence output (large; excluded from repo by default)
  └─ vertebrata_motifs/
    ├─ vertebrata_cne.fa
    └─ ame_output/
      ├─ ame.html
      ├─ ame.tsv
      ├─ ame.tsv.1
      └─ sequences.tsv          # derived sequence output (large; excluded from repo by default)
```

Notes:
- The files above are the generated outputs present in `output/`. Most are small analysis result tables and figures and are intended to be tracked in the repository. Large derived sequence tables under `output/meme_analysis/*/ame_output/sequences.tsv` remain excluded by default — use Git LFS or external archives if you want to version them.

Update about scripts
--------------------
- I will remove any deleted script from this README when you tell me which script you removed. Current scripts in `scripts/` are `cne_analysis.R`, `custom_functions.R`, and `meme_motif_enrichment.sh`.

Notes and recommendations
-------------------------
- Raw genome FASTA and large annotation files (for example [ancilliary_files/drer.fa](ancilliary_files/drer.fa) and [ancilliary_files/drer.gff](ancilliary_files/drer.gff)) are large and are excluded by the repository's `.gitignore`
- Derived large outputs under `output/meme_analysis/` are also excluded by default
