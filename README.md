# CNE_fin

Summary
-------
This repository contains analyses and results for comparative non-coding element (CNE) discovery and motif enrichment. It includes input data (genome and annotation fragments, peak sets, PWM files), analysis scripts (R and shell), and generated outputs (tables, motif analysis results and a Shiny app for exploration).

Definitions: Teleost vs Vertebrate CNEs
-------------------------------------
- Teleost CNEs: CNEs which are found in teleost species but not in tetrapods (teleost-specific).
- Vertebrate CNEs: CNEs which are shared by teleost species and also tetrapod species (conserved across vertebrates).

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
- [shiny_app/app.R](shiny_app/app.R) — Shiny app for result exploration.

List of outputs
-------------------------
```
output/
├─ cne_atac_overlaps_by_annotation.png                  # counts of CNE–ATAC overlaps by genomic annotation (promoter/enhancer/etc.).
├─ cne_atac_overlaps_fin_specific_by_annotation.png     # same plot restricted to fin-specific CNEs.
├─ cne_atac_overlaps.tsv                                # table listing overlaps between CNEs and ATAC peaks with annotation fields.
├─ cne_widths.png                                       # distribution plot of CNE lengths.
├─ comb_heatmap.pdf                                     # heatmap of ATAC-seq signal in Teleost and Vertebrata CNEs.
├─ comb_profile.pdf                                     # averaged ATAC-seq profile across Teleost and Vertebrata CNEs.
├─ GO_enrichment_results.csv                            # GO enrichment results for CNE-associated gene sets (CSV).
├─ GO_enrichment_results.tsv                            # GO enrichment results (TSV).
├─ GO-filtered_cne_annotations.pdf                      # PDF summarizing GO-filtered annotations for candidate CNE genes.
├─ teleost_cne_atac_overlaps_annotated.tsv              # annotated overlaps for teleost CNEs and ATAC peaks.
├─ teleost_cne_enrichment_dotplot.pdf                   # dotplot of enrichment statistics for teleost CNE-associated terms.
├─ teleost_GO-fin_specific_cne.csv                      # GO enrichment (CSV) for fin-specific teleost CNE-associated genes.
├─ teleost_GO-fin_specific_cne.tsv                      # GO enrichment (TSV) for fin-specific teleost CNE-associated genes.
├─ teleost_GO-fin-skeletal_specific_cne.csv             # GO enrichment (CSV) for skeletal-specific teleost CNE-associated genes.
├─ teleost_GO-fin-skeletal_specific_cne.tsv             # GO enrichment (TSV) for skeletal-specific teleost CNE-associated genes.
├─ teleost_hits_with_atac_annotated.tsv                 # teleost CNEs overlapping ATAC peaks with gene and peak metadata.
├─ teleost_specific_cne.csv                             # list of teleost-specific CNEs (CSV).
├─ teleost_specific_cne.tsv                             # list of teleost-specific CNEs (TSV).
├─ teleost_unique_cne_with_atac_annotated_active_genes.tsv  # teleost-unique CNEs annotated with active gene assignments.
├─ teleost_unique_cne_with_atac_annotated.tsv           # teleost-unique CNEs with ATAC overlap annotations.
├─ upset_overlaps.pdf                                   # UpSet plot summarizing overlaps between different CNE and peak sets.
└─ meme_analysis/
  ├─ teleostei_motifs/
  │  ├─ teleostei_cne.fa                                # FASTA of teleost CNE sequences used for motif discovery/enrichment.
  │  └─ ame_output/
  │     ├─ ame.html                                     # AME HTML report summarizing motif enrichment for teleost CNEs.
  │     ├─ ame.tsv                                      # Tabular AME results.
  │     └─ sequences.tsv                                # derived sequence output (large; excluded from repo by default).
  └─ vertebrata_motifs/
    ├─ vertebrata_cne.fa                                # FASTA of vertebrate-shared CNE sequences used for motif analyses.
    └─ ame_output/
      ├─ ame.html                                       # AME HTML report summarizing motif enrichment for vertebrate CNEs.
      ├─ ame.tsv                                        # Tabular AME results.
      └─ sequences.tsv                                  # derived sequence output (large; excluded from repo by default).
```

> **Note**
>The files above are the generated outputs present in `output/`. Most are small analysis result tables and figures and are intended to be tracked in the repository. Large derived sequence tables under `output/meme_analysis/*/ame_output/sequences.tsv` remain excluded by default.
