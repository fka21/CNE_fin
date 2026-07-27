# CNE_fin

Comparative analysis of conserved non-coding elements (CNEs) across vertebrates, identifying lineage-specific regulatory elements and their associations with gene function and chromatin accessibility in *Danio rerio*.

## Definitions

- **Actinopterygii-specific CNEs** — conserved in ray-finned fishes but absent in tetrapods; coordinates in GRCz11.
- **Gnathostomata-conserved CNEs** — conserved across jawed vertebrates; coordinates in GRCz11 (zebrafish reference) and GRCh38 (human reference).
- **Sarcopterygii-specific CNEs** — conserved in lobe-finned vertebrates but absent in the broader gnathostome set; coordinates in GRCh38.

## Repository layout

```
input/                        CNE BED files, ATAC peaks, RNA-seq TPM, external datasets
ancilliary_files/             Reference annotation (GFF, chrom sizes, sequence reports)
workflow/
  Snakefile                   End-to-end CNE discovery pipeline (Cactus → phastCons → classify)
  config/
    config.yml                Species, tree strings, paths, phastCons parameters
    cactus.config             Progressive Cactus alignment config
  envs/
    maf_manipulation.yml      Conda env: maf_sort, msa_view, phyloFit
    bedtools.yml              Conda env: bedtools
scripts/
  custom_functions.R          Shared helper functions (sourced by all R scripts)
  01_cne_preprocess.R         CNE import, ChIPseeker annotation, gene activity
  02_go_enrichment.R          rGREAT GO enrichment, hotspot scoring, TSS diagnostics
  03_post_process.R           Filtered exports, UpSet plots, ATAC overlaps
  04_visualizations.R         Circos plots, Gviz locus views, CNE width distributions
  figure_generators/
    Kagan2026-Figure_2.R      CNE set comparison (counts, widths, hotspots)
    Kagan2026-Figure_3.R      GO enrichment network panels
    Kagan2026-Figure_4.R      TF motif enrichment volcano panels
    Kagan2026-Supplementary_Figure_1.R   pcdh cluster locus views
    Kagan2026-Supplementary_Figure_5.R   CNE abundance vs gene structure
    Kagan2026-Supplementary_Figure_6.R   Top-N hotspot summary with GO annotation
    Kagan2026-Supplementary_Figures_7_8.R  GO sub-tree panels by biological process
  tf_enrichment/              Shell/Python scripts for AME motif enrichment
output/                       Generated figures and tables (see below)
```

## Reproducing the analysis

### Step 1 — CNE discovery (Snakemake)

Run from the repository root. Requires Singularity/Apptainer for the Cactus container and Conda for the MAF and bedtools environments.

```bash
snakemake --cores <N> --use-conda --use-singularity
```

Outputs four classified BED files under `classified/` that feed into Step 2.

### Step 2 — R analysis pipeline

Run in order from the `scripts/` directory.

```bash
Rscript 01_cne_preprocess.R
Rscript 02_go_enrichment.R
Rscript 03_post_process.R
Rscript 04_visualizations.R
```

The sarcopterygii / human-reference panels require three additional files in `ancilliary_files/` and `input/`; after running `01` check `output/reports/` for a summary of what was found.

### Step 3 — Figure scripts

Each script in `scripts/figure_generators/` is self-contained and reads from `output/`. Run from that directory or set paths at the top of each file.

```bash
cd output
Rscript ../scripts/figure_generators/Kagan2026-Figure_2.R
# etc.
```

TF enrichment scripts in `scripts/tf_enrichment/` are run independently; see the README in that subdirectory.

## Key outputs

```
output/
├── preprocessed/                          Serialised R objects shared between analysis scripts
├── great/
│   ├── great_*_GObp*.tsv                  rGREAT GO:BP enrichment (whole-genome and CNE-union backgrounds)
│   ├── great_strict_dualtest*.tsv         Dual-test (binomial + hypergeometric) significant terms
│   ├── cne_hotspots_domain_normalised.tsv Per-gene CNE density normalised to regulatory domain length
│   └── shared_fc_GO_great.tsv             Fold-enrichment comparison across CNE sets
├── reports/
│   ├── symbol_mapping_rate.tsv            SYMBOL match rate between GFF annotation and org.Dr.eg.db
│   ├── tss_distance_summary.tsv           Signed TSS-distance summary (rGREAT region–gene associations)
│   ├── annotation_composition*.tsv        ChIPseeker annotation class breakdown
│   └── hotspot_rank_raw_vs_normalised.tsv Effect of domain-length normalisation on hotspot ranking
├── circos_element_density.pdf             CNE density, zebrafish chromosomes (GRCz11)
├── circos_element_density_hsap.pdf        CNE density, human chromosomes (GRCh38)
├── cne_hotspots_domain_normalised.pdf     Hotspot loci ranked by obs/expected per regulatory domain
├── cne_annotations.pdf                    ChIPseeker annotation bar chart
├── power_law_like_fits*.pdf               CNE width distributions
├── upset_overlaps_*.pdf                   UpSet plots (ATAC, active genes, Chan enhancers, YueSong CNEs)
├── gviz_GO0007224/                        Gviz locus panels for smoothened-signalling genes
└── meme_analysis/                         AME motif enrichment results
```

## Dependencies

**Workflow:** Cactus v2.9.9 (container), maf_sort, msa_view, phyloFit (PHAST), bedtools.

**R packages:** `GenomicRanges`, `GenomicFeatures`, `rtracklayer`, `ChIPseeker`, `rGREAT` (≥ 2.0), `org.Dr.eg.db`, `simplifyEnrichment`, `ComplexHeatmap`, `circlize`, `CNEr`, `Gviz`, `zFPKM`, `SummarizedExperiment`, `GO.db`, `igraph`, `ggraph`, `tidygraph`, `tidyverse`, `tidyplots`, `readxl`, `poweRlaw`, `patchwork`, `ggrepel`, `ggpubr`.
