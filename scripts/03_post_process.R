### LIBRARIES ###
pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  readxl,
  ComplexHeatmap,
  circlize,
  GenomicFeatures,
  patchwork
)

setwd(this.path::here())
source("custom_functions.R")

preproc_dir <- "../output/preprocessed"
great_dir <- "../output/great"
report_dir <- "../output/reports"
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

### --- Load preprocessed objects -------------------------------------------
actinopteriigy_cne_gr <- readRDS(
  file.path(preproc_dir, "actinopteriigy_cne_gr.rds")
)
gnathostomata_cne_gr <- readRDS(
  file.path(preproc_dir, "gnathostomata_cne_gr.rds")
)
peak_anno_list <- readRDS(file.path(preproc_dir, "peak_anno_list.rds"))
gene_activity <- readRDS(file.path(preproc_dir, "gene_activity.rds"))

# Built once in script 1 and reused here rather than re-read with a different
# coordinate convention.
atac_peaks_gr <- readRDS(file.path(preproc_dir, "atac_peaks_gr.rds"))
enh_gr <- readRDS(file.path(preproc_dir, "chan_enhancers_gr.rds"))
actinopteriigy_cne_ov <- readRDS(
  file.path(preproc_dir, "actinopteriigy_cne_ov.rds")
)

### --- Load rGREAT-derived gene lists --------------------------------------
fg <- readRDS(file.path(great_dir, "filtered_genes_for_downstream.rds"))
filtered_genes_fin <- fg$fin
filtered_genes_skeletal_fin <- fg$fin_skeletal
filtered_genes <- fg$all_signif_actino # for the "any signif GO" filters

great_tbl_actino <- read_tsv(
  file.path(great_dir, "great_actinopterygii_GObp.tsv"),
  show_col_types = FALSE
)
hotspots <- readRDS(file.path(great_dir, "cne_hotspots.rds"))

drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)

### --- YueSong overlap set --------------------------------------------------
# Coordinates come straight from the liftover overlap computed in script 1; the
# earlier reconstruction of cne_name by pasting ".1" onto the chromosome assumed
# every accession carried version 1.
yuesong_gr <- GRanges(
  seqnames = sub("\\.[0-9]+$", "", actinopteriigy_cne_ov$cne_name),
  ranges = IRanges(
    start = actinopteriigy_cne_ov$start,
    end = actinopteriigy_cne_ov$end
  )
)

### --- GO-derived gene filtering on peak annotations -----------------------
anno_actino <- non_exon(peak_anno_list$actinopteriigy_CNE@anno)
anno_gnatho <- non_exon(peak_anno_list$gnathostomata_CNE@anno)

actinopteriigy_filtered_skeletal_fin <- as.data.frame(
  anno_actino[
    (anno_actino$geneId %in% filtered_genes_skeletal_fin)
  ]
)
actinopteriigy_filtered_fin <- as.data.frame(
  anno_actino[
    (anno_actino$geneId %in% filtered_genes_fin)
  ]
)

write.table(
  actinopteriigy_filtered_fin,
  "../output/actinopteriigy_GO-fin_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  actinopteriigy_filtered_skeletal_fin,
  "../output/actinopteriigy_GO-fin-skeletal_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
### CNE HOTSPOTS — DOMAIN-NORMALISED FIGURE                                 ###
###############################################################################
# Bar height is observed/expected against the GREAT regulatory-domain length, so
# loci are ranked by CNE density rather than by size. Raw counts and domain
# extents are printed on each bar so the two quantities can be compared
# directly, and 3R ohnologue pair members are shaded separately because each
# paralogue carries its own domain.

top_n_hotspots <- 20

actino_col <- "#CC79A7"
gnatho_col <- "#009E73"

p_hotspots <- plot_hotspot_panel(
  hotspots,
  set_name = "actinopterygii",
  fill_colour = actino_col,
  top_n = top_n_hotspots,
  title = "Actinopterygii-specific"
) +
  plot_hotspot_panel(
    hotspots,
    set_name = "gnathostomata",
    fill_colour = gnatho_col,
    top_n = top_n_hotspots,
    title = "Gnathostomata-conserved"
  ) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Bar labels give observed CNEs / regulatory domain length.",
      "Dashed line marks the domain-length expectation (obs/exp = 1)."
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  "../output/cne_hotspots_domain_normalised.pdf",
  p_hotspots,
  width = 12,
  height = 7
)


hotspot_size_correlation <- hotspots %>%
  filter(observed > 0) %>%
  group_by(set) %>%
  summarise(
    spearman_count_vs_domain_length = round(
      cor(observed, domain_width, method = "spearman"),
      3
    ),
    spearman_obsexp_vs_domain_length = round(
      cor(obs_exp, domain_width, method = "spearman"),
      3
    ),
    .groups = "drop"
  )
write_tsv(
  hotspot_size_correlation,
  file.path(report_dir, "hotspot_size_bias_correlation.tsv")
)
print(hotspot_size_correlation)

### --- Active-gene filtered exports ----------------------------------------
export_active <- as.data.frame(
  anno_actino[
    which(anno_actino$is_active)
  ]
)
write.table(
  export_active,
  "../output/actinopteriigy_cne_near_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  filter(actinopteriigy_filtered_fin, is_active),
  "../output/actinopteriigy_GO-fin_specific_cne_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(actinopteriigy_filtered_skeletal_fin, is_active),
  "../output/actinopteriigy_GO-fin-skeletal_specific_cne_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(anno_gnatho),
  "../output/gnathostomata_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- ATAC-seq overlap ----------------------------------------------------
cne_atac_annotated <- subsetByOverlaps(anno_actino, atac_peaks_gr)
cne_atac_annotated_df <- as.data.frame(cne_atac_annotated)

write.table(
  cne_atac_annotated_df,
  "../output/actinopteriigy_cne_atac_overlaps_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(cne_atac_annotated_df, is_active),
  "../output/actinopteriigy_cne_atac_overlaps_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ATAC × annotation bar chart
overlap_by_annotation <- cne_atac_annotated_df |>
  as_tibble() |>
  mutate(annotation = str_remove(annotation, " \\(.*$")) |>
  dplyr::count(annotation) |>
  mutate(perc = n / sum(n) * 100) |>
  arrange(desc(n))

p <- ggplot(
  overlap_by_annotation,
  aes(reorder(annotation, -n), n, fill = annotation)
) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(
    aes(label = paste0(n, "\n(", round(perc, 1), "%)")),
    vjust = -0.5,
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Annotation Type",
    y = "Count",
    title = "Distribution of ATAC-overlapping CNEs by Annotation"
  )
ggsave(
  "../output/cne_atac_overlaps_by_annotation.pdf",
  p,
  width = 8,
  height = 6,
  dpi = 300
)

# fin-specific ATAC subsets
write.table(
  filter(cne_atac_annotated_df, geneId %in% filtered_genes_fin),
  "../output/actinopteriigy_cne_atac_overlaps_fin_specific.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(
    cne_atac_annotated_df,
    geneId %in% filtered_genes_fin,
    is_active
  ),
  "../output/actinopteriigy_cne_atac_overlaps_fin_specific_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- actinopteriigy-unique CNEs ∩ ATAC -----------------------------------
actinopteriigy_only_gr <- subsetByOverlaps(
  actinopteriigy_cne_gr,
  gnathostomata_cne_gr,
  invert = TRUE
)
tel_only_with_atac <- subsetByOverlaps(actinopteriigy_only_gr, atac_peaks_gr)
final_df_tel <- as.data.frame(
  subsetByOverlaps(anno_actino, tel_only_with_atac)
)

write.table(
  final_df_tel,
  "../output/actinopteriigy_unique_cne_with_atac_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(final_df_tel, is_active),
  "../output/actinopteriigy_unique_cne_with_atac_annotated_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- GRanges RDS exports for GenometriCorr etc. --------------------------
export_rds_dir <- "../output/granges_rds"
dir.create(export_rds_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  actinopteriigy_cne_gr,
  file.path(export_rds_dir, "actinopteriigy_cne_gr.rds")
)
saveRDS(
  gnathostomata_cne_gr,
  file.path(export_rds_dir, "gnathostomata_cne_gr.rds")
)
saveRDS(atac_peaks_gr, file.path(export_rds_dir, "atac_peaks_gr.rds"))
saveRDS(peak_anno_list, file.path(export_rds_dir, "peak_annotations.rds"))
saveRDS(genes(drer_anno), file.path(export_rds_dir, "genome_genes_gr.rds"))
saveRDS(
  rtracklayer::import("../ancilliary_files/drer.gff"),
  file.path(export_rds_dir, "genome_annotation_gr.rds")
)

### --- UpSet sets ----------------------------------------------------------
# The Chan enhancer and YueSong CNE sets are still computed and written to the
# final tables, but the UpSet panels are restricted to the accessibility and
# gene-activity sets: those are the two independent lines of evidence for
# regulatory potential, and the external CNE/enhancer sets overlap the universe
# so sparsely that they compress the intersection bars.

SLACK <- 0L # 0 for native coordinates, 50-200 bp for lifted-over sets

UPSET_SETS <- c("ATAC peaks", "Active genes nearby")

res_actino <- analyse_cne_universe(
  anno_gr = unique(anno_actino),
  label = "actinopteriigy",
  atac_peaks_gr = atac_peaks_gr,
  enh_gr = enh_gr,
  yuesong_gr = yuesong_gr,
  fin_geneIds = filtered_genes_fin,
  slack = SLACK,
  upset_sets = UPSET_SETS
)

res_gnatho <- analyse_cne_universe(
  anno_gr = unique(anno_gnatho),
  label = "gnathostomata",
  atac_peaks_gr = atac_peaks_gr,
  enh_gr = enh_gr,
  fin_geneIds = filtered_genes_fin,
  slack = SLACK,
  upset_sets = UPSET_SETS
)

### --- Final tables for the Shiny app --------------------------------------

write.table(
  res_actino$final,
  "../output/actinopteriigy_cne_final_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  res_gnatho$final,
  "../output/gnathostomata_cne_final_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
