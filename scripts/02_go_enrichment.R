### LIBRARIES ###
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  GenomicRanges,
  GenomicFeatures,
  tidyverse,
  rGREAT,
  GenomeInfoDb,
  org.Dr.eg.db,
  AnnotationDbi,
  tidyplots,
  ggplot2,
  simplifyEnrichment
)

setwd(this.path::here())

source('custom_functions.R')

preproc_dir <- "../output/preprocessed"
great_dir <- "../output/great/"
gviz_dir <- "../output/great/gviz/"
report_dir <- "../output/reports"
dir.create(great_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gviz_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
### 1. LOAD PREPROCESSED INPUTS                                             ###
###############################################################################

actinopteriigy_cne_gr <- readRDS(file.path(
  preproc_dir,
  "actinopteriigy_cne_gr.rds"
))
seqlevels(actinopteriigy_cne_gr, pruning.mode = "coarse") <- str_remove(
  seqlevels(actinopteriigy_cne_gr),
  "\\.[0-9]+$"
)
gnathostomata_cne_gr <- readRDS(file.path(
  preproc_dir,
  "gnathostomata_cne_gr.rds"
))
seqlevels(gnathostomata_cne_gr, pruning.mode = "coarse") <- str_remove(
  seqlevels(gnathostomata_cne_gr),
  "\\.[0-9]+$"
)
drer_sizes <- readRDS(file.path(preproc_dir, "drer_sizes.rds"))
drer_anno <- AnnotationDbi::loadDb(
  file.path(preproc_dir, "drer_anno.sqlite")
)

###############################################################################
### 2. BUILD EXTENDED TSS FROM CUSTOM GFF-DERIVED TxDb                      ###
###############################################################################

gene_gr <- genes(drer_anno, single.strand.genes.only = TRUE)
mcols(gene_gr)$gene_id <- sub("^[a-z]+-", "", mcols(gene_gr)$gene_id)
seqlevels(gene_gr, pruning.mode = "coarse") <- str_remove(
  seqlevels(gene_gr),
  "\\.[0-9]+$"
)

sl <- setNames(drer_sizes$length, drer_sizes$chrom)

gene_gr <- gene_gr[as.character(seqnames(gene_gr)) %in% names(sl)]

gene_gr <- dropSeqlevels(gene_gr, 'NC_002333')
seqlengths(gene_gr) <- sl[seqlevels(gene_gr)]

et <- extendTSS(gene_gr, seqlengths = sl, gene_id_type = "SYMBOL")

###############################################################################
### 3. BUILD GO:BP GENE SETS FROM org.Dr.eg.db                              ###
###############################################################################
# Using GOALL/ONTOLOGYALL propagates annotations up the GO DAG.

go_map <- AnnotationDbi::select(
  org.Dr.eg.db,
  keys = keys(org.Dr.eg.db, "SYMBOL"),
  columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
  keytype = "SYMBOL"
)
go_bp <- go_map |>
  dplyr::filter(ONTOLOGYALL == "BP", !is.na(GOALL)) |>
  dplyr::distinct(GOALL, SYMBOL)

go_bp_list <- split(go_bp$SYMBOL, go_bp$GOALL)
go_bp_list <- go_bp_list[lengths(go_bp_list) >= 10 & lengths(go_bp_list) <= 500]

###############################################################################
### 4. SYMBOL MAPPING RATE                                                  ###
###############################################################################
# Genes in the custom GFF annotation that fail to match an org.Dr.eg.db SYMBOL
# are dropped from the gene sets without notice, so the rate is reported
# explicitly. The second pair of rows gives the fraction of GO gene-set members
# that receive a regulatory domain, which is the effective universe of the
# binomial region test.

symbol_mapping <- symbol_mapping_table(gene_gr, et, go_bp_list, org.Dr.eg.db)

write_tsv(symbol_mapping, file.path(report_dir, "symbol_mapping_rate.tsv"))
print(symbol_mapping)

###############################################################################
### 5. RUN GREAT — WHOLE-GENOME BACKGROUND (DESCRIPTIVE BASELINE)           ###
###############################################################################

res_actino <- great(
  gr = actinopteriigy_cne_gr,
  gene_sets = go_bp_list,
  extended_tss = et
)
res_gnatho <- great(
  gr = gnathostomata_cne_gr,
  gene_sets = go_bp_list,
  extended_tss = et
)

tbl_actino <- getEnrichmentTable(res_actino)
tbl_gnatho <- getEnrichmentTable(res_gnatho)

###############################################################################
### 6. RUN GREAT — UNION-CNE BACKGROUND (PRIMARY LINEAGE CONTRAST)          ###
###############################################################################

cne_background <- IRanges::reduce(
  c(actinopteriigy_cne_gr, gnathostomata_cne_gr),
  ignore.strand = TRUE
)

res_actino_vs_cne <- great(
  gr = actinopteriigy_cne_gr,
  gene_sets = go_bp_list,
  extended_tss = et,
  background = cne_background
)
res_gnatho_vs_cne <- great(
  gr = gnathostomata_cne_gr,
  gene_sets = go_bp_list,
  extended_tss = et,
  background = cne_background
)

tbl_actino_vs_cne <- getEnrichmentTable(res_actino_vs_cne)
tbl_gnatho_vs_cne <- getEnrichmentTable(res_gnatho_vs_cne)

# Background composition, for the Methods.
background_summary <- tibble(
  background = c("whole genome", "CNE union"),
  n_regions = c(NA_integer_, length(cne_background)),
  total_bp = c(
    sum(as.numeric(drer_sizes$length)),
    sum(as.numeric(width(cne_background)))
  )
) %>%
  mutate(
    fraction_of_genome = round(total_bp / sum(as.numeric(drer_sizes$length)), 5)
  )
write_tsv(background_summary, file.path(report_dir, "great_background.tsv"))

###############################################################################
### 7. COMBINED TABLES + DUAL-TEST HEADLINE CALLS                           ###
###############################################################################

combined <- bind_rows(
  combine_great(tbl_actino, "actinopterygii"),
  combine_great(tbl_gnatho, "gnathostomata")
)
combined_vs_cne <- bind_rows(
  combine_great(tbl_actino_vs_cne, "actinopterygii"),
  combine_great(tbl_gnatho_vs_cne, "gnathostomata")
)

combined_strict <- strict_call(combined)
combined_vs_cne_strict <- strict_call(combined_vs_cne)

### --- Significant rGREAT terms export (replaces the old signif_go dump) ---
go_signif <- combined |>
  filter(p.adjust <= 0.05 | p.adjust_hyper <= 0.05)
write.table(
  go_signif,
  "../output/GO_enrichment_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Shared-terms (per-cluster wide) table for the FE-vs-FE comparison
shared_fc <- combined %>%
  select(ID, Description, Cluster, FoldEnrichment, p.adjust, p.adjust_hyper) %>%
  pivot_wider(
    names_from = Cluster,
    values_from = c(FoldEnrichment, p.adjust, p.adjust_hyper),
    names_glue = "{.value}_{Cluster}"
  ) %>%
  drop_na(FoldEnrichment_actinopterygii, FoldEnrichment_gnathostomata) %>%
  mutate(
    delta_FE = FoldEnrichment_actinopterygii - FoldEnrichment_gnathostomata,
    log2ratio_FE = log2(
      FoldEnrichment_actinopterygii /
        FoldEnrichment_gnathostomata
    ),
    signif_both = p.adjust_actinopterygii <= 0.05 &
      p.adjust_gnathostomata <= 0.05 &
      p.adjust_hyper_actinopterygii <= 0.05 &
      p.adjust_hyper_gnathostomata <= 0.05
  ) %>%
  arrange(desc(abs(log2ratio_FE)))

write.table(
  shared_fc,
  file.path(great_dir, "shared_fc_GO_great.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.data.frame(tbl_actino),
  file.path(great_dir, "great_actinopterygii_GObp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.data.frame(tbl_gnatho),
  file.path(great_dir, "great_gnathostomata_GObp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.data.frame(tbl_actino_vs_cne),
  file.path(great_dir, "great_actinopterygii_GObp_vs_cne_bg.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.data.frame(tbl_gnatho_vs_cne),
  file.path(great_dir, "great_gnathostomata_GObp_vs_cne_bg.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
### 8. CNE HOTSPOTS ON GREAT REGULATORY DOMAINS                             ###
###############################################################################
# Nearest-gene assignment with raw counts ranks loci partly by size, because a
# large locus intercepts more elements for the same underlying density. Scoring
# each gene's regulatory domain with the same binomial region model used for the
# GO analysis puts domain length into the null expectation, so obs/exp and its
# FDR are comparable across loci of very different extents. Ohnologue pairs are
# flagged because 3R duplicates appear as two independent domains and their
# summed count is not a single inflated locus.

hotspots <- bind_rows(
  hotspot_domains(actinopteriigy_cne_gr, et, "actinopterygii"),
  hotspot_domains(gnathostomata_cne_gr, et, "gnathostomata")
) |>
  flag_ohnologue_pairs()

write_tsv(hotspots, file.path(great_dir, "cne_hotspots_domain_normalised.tsv"))
saveRDS(hotspots, file.path(great_dir, "cne_hotspots.rds"))

# Raw nearest-gene ranking retained for the direct comparison, so the effect of
# length normalisation on the ranking can be stated quantitatively.
raw_counts <- as_tibble(mcols(actinopteriigy_cne_gr)) %>%
  filter(!is.na(geneId)) %>%
  dplyr::count(geneId, name = "raw_count")

hotspot_rank_shift <- hotspots %>%
  filter(set == "actinopterygii", observed > 0) %>%
  left_join(raw_counts, by = c("gene" = "geneId")) %>%
  mutate(
    rank_domain_normalised = rank(-obs_exp, ties.method = "min"),
    rank_raw_domain_count = rank(-observed, ties.method = "min"),
    rank_change = rank_raw_domain_count - rank_domain_normalised
  ) %>%
  arrange(rank_domain_normalised) %>%
  select(
    gene,
    observed,
    raw_count,
    domain_width,
    expected,
    obs_exp,
    per_100kb,
    fdr,
    rank_raw_domain_count,
    rank_domain_normalised,
    rank_change,
    is_ohnologue
  )

write_tsv(
  hotspot_rank_shift,
  file.path(report_dir, "hotspot_rank_raw_vs_normalised.tsv")
)

# The two loci raised as size-driven, reported per paralogue.
hotspots %>%
  filter(str_detect(gene, "^(ryr|pcdh)")) %>%
  arrange(set, gene) %>%
  write_tsv(file.path(report_dir, "hotspots_ryr_pcdh.tsv"))

###############################################################################
### 9. RGREAT INHERENT VOLCANO PLOTS                                        ###
###############################################################################

actino_col <- "#CC79A7"
gnatho_col <- "#009E73"

pdf(
  file.path(great_dir, "plotVolcano_actinopterygii.pdf"),
  width = 6,
  height = 5
)
print(plotVolcano(res_actino))
dev.off()

pdf(
  file.path(great_dir, "plotVolcano_gnathostomata.pdf"),
  width = 6,
  height = 5
)
print(plotVolcano(res_gnatho))
dev.off()

pdf(
  file.path(great_dir, "plotVolcano_actinopterygii_vs_cne_bg.pdf"),
  width = 6,
  height = 5
)
print(plotVolcano(res_actino_vs_cne))
dev.off()

pdf(
  file.path(great_dir, "plotVolcano_gnathostomata_vs_cne_bg.pdf"),
  width = 6,
  height = 5
)
print(plotVolcano(res_gnatho_vs_cne))
dev.off()

###############################################################################
### 10. TSS DISTANCE DISTRIBUTION (rGREAT REGION-GENE ASSOCIATION PANEL)    ###
###############################################################################
# Distance to TSS is taken from rGREAT's own region-gene association plot
# (second panel), which measures the signed distance from each element to the
# TSS of the genes whose regulatory domain it falls in, in the gene's own
# orientation. This replaces the earlier distance summary and removes any
# dependence on an external distance convention.

pdf(
  file.path(great_dir, "regionGeneAssoc_distTSS_actinopterygii.pdf"),
  width = 5,
  height = 4.5
)
plotRegionGeneAssociations(res_actino, which_plot = 2)
dev.off()

pdf(
  file.path(great_dir, "regionGeneAssoc_distTSS_gnathostomata.pdf"),
  width = 5,
  height = 4.5
)
plotRegionGeneAssociations(res_gnatho, which_plot = 2)
dev.off()

# Numerical version of the same quantity, so the distribution can be described
# in the text rather than only shown.
tss_distances <- list(
  actinopterygii = getRegionGeneAssociations(res_actino),
  gnathostomata = getRegionGeneAssociations(res_gnatho)
)

tss_distance_summary <- imap_dfr(tss_distances, function(rga, nm) {
  d <- unlist(mcols(rga)$dist_to_TSS, use.names = FALSE)
  d <- d[!is.na(d)]
  tibble(
    set = nm,
    n_region_gene_pairs = length(d),
    fraction_upstream = round(mean(d < 0), 3),
    fraction_downstream = round(mean(d > 0), 3),
    median_kb = round(median(d) / 1e3, 1),
    median_abs_kb = round(median(abs(d)) / 1e3, 1)
  )
})

write_tsv(
  tss_distance_summary,
  file.path(report_dir, "tss_distance_summary.tsv")
)
print(tss_distance_summary)

###############################################################################
### 11. TIDYPLOT — FE vs FE SCATTER (SHARED TERMS, log2)                    ###
###############################################################################

fe_df <- shared_fc %>%
  mutate(
    actinopterygii = log2(FoldEnrichment_actinopterygii),
    gnathostomata = log2(FoldEnrichment_gnathostomata)
  ) %>%
  filter(is.finite(actinopterygii), is.finite(gnathostomata))

fe_df %>%
  tidyplot(x = actinopterygii, y = gnathostomata) |>
  add_data_points(alpha = 0.5, color = "black", shape = 21, size = 3) |>
  add(geom_abline(intercept = 0, slope = 1, linetype = "dashed")) |>
  adjust_theme_details(
    panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
    panel.grid = ggplot2::element_line(
      color = "gray80",
      linewidth = 0.5,
      linetype = 2
    )
  ) |>
  add_data_labels_repel(
    data = max_rows(actinopterygii, n = 8),
    label = Description,
    color = "black",
    background = TRUE,
    min.segment.length = 0
  ) |>
  adjust_font(fontsize = 12) |>
  adjust_size(width = 4, height = 3, unit = "in") |>
  adjust_x_axis_title("log2 FE (actinopterygii)") |>
  adjust_y_axis_title("log2 FE (gnathostomata)") |>
  save_plot(
    file.path(great_dir, "shared_GO_FE_vs_FE.pdf"),
    width = 6,
    height = 5,
    units = "in"
  )

###############################################################################
### 12. TIDYPLOT — SHARED-TERM VOLCANO (-log10 p.adj vs log2 FE)            ###
###############################################################################
# "Shared" here means significant in both sets, not merely tested in both.

shared_ids <- combined %>%
  filter(p.adjust <= 0.05) %>%
  dplyr::count(ID) %>%
  filter(n > 1) %>%
  pull(ID)

combined %>%
  filter(ID %in% shared_ids) %>%
  mutate(
    FoldEnrichment = log2(FoldEnrichment),
    p.adjust = -log10(p.adjust)
  ) %>%
  tidyplot(x = FoldEnrichment, y = p.adjust, color = Cluster) |>
  add_data_points(alpha = 0.5, color = "black", shape = 21, size = 3) |>
  adjust_theme_details(
    panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
    panel.grid = ggplot2::element_line(
      color = "gray80",
      linewidth = 0.5,
      linetype = 2
    )
  ) |>
  adjust_font(fontsize = 12) |>
  adjust_size(width = 4, height = 3, unit = "in") |>
  adjust_x_axis_title("log2(FoldEnrichment)") |>
  adjust_y_axis_title("-log10(adjusted p-value)") |>
  adjust_colors(new_colors = c(actino_col, gnatho_col)) |>
  save_plot(
    file.path(great_dir, "shared_GO_volcano_compare.pdf"),
    width = 6,
    height = 5,
    units = "in"
  )

###############################################################################
### 13. simplifyGO — SEPARATE TERM CLUSTERING PER CNE CLASS                 ###
###############################################################################

go_clusters_actino <- run_simplify(tbl_actino, "actinopterygii")
go_clusters_gnatho <- run_simplify(tbl_gnatho, "gnathostomata")

###############################################################################
### 14. REGION-GENE ASSOCIATIONS FOR GO:0007224 (SMOOTHENED SIGNALLING)     ###
###     + PER-TERM ASSOCIATIONS FOR DOWNSTREAM SCRIPT                       ###
###############################################################################

pdf(
  file.path(great_dir, "plotRegionGeneAssoc_GO0007224_actinopterygii.pdf"),
  width = 9,
  height = 4
)
plotRegionGeneAssociations(res_actino, term_id = "GO:0007224")
dev.off()

pdf(
  file.path(great_dir, "plotRegionGeneAssoc_GO0007224_gnathostomata.pdf"),
  width = 9,
  height = 4
)
plotRegionGeneAssociations(res_gnatho, term_id = "GO:0007224")
dev.off()

# Pull regions and gene names for Gviz later
rga_actino_shh <- getRegionGeneAssociations(res_actino, term_id = "GO:0007224")
rga_gnatho_shh <- getRegionGeneAssociations(res_gnatho, term_id = "GO:0007224")

shh_genes_actino <- unique(unlist(mcols(rga_actino_shh)$annotated_genes))
shh_genes_gnatho <- unique(unlist(mcols(rga_gnatho_shh)$annotated_genes))

saveRDS(
  list(
    term_id = "GO:0007224",
    description = "smoothened signaling pathway",
    actino_regions = rga_actino_shh,
    gnatho_regions = rga_gnatho_shh,
    actino_genes = shh_genes_actino,
    gnatho_genes = shh_genes_gnatho
  ),
  file.path(great_dir, "GO0007224_smoothened_regions_genes.rds")
)

# Also pre-compute region-gene associations for a list of biologically
# important terms so 03_post_process.R can read them without touching rGREAT.
key_descriptions <- c(
  "smoothened signaling pathway",
  "fin morphogenesis",
  "fin development",
  "pectoral fin development",
  "pectoral fin morphogenesis",
  "skeletal system development",
  "muscle structure development"
)

# Map descriptions to actual IDs present in the rGREAT results
desc_to_id <- combined %>%
  filter(Description %in% key_descriptions) %>%
  distinct(ID, Description)

extract_associations <- function(res, ids) {
  out <- list()
  for (tid in ids) {
    out[[tid]] <- tryCatch(
      getRegionGeneAssociations(res, term_id = tid),
      error = function(e) NULL
    )
  }
  out
}

term_associations <- list(
  descriptions = desc_to_id,
  actino = extract_associations(res_actino, desc_to_id$ID),
  gnatho = extract_associations(res_gnatho, desc_to_id$ID)
)
saveRDS(
  term_associations,
  file.path(great_dir, "term_region_gene_associations.rds")
)

###############################################################################
### 15. PERSIST CORE RESULTS FOR DOWNSTREAM SCRIPT                          ###
###############################################################################

saveRDS(combined, file.path(great_dir, "combined_GObp.rds"))
saveRDS(combined_vs_cne, file.path(great_dir, "combined_GObp_vs_cne_bg.rds"))
saveRDS(shared_fc, file.path(great_dir, "shared_fc.rds"))
saveRDS(et, file.path(preproc_dir, "extended_tss.rds"))

# Strict (dual-test) calls
write.table(
  combined_strict,
  file.path(great_dir, "great_strict_dualtest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_vs_cne_strict,
  file.path(great_dir, "great_strict_dualtest_vs_cne_bg.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fin_terms <- c("fin morphogenesis", "fin development")
fin_skeletal_terms <- c(
  "fin morphogenesis",
  "fin development",
  "skeletal system development",
  "muscle structure development"
)

filtered_genes_for_downstream <- list(
  fin = extract_genes_for_terms(res_actino, tbl_actino, fin_terms),
  fin_skeletal = extract_genes_for_terms(
    res_actino,
    tbl_actino,
    fin_skeletal_terms
  ),
  all_signif_actino = unique(unlist(lapply(
    tbl_actino$id[tbl_actino$p_adjust <= 0.05],
    function(tid) {
      unique(unlist(
        getRegionGeneAssociations(res_actino, term_id = tid)$annotated_genes
      ))
    }
  )))
)

saveRDS(
  filtered_genes_for_downstream,
  file.path(great_dir, "filtered_genes_for_downstream.rds")
)

### --- GO:0007224 (smoothened signaling) → Gviz-ready exports --------------
genes_smo <- unique(c(shh_genes_actino, shh_genes_gnatho))
gene_gr_smo <- gene_gr[mcols(gene_gr)$gene_id %in% genes_smo]

saveRDS(gene_gr_smo, file.path(gviz_dir, "GO0007224_genes.rds"))
saveRDS(rga_actino_shh, file.path(gviz_dir, "GO0007224_cnes_actino.rds"))
saveRDS(rga_gnatho_shh, file.path(gviz_dir, "GO0007224_cnes_gnatho.rds"))
