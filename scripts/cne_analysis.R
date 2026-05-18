### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  gprofiler2,
  rGREAT,
  regioneR,
  patchwork,
  RColorBrewer,
  Biostrings,
  ChIPseeker,
  GenomicFeatures,
  zFPKM,
  SummarizedExperiment,
  readxl,
  UpSetR,
  ComplexHeatmap,
  GenometriCorr
)
### SET WORKING DIRECTORY ###
setwd(this.path::here())

### READ IN DATA ###
# Read in CNE bed files
bed_colnames <- c(
  'chromosome',
  'start',
  'end',
  'cne_name',
  'phastcons_score',
  'strand'
)
actinopteriigy_cne <- read_tsv(
  "../input/actinopterygii_specific_drer.bed",
  col_names = F
)
gnathostomata_cne <- read_tsv(
  "../input/gnathostomata_conserved_drer.bed",
  col_names = F
)

colnames(actinopteriigy_cne) <- bed_colnames
colnames(gnathostomata_cne) <- bed_colnames

drer_sizes <- read_tsv(
  "../ancilliary_files/drer_chrom_info.txt",
  col_names = FALSE
)
colnames(drer_sizes) <- c("chrom", "length")

drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)

# Read in salmon TPM expression data
salmon_tpm <- read_tsv("../input/salmon.merged.gene_tpm.tsv")

### GENE ACTIVITY CLASSIFICATION USING zFPKM ###
# Extract TPM columns (replicates) as a matrix
tpm_matrix <- salmon_tpm %>%
  dplyr::select(starts_with("SRR")) %>%
  as.matrix()

# Set row names to gene IDs
rownames(tpm_matrix) <- salmon_tpm$gene_id

# Create SummarizedExperiment object with TPM data
se <- SummarizedExperiment(
  assays = SimpleList(fpkm = tpm_matrix),
  rowData = DataFrame(
    gene_id = salmon_tpm$gene_id,
    gene_name = salmon_tpm$gene_name
  )
)

# Compute zFPKM using the multivariate Gaussian mixture modeling approach
assay(se, "zfpkm") <- zFPKM(se)

# Calculate mean TPM across replicates
salmon_tpm <- salmon_tpm %>%
  mutate(
    mean_tpm = rowMeans(dplyr::select(., starts_with("SRR")), na.rm = TRUE)
  )

# Extract zFPKM values and create gene activity lookup table
# Genes with zFPKM > -3 are considered active
gene_activity <- salmon_tpm %>%
  mutate(
    z_fpkm = assay(se, "zfpkm")[, 1], # Use first zFPKM column (averaged across samples)
    is_active = z_fpkm > -3 # Standard zFPKM threshold for active genes
  ) %>%
  dplyr::select(gene_id, gene_name, mean_tpm, z_fpkm, is_active) %>%
  mutate(gene_id = str_remove_all(gene_id, "^[a-z]+-"))

saveRDS(gene_activity, '../output/gene_activity.RDS')

### EXPLORATORY DATA ANALYSIS (EDA) ###
# Convert CNEs and mask data to GRanges objects
actinopteriigy_cne_gr <- GRanges(
  seqnames = sub("\\.[0-9]+$", "", actinopteriigy_cne$cne_name),
  ranges = IRanges(
    start = actinopteriigy_cne$start,
    end = actinopteriigy_cne$end,
  ),
  phastcons = actinopteriigy_cne$phastcons_score
)

gnathostomata_cne_gr <- GRanges(
  seqnames = sub("\\.[0-9]+$", "", gnathostomata_cne$cne_name),
  ranges = IRanges(
    start = gnathostomata_cne$start,
    end = gnathostomata_cne$end,
  ),
  phastcons = gnathostomata_cne$phastcons_score
)
crossmapping <- data.frame(
  nonstrip = seqlevels(actinopteriigy_cne_gr),
  stripped = sub('\\.[0-9]+$', '', seqlevels(actinopteriigy_cne_gr))
)
drer_sizes$chrom <- crossmapping$nonstrip[match(
  drer_sizes$chrom,
  crossmapping$stripped
)]

# Strip version suffixes from both sides and match
drer_sizes$chrom_base <- sub("\\..*", "", drer_sizes$chrom)
sl_base <- sub("\\..*", "", names(seqlengths(actinopteriigy_cne_gr)))

sl <- seqlengths(actinopteriigy_cne_gr)
names(sl) <- sl_base # temporarily strip versions

sl[drer_sizes$chrom_base] <- drer_sizes$length

names(sl) <- seqlevels(actinopteriigy_cne_gr) # restore original names
seqlengths(actinopteriigy_cne_gr) <- sl[which(sl > 20000)]

### IMPORT SHEET 7 FROM YUESONG ET AL 2025 AND LIFTOVER ###
# Read sheet 7 from the Excel file (pseudo-bed format: chrom, start, end, and additional columns)
sheet7_data <- read_excel(
  "../input/YueSong-et-al_2025.xlsx",
  sheet = 7,
  col_names = T,
  skip = 1
)

alias_url <- "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/002/035/GCF_000002035.6/GCF_000002035.6.chromAlias.txt"

alias <- read_tsv(
  alias_url,
  comment = "#",
  show_col_types = FALSE,
  col_names = F
)

# Build RefSeq -> UCSC map (NC_.... -> chrN)
refseq2ucsc <- setNames(alias$X5, alias$X1)

sheet7_data$Chromosome <- alias$X5[match(sheet7_data$Chromosome, alias$X1)]

# Convert to GRanges object (assuming first 3 columns are chrom, start, end)
sheet7_gr <- GRanges(
  seqnames = sheet7_data[[1]],
  ranges = IRanges(start = sheet7_data[[2]], end = sheet7_data[[3]]),
  mcols = sheet7_data[, 4:ncol(sheet7_data)]
)

liftover <- read_tsv(
  '../input/ucsc_GRCz11-GRCz12_liftover_actinopteriigy_cne.bed',
  col_names = F
)
# IMPORTANT: attach X4 as metadata (cne_name)
liftover_gr <- GRanges(
  seqnames = liftover$X1,
  ranges = IRanges(start = liftover$X2 + 1, end = liftover$X3), # BED start is 0-based
  strand = "*",
  cne_name = liftover$X4
)

hits <- findOverlaps(liftover_gr, sheet7_gr, ignore.strand = TRUE)

overlapping_cne_names <- unique(mcols(liftover_gr)$cne_name[queryHits(hits)])

actinopteriigy_cne_ov <- actinopteriigy_cne %>%
  filter(cne_name %in% overlapping_cne_names)

# Annotate peaks and create visualizations
peak_anno_list <- lapply(
  list(actinopteriigy_cne_gr, gnathostomata_cne_gr),
  annotatePeak,
  overlap = 'all',
  TxDb = drer_anno,
  tssRegion = c(-3000, 3000),
  genomicAnnotationPriority = c(
    "Intergenic",
    "Downstream",
    "Promoter",
    "5UTR",
    "3UTR",
    "Intron",
    "Exon"
  )
)
names(peak_anno_list) <- c("actinopteriigy_CNE", "gnathostomata_CNE")

# Add gene expression information to annotations
# For each annotation, find the nearest gene and add its expression status
peak_anno_list <- lapply(peak_anno_list, function(anno_obj) {
  anno_df <- as_tibble(anno_obj@anno) %>%
    # Join with gene activity info by nearest gene (geneId)
    left_join(gene_activity, by = c("geneId" = "gene_id"))

  # Update the annotation object
  anno_obj@anno <- as(anno_df, "GRanges")
  anno_obj
})

#####################
### GO ernichment ###
#####################

# Load libraries
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(enrichplot)

gene_list_actinopteriigy <- peak_anno_list$actinopteriigy_CNE@anno$geneId[
  !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, "Exon")) &
    peak_anno_list$actinopteriigy_CNE@anno$is_active == T
]
gene_list_actinopteriigy <- unique(gene_list_actinopteriigy)

gene_list_verte <- peak_anno_list$gnathostomata_CNE@anno$geneId[
  !(str_detect(peak_anno_list$gnathostomata_CNE@anno$annotation, "Exon")) &
    peak_anno_list$gnathostomata_CNE@anno$is_active == T
]
gene_list_verte <- unique(gene_list_verte)

gene_list <- list(
  actinopteriigy = gene_list_actinopteriigy,
  gnathostomata = gene_list_verte
)

# Using gene symbols directly
ego_bp <- compareCluster(
  gene = gene_list,
  OrgDb = org.Dr.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

# Remove redundant GO terms
ego_bp <- simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(ego_bp, showCategory = 10) +
  theme_bw(base_size = 12)
ggsave('../output/go_enrich_dotplot.pdf', width = 6, height = 8, units = 'in')

# View results
head(ego_bp@compareClusterResult)

signif_go <- ego_bp@compareClusterResult[
  ego_bp@compareClusterResult$qvalue <= 0.05,
]
signif_go_filtered_fin <- signif_go[
  signif_go$Description %in%
    c(
      "fin morphogenesis",
      "fin development",
      "skeletal system development",
      "muscle structure development"
    ),
]

filtered_genes <- unique(unlist(strsplit(signif_go_filtered_fin$geneID, "/")))

actinopteriigy_filtered <- as.data.frame(peak_anno_list$actinopteriigy_CNE@anno[
  (peak_anno_list$actinopteriigy_CNE@anno$geneId %in% filtered_genes) &
    !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
])
vertebrate_filtered <- as.data.frame(peak_anno_list$gnathostomata_CNE@anno[
  (peak_anno_list$gnathostomata_CNE@anno$geneId %in% filtered_genes) &
    !(str_detect(peak_anno_list$gnathostomata_CNE@anno$annotation, 'Exon')),
])


df <- ego_bp@compareClusterResult

tele <- df %>%
  filter(Cluster == "actinopteriigy") %>%
  select(ID, Description, FE_actinopteriigy = FoldEnrichment)

vert <- df %>%
  filter(Cluster == "Vertebrate") %>%
  select(ID, Description, FE_Vertebrate = FoldEnrichment)

shared_fc <- inner_join(tele, vert, by = c("ID", "Description")) %>%
  mutate(
    delta_FE = FE_actinopteriigy - FE_Vertebrate,
    log2ratio_FE = log2(FE_actinopteriigy / FE_Vertebrate)
  ) %>%
  arrange(desc(abs(delta_FE)))

shared_fc

actinopteriigy_line <- "#1F78B4"
verte_line <- "skyblue"

library(tidyplots)

ego_bp@compareClusterResult %>%
  group_by(Description) %>%
  filter(n() > 1) %>%
  mutate(FoldEnrichment = log2(FoldEnrichment), p.adjust = -log10(p.adjust)) %>%
  tidyplot(x = FoldEnrichment, y = p.adjust, color = Cluster) |>
  add_data_points(alpha = 0.5, color = 'black', shape = 21, size = 3) |>
  adjust_theme_details(
    panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
    panel.grid = ggplot2::element_line(
      color = 'gray80',
      linewidth = 0.5,
      linetype = 2
    )
  ) |>
  adjust_font(fontsize = 12) |>
  adjust_size(width = 4, heigh = 3, unit = 'in') |>
  adjust_x_axis_title('log2(FoldEnrichment)') |>
  adjust_y_axis_title('-log10(adjusted p-value)') |>
  adjust_x_axis(limits = c(-2, 2)) |>
  adjust_colors(new_colors = c(actinopteriigy_line, verte_line)) |>
  save_plot(
    '../output/shared_GO_foldenrichment_difference.pdf',
    width = 6,
    height = 5,
    units = 'in'
  )

lf2lf <- ego_bp@compareClusterResult %>%
  mutate(FE = log2(FoldEnrichment)) %>%
  select(Description, Cluster, FE) %>%
  pivot_wider(names_from = Cluster, values_from = FE) %>%
  drop_na(actinopteriigy, gnathostomata) %>%
  mutate(delta = actinopteriigy - gnathostomata)

lf2lf %>%
  tidyplot(x = actinopteriigy, y = gnathostomata) |>
  add_data_points(alpha = 0.5, color = 'black', shape = 21, size = 3) |>
  add(geom_abline(intercept = 0, slope = 1, linetype = 'dashed')) |>
  adjust_theme_details(
    panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
    panel.grid = ggplot2::element_line(
      color = 'gray80',
      linewidth = 0.5,
      linetype = 2
    )
  ) |>
  add_data_labels_repel(
    data = max_rows(actinopteriigy, n = 5),
    label = Description,
    color = "black",
    background = TRUE,
    min.segment.length = 0
  ) |>
  adjust_y_axis(limits = c(0, 2)) |>
  adjust_x_axis(limits = c(0, 2)) |>
  adjust_font(fontsize = 12) |>
  adjust_size(width = 4, heigh = 3, unit = 'in') |>
  adjust_x_axis_title('Fold Enrichment in Actinopteriigy') |>
  adjust_y_axis_title('Fold Enrichment in Gnathostomata') |>
  save_plot(
    '../output/shared_GO_foldenrichment_comparisons.pdf',
    width = 6,
    height = 5,
    units = 'in'
  )

##############
### EXPORT ###
##############

export <- peak_anno_list$gnathostomata_CNE@anno[
  !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
]
write.table(
  export,
  '../output/gnathostomata_specific_cne.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)


export <- peak_anno_list$actinopteriigy_CNE@anno[
  !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
]
signif_go_filtered_skeletal_fin <- signif_go[
  signif_go$Description %in%
    c(
      "fin morphogenesis",
      "fin development",
      "skeletal system development",
      "muscle structure development"
    ),
]
signif_go_filtered_fin <- signif_go[
  signif_go$Description %in% c("fin morphogenesis", "fin development"),
]

filtered_genes_skeletal_fin <- unique(unlist(strsplit(
  signif_go_filtered_skeletal_fin$geneID,
  "/"
)))
filtered_genes_fin <- unique(unlist(strsplit(
  signif_go_filtered_fin$geneID,
  "/"
)))

actinopteriigy_filtered_skeletal_fin <- as.data.frame(peak_anno_list$actinopteriigy_CNE@anno[
  (peak_anno_list$actinopteriigy_CNE@anno$geneId %in%
    filtered_genes_skeletal_fin) &
    !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
])
actinopteriigy_filtered_fin <- as.data.frame(peak_anno_list$actinopteriigy_CNE@anno[
  (peak_anno_list$actinopteriigy_CNE@anno$geneId %in% filtered_genes_fin) &
    !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
])

write.table(
  actinopteriigy_filtered_fin,
  '../output/actinopteriigy_GO-fin_specific_cne.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
write.table(
  actinopteriigy_filtered_skeletal_fin,
  '../output/actinopteriigy_GO-fin-skeletal_specific_cne.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)


go_signif <- ego_bp@compareClusterResult[
  ego_bp@compareClusterResult$qvalue <= 0.05,
]

write.table(
  go_signif,
  '../output/GO_enrichment_results.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)

### ACTIVE GENE FILTERING ###

# Export CNEs near active genes (all annotations, non-exonic)
export_active <- as.data.frame(peak_anno_list$actinopteriigy_CNE@anno[
  (peak_anno_list$actinopteriigy_CNE@anno$nearest_gene_active == TRUE) &
    !(str_detect(peak_anno_list$actinopteriigy_CNE@anno$annotation, 'Exon')),
])

write.table(
  export_active,
  '../output/actinopteriigy_cne_near_active_genes.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
# Export GO-filtered CNEs that are near active genes (fin-specific)
actinopteriigy_filtered_fin_active <- actinopteriigy_filtered_fin %>%
  filter(is_active == TRUE)

write.table(
  actinopteriigy_filtered_fin_active,
  '../output/actinopteriigy_GO-fin_specific_cne_active_genes.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
# Export GO-filtered CNEs that are near active genes (fin+skeletal-specific)
actinopteriigy_filtered_skeletal_fin_active <- actinopteriigy_filtered_skeletal_fin %>%
  filter(is_active == TRUE)

write.table(
  actinopteriigy_filtered_skeletal_fin_active,
  '../output/actinopteriigy_GO-fin-skeletal_specific_cne_active_genes.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
############################
### ATACseq Peak Overlap ###
############################

# Read in ATACseq consensus peaks BED file
atac_peaks <- read_tsv(
  "../input/consensus_peaks.mLb.clN.bed",
  col_names = FALSE
)
colnames(atac_peaks) <- c(
  'chromosome',
  'start',
  'end',
  'peak_name',
  'score',
  'strand'
)

# Convert ATACseq peaks to GRanges object
atac_peaks_gr <- GRanges(
  seqnames = atac_peaks$chromosome,
  ranges = IRanges(start = atac_peaks$start, end = atac_peaks$end),
  peak_name = atac_peaks$peak_name,
  score = atac_peaks$score
)

# Find overlaps between actinopteriigy CNEs and ATACseq peaks
cne_atac_overlaps <- findOverlaps(actinopteriigy_cne_gr, atac_peaks_gr)
verte_atac_overlaps <- findOverlaps(gnathostomata_cne_gr, atac_peaks_gr)

# Extract overlapping CNEs and peaks
overlapping_cnes <- actinopteriigy_cne_gr[queryHits(cne_atac_overlaps)]
overlapping_peaks <- atac_peaks_gr[subjectHits(cne_atac_overlaps)]
overlapping_verte <- gnathostomata_cne_gr[queryHits(verte_atac_overlaps)]

# Export overlapping CNEs with annotation (non-exonic) and gene expression info using proper overlap logic
anno_gr <- peak_anno_list$actinopteriigy_CNE@anno
anno_atac_overlaps <- findOverlaps(anno_gr, atac_peaks_gr)

# Subset annotated CNEs that overlap ATACseq peaks and are non-exonic
cne_atac_annotated <- anno_gr[queryHits(anno_atac_overlaps)]
cne_atac_annotated <- cne_atac_annotated[
  !(str_detect(cne_atac_annotated$annotation, 'Exon'))
]

cne_atac_annotated_df <- as.data.frame(cne_atac_annotated)
write.table(
  cne_atac_annotated_df,
  '../output/actinopteriigy_cne_atac_overlaps_annotated.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)

# Export ATACseq-overlapping CNEs from ACTIVE genes only
cne_atac_annotated_active <- cne_atac_annotated_df %>%
  filter(is_active == TRUE)

write.table(
  cne_atac_annotated_active,
  '../output/actinopteriigy_cne_atac_overlaps_active_genes.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
# Create visualization of overlap distribution
overlap_by_annotation <- cne_atac_annotated_df %>%
  as_tibble() %>%
  mutate(annotation = str_remove(annotation, " \\(.*$")) %>%
  dplyr::count(annotation) %>%
  mutate(perc = n / sum(n) * 100) %>%
  arrange(desc(n))

overlap_plot <- ggplot(
  overlap_by_annotation,
  aes(x = reorder(annotation, -n), y = n, fill = annotation)
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
    title = "Distribution of ATACseq Peak-overlapping CNEs by Annotation"
  )

ggsave(
  "../output/cne_atac_overlaps_by_annotation.pdf",
  overlap_plot,
  width = 8,
  height = 6,
  units = 'in',
  dpi = 300
)

# Filter ATACseq-overlapping CNEs to fin-specific genes only
cne_atac_fin_specific <- cne_atac_annotated_df %>%
  filter(geneId %in% filtered_genes_fin)

# Export fin-specific ATACseq overlap results
write.table(
  cne_atac_fin_specific,
  '../output/actinopteriigy_cne_atac_overlaps_fin_specific.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)

# Filter to fin-specific AND active genes
cne_atac_fin_specific_active <- cne_atac_fin_specific %>%
  filter(is_active == TRUE)

write.table(
  cne_atac_fin_specific_active,
  '../output/actinopteriigy_cne_atac_overlaps_fin_specific_active_genes.tsv',
  sep = '\t',
  quote = F,
  col.names = T,
  row.names = F
)
# Calculate fin-specific overlap statistics
n_atac_fin_specific <- nrow(cne_atac_fin_specific)
pct_atac_fin_specific <- (n_atac_fin_specific / nrow(cne_atac_annotated_df)) *
  100


#############################
### actinopteriigy-unique + ATAC ###
#############################

# Select CNEs that are unique to actinopteriigys (no overlap with gnathostomata CNEs)
ov_tel_verte <- findOverlaps(actinopteriigy_cne_gr, gnathostomata_cne_gr)
actinopteriigy_only_idx <- setdiff(
  seq_along(actinopteriigy_cne_gr),
  unique(queryHits(ov_tel_verte))
)
actinopteriigy_only_gr <- actinopteriigy_cne_gr[actinopteriigy_only_idx]

# Find actinopteriigy-unique CNEs that overlap ATAC peaks
ov_tel_atac <- findOverlaps(actinopteriigy_only_gr, atac_peaks_gr)
tel_only_with_atac <- actinopteriigy_only_gr[unique(queryHits(ov_tel_atac))]

# Map annotated entries to these ranges
anno_gr <- peak_anno_list$actinopteriigy_CNE@anno
ov_anno_telonly <- findOverlaps(anno_gr, tel_only_with_atac)
final_anno_tel <- anno_gr[unique(queryHits(ov_anno_telonly))]
final_df_tel <- as.data.frame(final_anno_tel)

out_path <- '../output/actinopteriigy_unique_cne_with_atac_annotated.tsv'
write.table(
  final_df_tel,
  out_path,
  sep = '	',
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

# Filter to only those near active genes
final_df_tel_active <- final_df_tel %>%
  filter(is_active == TRUE)

out_path_active <- '../output/actinopteriigy_unique_cne_with_atac_annotated_active_genes.tsv'
write.table(
  final_df_tel_active,
  out_path_active,
  sep = '	',
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)


#############################
### EXPORT: CNE GRanges   ###
#############################

# These objects are exported as .rds for downstream visualization/statistical overlap testing (e.g. GenometriCorr)
export_rds_dir <- "../output/granges_rds"
dir.create(export_rds_dir, showWarnings = FALSE, recursive = TRUE)

# Core inputs
saveRDS(
  actinopteriigy_cne_gr,
  file = file.path(export_rds_dir, "actinopteriigy_cne_gr.rds")
)
saveRDS(
  gnathostomata_cne_gr,
  file = file.path(export_rds_dir, "gnathostomata_cne_gr.rds")
)
saveRDS(atac_peaks_gr, file = file.path(export_rds_dir, "atac_peaks_gr.rds"))

# CNE ↔ ATAC overlaps
saveRDS(
  overlapping_cnes,
  file = file.path(export_rds_dir, "actinopteriigy_cne_overlapping_atac_gr.rds")
)
saveRDS(
  overlapping_peaks,
  file = file.path(
    export_rds_dir,
    "atac_peaks_overlapping_actinopteriigy_cne_gr.rds"
  )
)
saveRDS(
  cne_atac_annotated,
  file = file.path(
    export_rds_dir,
    "actinopteriigy_cne_atac_overlaps_annotated_gr.rds"
  )
)

# Vertebrate shared CNEs
saveRDS(
  overlapping_verte,
  file = file.path(export_rds_dir, "gnathostomata_cne_overlapping_atac_gr.rds")
)

# CNE annotations
saveRDS(
  peak_anno_list,
  file = file.path(export_rds_dir, 'peak_annotations.rds')
)

# --- CNE subsets for downstream overlap/statistics ---

# 1) actinopteriigy CNEs with (any) active flanking gene nearby
anno_tel_gr <- peak_anno_list$actinopteriigy_CNE@anno
anno_ver_gr <- peak_anno_list$gnathostomata_CNE@anno

# Prefer flank activity flag if present; otherwise fall back to nearest gene activity flag ("is_active")
if ("is_active" %in% colnames(mcols(anno_tel_gr))) {
  cne_with_active_flanks_gr <- anno_tel_gr[which(
    !is.na(anno_tel_gr$is_active) & anno_tel_gr$is_active
  )]
  cne_with_active_flanks_ver_gr <- anno_ver_gr[which(
    !is.na(anno_ver_gr$is_active) & anno_ver_gr$is_active
  )]
} else if ("is_active" %in% colnames(mcols(anno_tel_gr))) {
  cne_with_active_flanks_gr <- anno_tel_gr[which(
    !is.na(anno_tel_gr$is_active) & anno_tel_gr$is_active
  )]
  cne_with_active_flanks_ver_gr <- anno_ver_gr[which(
    !is.na(anno_ver_gr$is_active) & anno_ver_gr$is_active
  )]
} else {
  cne_with_active_flanks_gr <- anno_tel_gr[0]
  cne_with_active_flanks_ver_gr <- anno_ver_gr[0]
}

# 2) actinopteriigy CNEs with active nearby genes AND those genes have GO annotation (fin/skeletal)
# (uses the precomputed GO-filtered gene sets in this script)
if (exists("filtered_genes")) {
  idx_go <- which(anno_tel_gr$geneId %in% filtered_genes)
  idx_go_ver <- which(anno_ver_gr$geneId %in% filtered_genes)
} else {
  idx_go <- integer(0)
  idx_go_ver <- integer(0)
}

# nearest-gene activity flag is "is_active" from gene_activity join; if absent, use is_active
if ("is_active" %in% colnames(mcols(anno_tel_gr))) {
  idx_active_nearest <- which(
    !is.na(anno_tel_gr$is_active) & anno_tel_gr$is_active
  )
  idx_active_nearest_ver <- which(
    !is.na(anno_ver_gr$is_active) & anno_ver_gr$is_active
  )
} else if ("is_active" %in% colnames(mcols(anno_tel_gr))) {
  idx_active_nearest <- which(
    !is.na(anno_tel_gr$is_active) & anno_tel_gr$is_active
  )
  idx_active_nearest_ver <- which(
    !is.na(anno_ver_gr$is_active) & anno_ver_gr$is_active
  )
} else {
  idx_active_nearest <- integer(0)
  idx_active_nearest_ver <- integer(0)
}

cne_active_go_fin_skeletal_gr <- anno_tel_gr[intersect(
  idx_go,
  idx_active_nearest
)]
cne_active_go_fin_skeletal_ver_gr <- anno_ver_gr[intersect(
  idx_go_ver,
  idx_active_nearest
)]

# 3) actinopteriigy CNEs with active nearby genes AND ATAC overlap
hits_active_atac <- findOverlaps(
  cne_with_active_flanks_gr,
  atac_peaks_gr,
  ignore.strand = TRUE
)
cne_active_with_atac_gr <- cne_with_active_flanks_gr[unique(queryHits(
  hits_active_atac
))]

hits_active_atac <- findOverlaps(
  cne_with_active_flanks_ver_gr,
  atac_peaks_gr,
  ignore.strand = TRUE
)
cne_active_with_atac_ver_gr <- cne_with_active_flanks_ver_gr[unique(queryHits(
  hits_active_atac
))]


# 4) Genome annotation as GRanges (full GFF import + genes-only convenience export)
genome_annotation_gr <- rtracklayer::import("../ancilliary_files/drer.gff")
genome_genes_gr <- genes(drer_anno)

# Export new GRanges
saveRDS(
  cne_with_active_flanks_gr,
  file = file.path(
    export_rds_dir,
    "actinopteriigy_cne_with_active_flanking_genes_gr.rds"
  )
)
saveRDS(
  cne_active_go_fin_skeletal_gr,
  file = file.path(
    export_rds_dir,
    "actinopteriigy_cne_active_go_fin_skeletal_gr.rds"
  )
)
saveRDS(
  cne_active_with_atac_gr,
  file = file.path(export_rds_dir, "actinopteriigy_cne_active_with_atac_gr.rds")
)
saveRDS(
  cne_with_active_flanks_ver_gr,
  file = file.path(
    export_rds_dir,
    "gnathostomata_cne_with_active_flanking_genes_gr.rds"
  )
)
saveRDS(
  cne_active_go_fin_skeletal_ver_gr,
  file = file.path(
    export_rds_dir,
    "gnathostomata_cne_active_go_fin_skeletal_gr.rds"
  )
)
saveRDS(
  cne_active_with_atac_ver_gr,
  file = file.path(export_rds_dir, "gnathostomata_cne_active_with_atac_gr.rds")
)
saveRDS(
  genome_annotation_gr,
  file = file.path(export_rds_dir, "genome_annotation_gr.rds")
)
saveRDS(
  genome_genes_gr,
  file = file.path(export_rds_dir, "genome_genes_gr.rds")
)

###################################
### UpSet Plot - CNE Categories ###
###################################
library(GenomicRanges)
library(dplyr)
library(stringr)
library(ComplexHeatmap)

source('custom_functions.R')

# ---------- Universe ----------
anno_gr <- unique(peak_anno_list$actinopteriigy_CNE@anno)

# optional filter
anno_gr2 <- anno_gr[!str_detect(anno_gr$annotation, "Exon")]

U <- unique(gr_id(anno_gr2)) # universe IDs

# ---------- Set 1: actinopteriigy CNEs ----------
S_actinopteriigy <- U

# ---------- Set 2: ATAC peaks overlap ----------
hits_atac <- findOverlaps(anno_gr2, atac_peaks_gr, ignore.strand = TRUE)
S_atac <- unique(gr_id(anno_gr2[unique(queryHits(hits_atac))]))

# ---------- Set 3: active genes nearby ----------
# (uses your metadata flag, but produces an ID set)
subset_gr_true <- function(gr, col) {
  stopifnot(col %in% colnames(mcols(gr)))
  idx <- which(!is.na(mcols(gr)[[col]]) & mcols(gr)[[col]])
  gr[idx]
}

anno_gr_active <- subset_gr_true(anno_gr2, "is_active")
S_active <- unique(gr_id(anno_gr_active))

# ---------- Set 4: fin developmental genes nearby ----------
S_fin <- unique(gr_id(anno_gr2[anno_gr2$geneId %in% filtered_genes_fin]))

# ---------- Set 5: YueSong Sheet7 overlap (actinopteriigy_cne_ov) ----------
S_sheet7 <- paste0(
  as.character(actinopteriigy_cne_ov$chromosome),
  ":",
  as.character(actinopteriigy_cne_ov$start),
  "-",
  as.character(actinopteriigy_cne_ov$end)
)

# ---------- Set 5: Chan el al.  overlap (actinopteriigy_cne_ov) ----------
ep <- read_tsv('../input/Chan_et_al_EP-loops.tsv')
en <- read_tsv('../ancilliary_files/enhancer.grcz12.bed', col_names = F)
pr <- read_tsv('../ancilliary_files/promoter.grcz12.bed', col_names = F)

enh_gr <- GRanges(
  seqnames = en$X1,
  ranges = IRanges(start = en$X2 + 1, end = en$X3), # BED -> GRanges
  ep_id = en$X4
)

hits_enh <- findOverlaps(anno_gr2, enh_gr, ignore.strand = TRUE)

S_chan_enh <- unique(gr_id(anno_gr2[unique(queryHits(hits_enh))]))


actinopteriigy_cne_ov <- actinopteriigy_cne_ov %>%
  mutate(chromosome = paste0(chromosome, ".1"))
S_sheet7 <- paste0(
  as.character(actinopteriigy_cne_ov$chromosome),
  ":",
  as.character(actinopteriigy_cne_ov$start),
  "-",
  as.character(actinopteriigy_cne_ov$end)
)

# ---------- assemble sets, and (optionally) clamp to universe ----------
set_list_ids <- list(
  `actinopteriigy CNEs` = S_actinopteriigy,
  `ATAC Peaks` = S_atac,
  `Active genes nearby` = S_active,
  `YueSong overlap` = S_sheet7,
  `Chan enhancers overlap` = S_chan_enh
)

# ensure all sets are subset of the same universe (important for sanity)
set_list_ids <- lapply(set_list_ids, function(s) intersect(s, U))

# ---------- UpSet ----------
comb <- make_comb_mat(set_list_ids)
pdf('../output/upset_overlaps.pdf', width = 8, height = 6)
UpSet(
  comb,
  set_order = names(set_list_ids),
  top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
  right_annotation = upset_right_annotation(comb, add_numbers = TRUE)
)
dev.off()

# ---------- assemble sets, and (optionally) clamp to universe ----------
set_list_ids <- list(
  `actinopteriigy CNEs` = S_actinopteriigy,
  `YueSong overlap` = S_sheet7,
  `Chan enhancers overlap` = S_chan_enh
)

# ensure all sets are subset of the same universe (important for sanity)
set_list_ids <- lapply(set_list_ids, function(s) intersect(s, U))

# ---------- UpSet ----------
comb <- make_comb_mat(set_list_ids)
pdf('../output/upset_external_studies.pdf', width = 8, height = 6)
UpSet(
  comb,
  set_order = names(set_list_ids),
  top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
  right_annotation = upset_right_annotation(comb, add_numbers = TRUE)
)
dev.off()

# ---------- assemble sets, and (optionally) clamp to universe ----------
set_list_ids <- list(
  `actinopteriigy CNEs` = S_actinopteriigy,
  `ATAC Peaks` = S_atac,
  `Active genes nearby` = S_active
)

# ensure all sets are subset of the same universe (important for sanity)
set_list_ids <- lapply(set_list_ids, function(s) intersect(s, U))

# ---------- UpSet ----------
comb <- make_comb_mat(set_list_ids)
pdf('../output/upset_external_studies.pdf', width = 8, height = 6)
UpSet(
  comb,
  set_order = names(set_list_ids),
  top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
  right_annotation = upset_right_annotation(comb, add_numbers = TRUE)
)
dev.off()


library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)

# ── 1. Define a tolerance-aware overlap function ──────────────────────────────
# Returns the indices in `query` that overlap anything in `subject`,
# allowing up to `slack` bp of coordinate shift on either side.

# ── 2. Base GRanges (your universe) ───────────────────────────────────────────
# All set membership is determined by overlap with anno_gr2, not string matching.
base_gr <- unique(anno_gr2) # GRanges, one range per CNE

# ── 3. Build each set as a logical index vector over base_gr ──────────────────
SLACK <- 0L # adjust as needed — 0 for exact, 50-200 for lifted coords

in_atac <- overlapping_idx(base_gr, atac_peaks_gr, slack = SLACK)
in_active <- which(
  !is.na(mcols(base_gr)$is_active) &
    mcols(base_gr)$is_active
)
in_fin <- which(base_gr$geneId %in% filtered_genes_fin)
in_sheet7 <- overlapping_idx(
  base_gr,
  GRanges(
    seqnames = str_remove_all(actinopteriigy_cne_ov$chromosome, "\\.1$"),
    ranges = IRanges(
      start = actinopteriigy_cne_ov$start,
      end = actinopteriigy_cne_ov$end
    )
  ),
  slack = SLACK
)
in_chan_enh <- overlapping_idx(base_gr, enh_gr, slack = SLACK)

# ── Run analyses ─────────────────────────────────────────────────────────────

res_actino <- analyse_cne_universe(
  anno_gr = unique(peak_anno_list$actinopteriigy_CNE@anno),
  label = "actinopteriigy",
  atac_peaks_gr = atac_peaks_gr,
  enh_gr = enh_gr,
  yuesong_gr = yuesong_gr,
  fin_geneIds = filtered_genes_fin,
  slack = 0L
)

res_gnatho <- analyse_cne_universe(
  anno_gr = unique(peak_anno_list$gnathostomata_CNE@anno),
  label = "gnathostomata",
  atac_peaks_gr = atac_peaks_gr,
  enh_gr = enh_gr,
  yuesong_gr = yuesong_gr,
  fin_geneIds = filtered_genes_fin,
  slack = 0L
)

comb <- make_comb_mat(res_gnatho$final)
pdf('../output/upset_overlaps.pdf', width = 8, height = 6)
UpSet(
  comb,
  set_order = names(set_list_ids),
  top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
  right_annotation = upset_right_annotation(comb, add_numbers = TRUE)
)
dev.off()

##################################################
### FINAL ACTINOPTERIIGY TABLE (input to app.R) ###
##################################################
# Build a single tidy table from the universe used in the UpSet plot
# (`anno_gr2`, i.e. all non-exonic actinopteriigy CNE annotations).
#
# The column layout is identical to the existing teleost / actinopteriigy
# annotated tables (seqnames, start, end, width, strand, annotation, geneId,
# flank_geneIds, is_active, phastcons, ...) so the Shiny app can read
# the file directly. Four extra binary columns encode the upset-plot set
# memberships and let the app filter on them independently.

# 1. Coerce the universe GRanges to a data.frame (same layout as the existing
#    teleost-input file)
final_actinopteriigy_df <- as.data.frame(unique(anno_gr2))

# 2. Add the four binary membership columns (1 = member, 0 = not).
#    `in_atac`, `in_active`, `in_chan_enh`, `in_sheet7` are positional indices
#    into `anno_gr2`, defined in the upset section above.
n_anno <- length(unique(anno_gr2))

final_actinopteriigy_df$in_atac_peak <- as.integer(seq_len(n_anno) %in% in_atac)
final_actinopteriigy_df$nearby_gene_active <- as.integer(
  seq_len(n_anno) %in% in_active
)
final_actinopteriigy_df$overlaps_chan_enhancer <- as.integer(
  seq_len(n_anno) %in% in_chan_enh
)
final_actinopteriigy_df$overlaps_yuesong_cne <- as.integer(
  seq_len(n_anno) %in% in_sheet7
)

# 3. Sanity check: column sums must agree with the upset-plot index vectors
stopifnot(
  sum(final_actinopteriigy_df$in_atac_peak) == length(in_atac),
  sum(final_actinopteriigy_df$nearby_gene_active) == length(in_active),
  sum(final_actinopteriigy_df$overlaps_chan_enhancer) == length(in_chan_enh),
  sum(final_actinopteriigy_df$overlaps_yuesong_cne) == length(in_sheet7)
)

# 4. Write final table as TSV (read by app.R)
out_final_path <- '../output/actinopteriigy_cne_final_table.tsv'
write.table(
  final_actinopteriigy_df,
  out_final_path,
  sep = '\t',
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

base_gr <- unique(peak_anno_list$gnathostomata_CNE@anno[
  !str_detect(peak_anno_list$gnathostomata_CNE@anno$annotation, "Exon")
])

final_gnathostomata_df <- as.data.frame(base_gr)
n_anno <- length(base_gr)


# ── 3. Build each set as a logical index vector over base_gr ──────────────────

in_atac <- overlapping_idx(base_gr, atac_peaks_gr, slack = SLACK)
in_active <- which(
  !is.na(mcols(base_gr)$is_active) &
    mcols(base_gr)$is_active
)
in_fin <- which(base_gr$geneId %in% filtered_genes_fin)
in_sheet7 <- overlapping_idx(
  base_gr,
  GRanges(
    seqnames = str_remove_all(actinopteriigy_cne_ov$chromosome, "\\.1$"),
    ranges = IRanges(
      start = actinopteriigy_cne_ov$start,
      end = actinopteriigy_cne_ov$end
    )
  ),
  slack = SLACK
)
in_chan_enh <- overlapping_idx(base_gr, enh_gr, slack = SLACK)

final_gnathostomata_df$in_atac_peak <- as.integer(seq_len(n_anno) %in% in_atac)
final_gnathostomata_df$nearby_gene_active <- as.integer(
  seq_len(n_anno) %in% in_active
)
final_gnathostomata_df$overlaps_chan_enhancer <- as.integer(
  seq_len(n_anno) %in% in_chan_enh
)
final_gnathostomata_df$overlaps_yuesong_cne <- as.integer(
  seq_len(n_anno) %in% in_sheet7
)

# 4. Write final table as TSV (read by app.R)
out_final_path <- '../output/gnathostomata_cne_final_table.tsv'
write.table(
  final_gnathostomata_df,
  out_final_path,
  sep = '\t',
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
