### LIBRARIES ###
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  ChIPseeker,
  GenomicFeatures,
  zFPKM,
  SummarizedExperiment,
  readxl
)

setwd(this.path::here())

source('custom_functions.R')

### OUTPUT DIRECTORY ###
preproc_dir <- "../output/preprocessed"
report_dir <- "../output/reports"
dir.create(preproc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
### 1. READ CNE BED FILES                                                   ###
###############################################################################

actinopteriigy_cne_gr <- read_cne_bed(
  "../input/actinopterygii_specific_drer.bed"
)
gnathostomata_cne_gr <- read_cne_bed(
  "../input/gnathostomata_conserved_drer.bed"
)

# Raw tables are still needed for the YueSong key, which is built on cne_name.
actinopteriigy_cne <- read_tsv(
  "../input/actinopterygii_specific_drer.bed",
  col_names = c(
    "chromosome",
    "start",
    "end",
    "cne_name",
    "phastcons_score",
    "strand"
  ),
  show_col_types = FALSE
)

###############################################################################
### 2. READ CHROMOSOME SIZES AND BUILD TxDb                                 ###
###############################################################################

drer_sizes <- read_tsv(
  "../ancilliary_files/drer_chrom_info.txt",
  col_names = c("chrom", "length"),
  show_col_types = FALSE
) %>%
  filter(length > 20000) %>%
  mutate(chrom = sub("\\.[0-9]+$", "", chrom))

drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)

###############################################################################
### 3. READ SALMON TPM AND COMPUTE zFPKM-BASED GENE ACTIVITY                ###
###############################################################################

salmon_tpm <- read_tsv("../input/salmon.merged.gene_tpm.tsv")

tpm_matrix <- salmon_tpm %>%
  dplyr::select(starts_with("SRR")) %>%
  as.matrix()
rownames(tpm_matrix) <- salmon_tpm$gene_id

se <- SummarizedExperiment(
  assays = SimpleList(fpkm = tpm_matrix),
  rowData = DataFrame(
    gene_id = salmon_tpm$gene_id,
    gene_name = salmon_tpm$gene_name
  )
)
assay(se, "zfpkm") <- zFPKM(se)

salmon_tpm <- salmon_tpm %>%
  mutate(
    mean_tpm = rowMeans(dplyr::select(., starts_with("SRR")), na.rm = TRUE)
  )

gene_activity <- salmon_tpm %>%
  mutate(
    z_fpkm = assay(se, "zfpkm")[, 1],
    is_active = z_fpkm > -3
  ) %>%
  dplyr::select(gene_id, gene_name, mean_tpm, z_fpkm, is_active) %>%
  mutate(gene_id = str_remove_all(gene_id, "^[a-z]+-"))

saveRDS(gene_activity, file.path(preproc_dir, "gene_activity.rds"))

###############################################################################
### 4. INJECT SEQLENGTHS                                                    ###
###############################################################################

inject_seqlengths <- function(gr, sizes) {
  sl <- setNames(rep(NA_real_, length(seqlevels(gr))), seqlevels(gr))
  m <- match(names(sl), sizes$chrom)
  sl[!is.na(m)] <- sizes$length[m[!is.na(m)]]

  if (anyNA(sl)) {
    warning(
      "No length available for: ",
      paste(names(sl)[is.na(sl)], collapse = ", ")
    )
  }
  seqlengths(gr) <- sl
  gr
}

actinopteriigy_cne_gr <- inject_seqlengths(actinopteriigy_cne_gr, drer_sizes)
gnathostomata_cne_gr <- inject_seqlengths(gnathostomata_cne_gr, drer_sizes)

###############################################################################
### 5. IMPORT EXTERNAL DATASETS                                             ###
###############################################################################

# --- YueSong et al 2025, sheet 7 + RefSeq->UCSC alias mapping ----------------
sheet7_data <- read_excel(
  "../input/YueSong-et-al_2025.xlsx",
  sheet = 7,
  col_names = TRUE,
  skip = 1
)

alias_url <- paste0(
  "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/002/035/GCF_000002035.6/",
  "GCF_000002035.6.chromAlias.txt"
)
alias <- read_tsv(
  alias_url,
  comment = "#",
  show_col_types = FALSE,
  col_names = FALSE
)
refseq2ucsc <- setNames(alias$X5, alias$X1)
sheet7_data$Chromosome <- alias$X5[match(sheet7_data$Chromosome, alias$X1)]

sheet7_gr <- GRanges(
  seqnames = sheet7_data[[1]],
  ranges = IRanges(start = sheet7_data[[2]], end = sheet7_data[[3]]),
  mcols = sheet7_data[, 4:ncol(sheet7_data)]
)

# --- Liftover of actinopteriigy CNEs to GRCz12 for overlap with sheet 7 ------
liftover <- read_tsv(
  "../input/ucsc_GRCz11-GRCz12_liftover_actinopteriigy_cne.bed",
  col_names = FALSE,
  show_col_types = FALSE
)
liftover_gr <- GRanges(
  seqnames = liftover$X1,
  ranges = IRanges(start = liftover$X2 + 1, end = liftover$X3), # BED -> 1-based
  strand = "*",
  cne_name = liftover$X4
)
hits <- findOverlaps(liftover_gr, sheet7_gr, ignore.strand = TRUE)
overlapping_cne_names <- unique(mcols(liftover_gr)$cne_name[queryHits(hits)])
actinopteriigy_cne_ov <- actinopteriigy_cne %>%
  filter(cne_name %in% overlapping_cne_names) %>%
  mutate(start = start + 1L)

# --- Chan et al EP-loops + enhancer / promoter calls (GRCz12) ----------------
ep <- read_tsv("../input/Chan_et_al_EP-loops.tsv", show_col_types = FALSE)
en <- read_tsv(
  "../ancilliary_files/enhancer.grcz12.bed",
  col_names = FALSE,
  show_col_types = FALSE
)
pr <- read_tsv(
  "../ancilliary_files/promoter.grcz12.bed",
  col_names = FALSE,
  show_col_types = FALSE
)

enh_gr <- GRanges(
  seqnames = en$X1,
  ranges = IRanges(start = en$X2 + 1, end = en$X3),
  ep_id = en$X4
)

# --- ATAC consensus peaks ----------------------------------------------------
atac_peaks <- read_tsv(
  "../input/consensus_peaks.mLb.clN.bed",
  col_names = FALSE,
  show_col_types = FALSE
)
colnames(atac_peaks) <- c(
  "chromosome",
  "start",
  "end",
  "peak_name",
  "score",
  "strand"
)
atac_peaks_gr <- GRanges(
  seqnames = atac_peaks$chromosome,
  ranges = IRanges(start = atac_peaks$start + 1L, end = atac_peaks$end),
  peak_name = atac_peaks$peak_name,
  score = atac_peaks$score
)

###############################################################################
### 6. ChIPseeker ANNOTATION + GENE-ACTIVITY JOIN                           ###
###############################################################################

peak_anno_list <- lapply(
  list(actinopteriigy_cne_gr, gnathostomata_cne_gr),
  annotatePeak,
  overlap = "all",
  TxDb = drer_anno,
  tssRegion = c(-3000, 3000)
)
names(peak_anno_list) <- c("actinopteriigy_CNE", "gnathostomata_CNE")

# Annotation-class composition under the default priority. Reported because the
# switch away from the custom ordering moves elements between classes, and
# because Exon-annotated elements are dropped by non_exon() further down.
annotation_composition <- imap_dfr(peak_anno_list, function(obj, nm) {
  as_tibble(obj@anno) %>%
    mutate(annotation_class = str_remove(annotation, " \\(.*$")) %>%
    dplyr::count(annotation_class, name = "n") %>%
    mutate(set = nm, percent = round(100 * n / sum(n), 1))
}) %>%
  arrange(set, desc(n))

write_tsv(
  annotation_composition,
  file.path(report_dir, "annotation_composition_default_priority.tsv")
)

peak_anno_list <- lapply(peak_anno_list, function(anno_obj) {
  anno_df <- as_tibble(anno_obj@anno) %>%
    left_join(gene_activity, by = c("geneId" = "gene_id"))
  anno_obj@anno <- as(anno_df, "GRanges")
  anno_obj@anno <- non_exon(anno_obj@anno)
  anno_obj
})

# Element counts surviving the exon / expression filter, for the Methods.
filter_summary <- tibble(
  set = c("actinopterygii_specific", "gnathostomata_conserved"),
  n_input = c(length(actinopteriigy_cne_gr), length(gnathostomata_cne_gr)),
  n_retained = c(
    length(peak_anno_list$actinopteriigy_CNE@anno),
    length(peak_anno_list$gnathostomata_CNE@anno)
  )
) %>%
  mutate(percent_retained = round(100 * n_retained / n_input, 1))

write_tsv(filter_summary, file.path(report_dir, "cne_filter_summary.tsv"))
print(filter_summary)

write.table(
  as.data.frame(non_exon(peak_anno_list$gnathostomata_CNE@anno)),
  "../output/gnathostomata_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
write.table(
  as.data.frame(non_exon(peak_anno_list$actinopteriigy_CNE@anno)),
  "../output/actinopteriigy_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

###############################################################################
### 7. SARCOPTERYGII-SPECIFIC SET (HUMAN COORDINATES)                       ###
###############################################################################
# The third phastCons round is referenced to human, so this set cannot be
# projected onto the zebrafish annotation used above: tetrapod-specific
# elements are by definition largely absent from the zebrafish assembly. It is
# therefore carried in its own coordinate space and reported alongside the
# gnathostome set called on the same human reference, which is the symmetric
# comparison. Only the descriptive products (element counts, widths, chromosomal
# distribution) are shared between the two coordinate spaces.

hsap_files <- list(
  sarcopterygii = "../input/sarcopterygii_specific_hsap.bed",
  gnathostomata = "../input/gnathostomata_conserved_hsap.bed"
)
hsap_sizes_file <- "../ancilliary_files/hsap_chrom_info.txt"

if (all(file.exists(unlist(hsap_files), hsap_sizes_file))) {
  hsap_sizes <- read_tsv(
    hsap_sizes_file,
    col_names = c("chrom", "length"),
    show_col_types = FALSE
  ) %>%
    filter(length > 20000) %>%
    mutate(chrom = sub("\\.[0-9]+$", "", chrom))

  hsap_cne_gr <- lapply(hsap_files, read_cne_bed)
  hsap_cne_gr <- lapply(hsap_cne_gr, inject_seqlengths, sizes = hsap_sizes)

  saveRDS(hsap_cne_gr, file.path(preproc_dir, "hsap_cne_gr.rds"))
  saveRDS(hsap_sizes, file.path(preproc_dir, "hsap_sizes.rds"))

  hsap_summary <- imap_dfr(hsap_cne_gr, function(gr, nm) {
    tibble(
      set = nm,
      reference = "hsap",
      n_elements = length(gr),
      total_bp = sum(as.numeric(width(gr))),
      median_width = median(width(gr)),
      mean_width = round(mean(width(gr)), 1)
    )
  })
  write_tsv(hsap_summary, file.path(report_dir, "hsap_cne_summary.tsv"))
  print(hsap_summary)
} else {
  message(
    "Human-coordinate CNE BEDs not found; skipping the sarcopterygii set.\n",
    "  Expected: ",
    paste(c(unlist(hsap_files), hsap_sizes_file), collapse = ", ")
  )
}

###############################################################################
### 8. PERSIST EVERYTHING DOWNSTREAM SCRIPTS NEED                           ###
###############################################################################

# CNE GRanges (consumed by both scripts 2 and 3)
saveRDS(
  non_exon(peak_anno_list$actinopteriigy_CNE@anno),
  file.path(preproc_dir, "actinopteriigy_cne_gr.rds")
)
saveRDS(
  non_exon(peak_anno_list$gnathostomata_CNE@anno),
  file.path(preproc_dir, "gnathostomata_cne_gr.rds")
)

# Chrom-size table (used by script 2 for extendTSS)
saveRDS(drer_sizes, file.path(preproc_dir, "drer_sizes.rds"))

# TxDb (saved via SQLite path; reload with loadDb())
txdb_path <- file.path(preproc_dir, "drer_anno.sqlite")
if (file.exists(txdb_path)) {
  file.remove(txdb_path)
}
AnnotationDbi::saveDb(drer_anno, txdb_path)

# ChIPseeker output + per-CNE activity (used by script 3)
saveRDS(peak_anno_list, file.path(preproc_dir, "peak_anno_list.rds"))

# External datasets (used by script 3)
saveRDS(sheet7_gr, file.path(preproc_dir, "sheet7_gr.rds"))
saveRDS(
  actinopteriigy_cne_ov,
  file.path(preproc_dir, "actinopteriigy_cne_ov.rds")
)
saveRDS(enh_gr, file.path(preproc_dir, "chan_enhancers_gr.rds"))
saveRDS(atac_peaks_gr, file.path(preproc_dir, "atac_peaks_gr.rds"))
