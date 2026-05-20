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
dir.create(preproc_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
### 1. READ CNE BED FILES                                                   ###
###############################################################################

bed_colnames <- c(
  "chromosome",
  "start",
  "end",
  "cne_name",
  "phastcons_score",
  "strand"
)

actinopteriigy_cne <- read_tsv(
  "../input/actinopterygii_specific_drer.bed",
  col_names = FALSE
)
gnathostomata_cne <- read_tsv(
  "../input/gnathostomata_conserved_drer.bed",
  col_names = FALSE
)
colnames(actinopteriigy_cne) <- bed_colnames
colnames(gnathostomata_cne) <- bed_colnames

###############################################################################
### 2. READ CHROMOSOME SIZES AND BUILD TxDb                                 ###
###############################################################################

drer_sizes <- read_tsv(
  "../ancilliary_files/drer_chrom_info.txt",
  col_names = FALSE
) %>%
  filter(X2 > 20000)
colnames(drer_sizes) <- c("chrom", "length")

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
### 4. BUILD CNE GRanges WITH PROPER SEQINFO                                ###
###############################################################################

# Build CNE GRanges. Seqnames come from `cne_name` with version stripped,
# matching downstream usage.
actinopteriigy_cne_gr <- GRanges(
  seqnames = sub("\\.[0-9]+$", "", actinopteriigy_cne$cne_name),
  ranges = IRanges(actinopteriigy_cne$start, actinopteriigy_cne$end),
  phastcons = actinopteriigy_cne$phastcons_score
)
gnathostomata_cne_gr <- GRanges(
  seqnames = sub("\\.[0-9]+$", "", gnathostomata_cne$cne_name),
  ranges = IRanges(gnathostomata_cne$start, gnathostomata_cne$end),
  phastcons = gnathostomata_cne$phastcons_score
)

# Width filter
actinopteriigy_cne_gr <- actinopteriigy_cne_gr[
  width(actinopteriigy_cne_gr) > 25
]
gnathostomata_cne_gr <- gnathostomata_cne_gr[
  width(gnathostomata_cne_gr) > 25
]

# --- Reconcile chromosome naming and inject seqlengths -----------------------
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
seqlengths(actinopteriigy_cne_gr) <- sl

sl_base <- sub("\\..*", "", names(seqlengths(gnathostomata_cne_gr)))

sl <- seqlengths(gnathostomata_cne_gr)
names(sl) <- sl_base # temporarily strip versions

sl[drer_sizes$chrom_base] <- drer_sizes$length

names(sl) <- seqlevels(gnathostomata_cne_gr) # restore original names
seqlengths(gnathostomata_cne_gr) <- sl

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
  col_names = FALSE
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
  filter(cne_name %in% overlapping_cne_names)

# --- Chan et al EP-loops + enhancer / promoter calls (GRCz12) ----------------
ep <- read_tsv("../input/Chan_et_al_EP-loops.tsv")
en <- read_tsv("../ancilliary_files/enhancer.grcz12.bed", col_names = FALSE)
pr <- read_tsv("../ancilliary_files/promoter.grcz12.bed", col_names = FALSE)

enh_gr <- GRanges(
  seqnames = en$X1,
  ranges = IRanges(start = en$X2 + 1, end = en$X3),
  ep_id = en$X4
)

# --- ATAC consensus peaks ----------------------------------------------------
atac_peaks <- read_tsv(
  "../input/consensus_peaks.mLb.clN.bed",
  col_names = FALSE
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
  ranges = IRanges(start = atac_peaks$start, end = atac_peaks$end),
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

peak_anno_list <- lapply(peak_anno_list, function(anno_obj) {
  anno_df <- as_tibble(anno_obj@anno) %>%
    left_join(gene_activity, by = c("geneId" = "gene_id"))
  anno_obj@anno <- as(anno_df, "GRanges")
  anno_obj@anno <- non_exon(anno_obj@anno)
  anno_obj
})

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
### 7. PERSIST EVERYTHING DOWNSTREAM SCRIPTS NEED                           ###
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
