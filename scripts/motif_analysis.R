### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(
  GenomicRanges, tidyverse, rtracklayer, gprofiler2, rGREAT, regioneR,
  patchwork, RColorBrewer, Biostrings, ChIPseeker, GenomicFeatures, universalmotif,
  readxl, UpSetR, ComplexHeatmap, GenometriCorr
)
### SET WORKING DIRECTORY ###
setwd(this.path::here())

source('custom_functions.R')

### READ IN DATA ###

gene_activity <- readRDS('../output/gene_activity.RDS')

files <- list.files("../output/granges_rds/", full.names = TRUE)

gr <- setNames(
  lapply(files, readRDS),
  stringr::str_remove(basename(files), "\\.rds$")
)

teleost <- read_tsv('../output/meme_analysis/teleostei_motifs/ame_output/ame.tsv')
vertebrata <- read_tsv('../output/meme_analysis/vertebrata_motifs/ame_output/ame.tsv')

teleost_seq <- read_tsv('../output/meme_analysis/teleostei_motifs/ame_output/sequences.tsv')
vertebrata_seq <- read_tsv('../output/meme_analysis/vertebrata_motifs/ame_output/sequences.tsv')

tf_link <- read_table('../ancilliary_files/tflink.tsv', col_names = T)

###############
### GRANGES ###
###############

teleost_anno_gr <- gr$peak_annotations$Teleost_CNE@anno

teleost_anno_df <- as_tibble(teleost_anno_gr) %>%
  transmute(
    cne_chr   = as.character(seqnames),
    cne_start = start,
    cne_end   = end,
    
    annotation,
    distanceToTSS,
    
    geneId,
    transcriptId,
    gene_name,
    mean_tpm,
    z_fpkm,
    is_active,
    
    flank_txIds,
    flank_geneIds,
    flank_gene_distances,
    flank_mean_tpm,
    flank_z_fpkm,
    flank_is_active
  ) %>%
  distinct()

teleost_seq_tp <- teleost_seq %>%
  filter(!str_detect(seq_ID, "shuf"), class == "tp") %>%
  mutate(
    chr   = paste0(str_split_i(seq_ID, ":", 1), ".1"),
    start = as.integer(str_split_i(str_split_i(seq_ID, ":", 2), "-", 1)),
    end   = as.integer(str_split_i(str_split_i(seq_ID, ":", 2), "-", 2)),
    hit_id = paste0("teleost_hit_", row_number())
  )

teleost_hits_gr <- GRanges(
  seqnames = teleost_seq_tp$chr,
  ranges   = IRanges(start = teleost_seq_tp$start, end = teleost_seq_tp$end),
  hit_id   = teleost_seq_tp$hit_id,
  seq_id   = teleost_seq_tp$seq_ID,
  motif_id = teleost_seq_tp$motif_ID,
  motif_alt= teleost_seq_tp$motif_ALT_ID,
  fasta_score = teleost_seq_tp$FASTA_score,
  pwm_score   = teleost_seq_tp$PWM_score
)

hits_ov <- findOverlaps(gr$teleost_cne_gr, teleost_hits_gr, ignore.strand = TRUE)

teleost_cne_motif_df <- tibble(
  cne_idx = queryHits(hits_ov),
  hit_idx = subjectHits(hits_ov)
) %>%
  mutate(
    cne_id   = gr$teleost_cne_gr$cne_id[cne_idx],
    cne_chr  = as.character(seqnames(gr$teleost_cne_gr))[cne_idx],
    cne_start= start(gr$teleost_cne_gr)[cne_idx],
    cne_end  = end(gr$teleost_cne_gr)[cne_idx],
    
    motif_id = teleost_hits_gr$motif_id[hit_idx],
    motif_alt= teleost_hits_gr$motif_alt[hit_idx],
    fasta_score = teleost_hits_gr$fasta_score[hit_idx],
    pwm_score   = teleost_hits_gr$pwm_score[hit_idx]
  ) %>%
  dplyr::select(-cne_idx, -hit_idx) %>%
  distinct()

teleost_cne_motif_gene_df <- teleost_cne_motif_df %>%
  left_join(
    teleost_anno_df,
    by = c("cne_chr", "cne_start", "cne_end")
  )

write_csv(teleost_cne_motif_gene_df, '../output/cne-x-motif-x-gene.csv')

##########################
### TF LINK VALIDATION ###
##########################

norm_tf <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "") %>%   # remove punctuation/underscores/dashes
    str_replace("_[a-z]+$", "")             # if you still have suffixes like _a, _b
}

tf_link2 <- tf_link %>%
  mutate(
    tf = norm_tf(Name.TF),
    target_sym = norm_tf(Name.Target),
    target_ncbi = na_if(NCBI.GeneID.Target, "-")
  )

teleost_edges <- teleost_cne_motif_gene_df %>%
  transmute(
    tf = norm_tf(motif_alt),
    target_sym = norm_tf(geneId),  # if geneId is symbol-like in your table; else use gene_name
    cne_chr, cne_start, cne_end,
    distanceToTSS,
    is_active,
    pwm_score
  )

teleost_edges_f <- teleost_edges %>%
  filter(is_active) %>%
  distinct(tf, target_sym, .keep_all = FALSE)

tfl_edges <- tf_link2 %>%
  distinct(tf, target_sym)

validated_hits <- inner_join(
  teleost_edges_f,
  tfl_edges,
  by = c("tf", "target_sym")
)

n_pred <- nrow(teleost_edges_f)
n_val  <- nrow(tfl_edges)
n_hit  <- nrow(validated_hits)

c(n_pred = n_pred, n_val = n_val, n_overlap = n_hit)

##################
### VERTEBRATE ###
##################

vertebrata_anno_gr <- gr$peak_annotations$Vertebrata_CNE@anno

vertebrata_anno_df <- as_tibble(vertebrata_anno_gr) %>%
  transmute(
    cne_chr   = as.character(seqnames),
    cne_start = start,
    cne_end   = end,
    
    annotation,
    distanceToTSS,
    
    geneId,
    transcriptId,
    gene_name,
    mean_tpm,
    z_fpkm,
    is_active,
    
    flank_txIds,
    flank_geneIds,
    flank_gene_distances,
    flank_mean_tpm,
    flank_z_fpkm,
    flank_is_active
  ) %>%
  distinct()

vertebrata_seq_tp <- vertebrata_seq %>%
  filter(!str_detect(seq_ID, "shuf"), class == "tp") %>%
  mutate(
    chr   = paste0(str_split_i(seq_ID, ":", 1), ".1"),
    start = as.integer(str_split_i(str_split_i(seq_ID, ":", 2), "-", 1)),
    end   = as.integer(str_split_i(str_split_i(seq_ID, ":", 2), "-", 2)),
    hit_id = paste0("vertebrata_hit_", row_number())
  )

vertebrata_hits_gr <- GRanges(
  seqnames = vertebrata_seq_tp$chr,
  ranges   = IRanges(start = vertebrata_seq_tp$start, end = vertebrata_seq_tp$end),
  hit_id   = vertebrata_seq_tp$hit_id,
  seq_id   = vertebrata_seq_tp$seq_ID,
  motif_id = vertebrata_seq_tp$motif_ID,
  motif_alt= vertebrata_seq_tp$motif_ALT_ID,
  fasta_score = vertebrata_seq_tp$FASTA_score,
  pwm_score   = vertebrata_seq_tp$PWM_score
)

hits_ov <- findOverlaps(gr$vertebrata_cne_gr, vertebrata_hits_gr, ignore.strand = TRUE)

vertebrata_cne_motif_df <- tibble(
  cne_idx = queryHits(hits_ov),
  hit_idx = subjectHits(hits_ov)
) %>%
  mutate(
    cne_id   = gr$vertebrata_cne_gr$cne_id[cne_idx],
    cne_chr  = as.character(seqnames(gr$vertebrata_cne_gr))[cne_idx],
    cne_start= start(gr$vertebrata_cne_gr)[cne_idx],
    cne_end  = end(gr$vertebrata_cne_gr)[cne_idx],
    
    motif_id = vertebrata_hits_gr$motif_id[hit_idx],
    motif_alt= vertebrata_hits_gr$motif_alt[hit_idx],
    fasta_score = vertebrata_hits_gr$fasta_score[hit_idx],
    pwm_score   = vertebrata_hits_gr$pwm_score[hit_idx]
  ) %>%
  dplyr::select(-cne_idx, -hit_idx) %>%
  distinct()

vertebrata_cne_motif_gene_df <- vertebrata_cne_motif_df %>%
  left_join(
    vertebrata_anno_df,
    by = c("cne_chr", "cne_start", "cne_end")
  )

##########################
### TF LINK VALIDATION ###
##########################

vertebrata_edges <- vertebrata_cne_motif_gene_df %>%
  transmute(
    tf = norm_tf(motif_alt),
    target_sym = norm_tf(geneId),  # if geneId is symbol-like in your table; else use gene_name
    cne_chr, cne_start, cne_end,
    distanceToTSS,
    is_active,
    pwm_score
  )

vertebrata_edges_f <- vertebrata_edges %>%
  filter(is_active) %>%
  distinct(tf, target_sym, .keep_all = FALSE)

validated_hits <- inner_join(
  vertebrata_edges_f,
  tfl_edges,
  by = c("tf", "target_sym")
)

n_pred <- nrow(vertebrata_edges_f)
n_val  <- nrow(tfl_edges)
n_hit  <- nrow(validated_hits)

c(n_pred = n_pred, n_val = n_val, n_overlap = n_hit)

############
### MISC ###
############
library(JASPAR)
library(RSQLite)
library(TFBSTools)
library(universalmotif)

JASPAR <- JASPAR()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR))
siteList <- TFBSTools::getMatrixSet(JASPARConnect, 
                                    list(tax_group = "vertebrates"))
siteList
dbDisconnect(JASPARConnect)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

jaspar_meta <- imap_dfr(siteList, ~{
  m <- .x
  tg <- TFBSTools::tags(m)
  
  tibble(
    motif_id    = TFBSTools::ID(m),          # "MA0004.1"
    jaspar_name = TFBSTools::name(m),        # "Arnt"
    symbol      = tg$symbol %||% NA_character_,  # "ARNT"
    family      = tg$family %||% NA_character_,  # "PAS domain factors"
    class       = tryCatch(TFBSTools::MatrixClass(m), error = function(e) NA_character_),
    collection  = tg$collection %||% NA_character_,
    type        = tg$type %||% NA_character_,
    tax_group   = tg$tax_group %||% NA_character_,
    acc         = tg$acc %||% NA_character_,
    alias       = tg$alias %||% NA_character_
  )
})

jaspar_meta

jasp <- read_meme('../ancilliary_files/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme')
jaspar_family_map <- setNames(jaspar_meta$family, jaspar_meta$motif_id)

jasp_annot <- lapply(jasp, function(m) {
  id <- m@name
  
  fam <- jaspar_family_map[id]          # single-element named vector or NA
  fam <- unname(fam)
  
  if (!is.na(fam) && length(fam) == 1) {
    slot(m, "family") <- fam
  }
  m
})

teleost_motif <- filter_motifs(jasp_annot, name = teleost$motif_ID)
vertebrate_motif <- filter_motifs(jasp_annot, name = vertebrata$motif_ID)

teleost_tree <- motif_tree(teleost_motif, layout = 'fan', labels = 'altname')
teleost_tree

vertebrate_tree <- motif_tree(vertebrate_motif, layout = 'fan', labels = 'altname')
vertebrate_tree
