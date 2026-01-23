### SCRIPT METADATA ###
# Title: Genomic Range Analysis and Annotation
# Author: Ferenc Kagan
# Date: 2024.08.24
# Description: This script performs a detailed analysis of conserved non-coding elements (CNEs) and associated genomic annotations
#              in the Macropodus opercularis genome, including data cleaning, exploratory data analysis, and enrichment analysis.

### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(
  GenomicRanges, tidyverse, rtracklayer, gprofiler2, rGREAT, regioneR,
  patchwork, RColorBrewer, Biostrings, ChIPseeker, GenomicFeatures, zFPKM, SummarizedExperiment,
  readxl, UpSetR, ComplexHeatmap
)
### SET WORKING DIRECTORY ###
setwd(this.path::here())

### READ IN DATA ###
# Read in CNE bed files
bed_colnames <- c('chromosome', 'start', 'end', 'cne_name', 'phastcons_score', 'strand')
teleost_cne <- read_tsv("../input/teleostei_cne.bed", col_names = F)
vertebrata_cne <- read_tsv("../input/vertebrata_cne.bed", col_names = F)

colnames(teleost_cne) <- bed_colnames
colnames(vertebrata_cne) <- bed_colnames

# Read in phyloP score calculations
# List files with "phyloP" suffix (adjust pattern as needed)
files <- list.files(path = '../input/phyloP/', pattern = "phyloP$", full.names = TRUE)

# Read all files into a list of data frames
data_list <- lapply(files, read.table, header = T, sep = "", stringsAsFactors = FALSE)

# Concatenate all data frames row-wise
combined_data <- do.call(rbind, data_list)

drer_sizes <- read_tsv("../ancilliary_files/drer_chrom_info.txt", col_names = FALSE)
colnames(drer_sizes) <- c("chrom", "length")

drer_anno <- txdbmaker::makeTxDbFromGFF("../ancilliary_files/drer.gff",
                                              organism = "Danio rerio")

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
  rowData = DataFrame(gene_id = salmon_tpm$gene_id, gene_name = salmon_tpm$gene_name)
)

# Compute zFPKM using the multivariate Gaussian mixture modeling approach
assay(se, "zfpkm") <- zFPKM(se)

# Calculate mean TPM across replicates
salmon_tpm <- salmon_tpm %>%
  mutate(mean_tpm = rowMeans(dplyr::select(., starts_with("SRR")), na.rm = TRUE))

# Extract zFPKM values and create gene activity lookup table
# Genes with zFPKM > -3 are considered active
gene_activity <- salmon_tpm %>%
  mutate(
    z_fpkm = assay(se, "zfpkm")[, 1],  # Use first zFPKM column (averaged across samples)
    is_active = z_fpkm > -3  # Standard zFPKM threshold for active genes
  ) %>%
  dplyr::select(gene_id, gene_name, mean_tpm, z_fpkm, is_active) %>%
  mutate(gene_id = str_remove_all(gene_id, "^[a-z]+-"))


### EXPLORATORY DATA ANALYSIS (EDA) ###
# Convert CNEs and mask data to GRanges objects
teleost_cne_gr <- GRanges(seqnames =  sub("\\.[0-9]+$", "", teleost_cne$cne_name),
                       ranges = IRanges(start = teleost_cne$start, 
                                        end = teleost_cne$end,),
                       phastcons = teleost_cne$phastcons_score)

vertebrata_cne_gr <- GRanges(seqnames = sub("\\.[0-9]+$", "", vertebrata_cne$cne_name),
                          ranges = IRanges(start = vertebrata_cne$start, 
                                           end = vertebrata_cne$end,),
                          phastcons = vertebrata_cne$phastcons_score)
crossmapping <- data.frame(nonstrip = seqlevels(teleost_cne_gr), stripped = sub('\\.[0-9]+$', '', seqlevels(teleost_cne_gr)))
drer_sizes$chrom <- crossmapping$nonstrip[match(drer_sizes$chrom, crossmapping$stripped)]

seqlengths(teleost_cne_gr)[drer_sizes$chrom] <- drer_sizes$length

### IMPORT SHEET 7 FROM YUESONG ET AL 2025 AND LIFTOVER ###
# Read sheet 7 from the Excel file (pseudo-bed format: chrom, start, end, and additional columns)
sheet7_data <- read_excel("../input/YueSong-et-al_2025.xlsx", sheet = 7, col_names = T, skip = 1)

alias_url <- "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/002/035/GCF_000002035.6/GCF_000002035.6.chromAlias.txt"

alias <- read_tsv(alias_url, comment = "#", show_col_types = FALSE, col_names = F)

# Build RefSeq -> UCSC map (NC_.... -> chrN)
refseq2ucsc <- setNames(alias$X5, alias$X1)

sheet7_data$Chromosome <- alias$X5[match(sheet7_data$Chromosome, alias$X1)]

# Convert to GRanges object (assuming first 3 columns are chrom, start, end)
sheet7_gr <- GRanges(
  seqnames = sheet7_data[[1]],
  ranges = IRanges(start = sheet7_data[[2]], end = sheet7_data[[3]]),
  mcols = sheet7_data[, 4:ncol(sheet7_data)]
)

liftover <- read_tsv('../input/ucsc_GRCz11-GRCz12_liftover_teleostei_cne.bed', col_names = F)
# IMPORTANT: attach X4 as metadata (cne_name)
liftover_gr <- GRanges(
  seqnames = liftover$X1,
  ranges   = IRanges(start = liftover$X2 + 1, end = liftover$X3),  # BED start is 0-based
  strand   = "*",
  cne_name = liftover$X4
)

hits <- findOverlaps(liftover_gr, sheet7_gr, ignore.strand = TRUE)

overlapping_cne_names <- unique(mcols(liftover_gr)$cne_name[queryHits(hits)])

teleost_cne_ov <- teleost_cne %>%
  filter(cne_name %in% overlapping_cne_names)

# Annotate peaks and create visualizations
peak_anno_list <- lapply(list(teleost_cne_gr, vertebrata_cne_gr), annotatePeak, TxDb = drer_anno,
                         tssRegion = c(-3000, 3000),
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))
names(peak_anno_list) <- c("Teleost_CNE", "Vertebrata_CNE")

# Add gene expression information to annotations
# For each annotation, find the nearest gene and add its expression status
peak_anno_list <- lapply(peak_anno_list, function(anno_obj) {
  anno_df <- as_tibble(anno_obj@anno) %>%
    # Join with gene activity info by nearest gene (geneId)
    left_join(gene_activity, by = c("geneId" = "gene_id")) %>%
    # For flank genes, also get their activity status
    separate_rows(flank_geneIds, sep = ";") %>%
    left_join(gene_activity %>% dplyr::select(gene_id, mean_tpm, z_fpkm, is_active) %>%
              rename_with(~paste0("flank_", .), starts_with(c("mean", "z_f", "is_"))),
              by = c("flank_geneIds" = "gene_id")) #%>%
    # Flag if nearest gene is active or any flank gene is active
    #mutate(
    #  nearest_gene_active = coalesce(is_active.y, FALSE),
    #  any_flank_gene_active = coalesce(flank_is_active.y, FALSE)
    #)
  
  # Update the annotation object
  anno_obj@anno <- as(anno_df, "GRanges")
  anno_obj
})


# Plot and save peak annotation distributions
p1 <- plotAnnoBar(peak_anno_list) + ggtitle(NULL)
p2 <- plotDistToTSS(peak_anno_list) + ggtitle(NULL) + ylab("CNE (%) (5' -> 3')")
combined_plot <- p1 / p2 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
ggsave("output/peak_annotation_distribution.png", combined_plot, units = 'in', width = 10, height = 10, dpi = 320)

#########################################################
### OVERLAP WITH DIFFERENTIALLY EXPRESSED GENES (DGE) ###
#########################################################

dge <- read_tsv("input/20240621-edgeR_LabOrg-Gill.tabular") %>%
  filter(FDR <= 0.05)
crossref <- read_tsv("input/pfish_macOpe2_geneID_UCSC_ann.gff", col_names = FALSE)

# Merge ANAR annotation with DGE and cross-reference data
anar_tibble <- peak_anno_list$ANAR@anno %>%
  as_tibble() %>%
  separate_rows(flank_geneIds, sep = ";") 

merged_data <- anar_tibble %>%
  left_join(dge, by = c("flank_geneIds" = "GeneID")) %>%
  left_join(crossref, by = c("flank_geneIds" = "X1")) %>%
  dplyr::select(-c("transcriptId", "flank_txIds")) %>%
  distinct()

# Write output data to file
write.table(merged_data, file = "output/anar_annotations.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Create BED file of unique annotations
asd_bed <- merged_data %>% distinct() %>% dplyr::select(1:3)
write.table(asd_bed, "output/ANAR.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Calculate genomic coverage percentages
cne_coverage <- sum(width(mope_cne_gr)) / sum(as.integer(mope_sizes$length)) * 100
anar_coverage <- sum(width(anar_gr)) / sum(as.integer(mope_sizes$length)) * 100

# Set sequence lengths for ANAR GRanges object
seqlengths(anar_gr) <- as.integer(mope_sizes$length[match(names(seqlengths(anar_gr)), mope_sizes$chrom)])

# Plot ANAR chromosomal distribution
cne_distribution_plot <- CNEr::plotCNEDistribution(teleost_cne_gr, chrScale = "Mb") +
  geom_point(size = 0.01, alpha = 0.5) +
  theme(legend.key = element_blank(), 
        text = element_text(size = 7),
        strip.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("output/ANAR_chromosomal_distribution.png", cne_distribution_plot, width = 10, height = 10, units = 'in', dpi = 300)

### ANNOTATION WIDTH DISTRIBUTION PLOT ###
# Prepare data for plotting annotation widths
plot_data1 <- peak_anno_list$Vertebrata_CNE@anno %>%
  as_tibble() %>%
  mutate(Annotation = str_replace_all(annotation, "Intron .*", "Intron"),
         Annotation = str_replace_all(annotation, "Promoter .*", "Promoter"),
         Category = "ANAR") 

plot_data2 <- peak_anno_list$Teleost_CNE@anno %>%
  as_tibble() %>%
  mutate(Annotation = str_replace_all(annotation, "Intron .*", "Intron"),
         Annotation = str_replace_all(annotation, "Promoter .*", "Promoter"), 
         Category = "CNE")

# Create and save density plot
density_plot <- bind_rows(plot_data1, plot_data2) %>%
  ggplot(aes(x = width, color = Category)) +
  geom_density(size = 2, alpha = 0.5) +
  facet_wrap(~Annotation) +
  scale_x_log10() +
  theme_bw() +
  ylab("Density") +
  xlab("log10(nt)") +
  theme(text = element_text(size = 15),
        legend.position = c(0.065, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))
ggsave("../output/cne_widths.png", density_plot, height = 10, width = 10, units = 'in', dpi = 300)

### ENRICHMENT ANALYSIS ###
cnes <- peak_anno_list$Teleost_CNE@anno %>%
  as_tibble() %>%
  group_by(geneId) %>%
  summarise(CNE = n())

anars <- peak_anno_list$Vertebrata_CNE@anno %>%
  as_tibble() %>%
  group_by(geneId) %>%
  summarise(ANAR = n())

enrichment_results <- full_join(cnes, anars) %>%
  mutate(CNE_genome = length(teleost_cne_gr),
         ANAR_genome = length(vertebrata_cne_gr)) %>%
  mutate(enrich = phyper(ANAR - 1, ANAR_genome, CNE_genome - ANAR_genome, CNE, lower.tail = FALSE)) %>%
  mutate(enrich = p.adjust(enrich)) %>%
  filter(enrich <= 0.05) %>%
  left_join(crossref, by = c("geneId" = "X1")) %>%
  left_join(dge, by = c("geneId" = "GeneID"))

write.table(enrichment_results, "output/genes_enriched_in_ANARs.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#####################
### GO ernichment ###
#####################

# Load libraries
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(enrichplot)

gene_list <- peak_anno_list$Teleost_CNE@anno$geneId[!(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, "Exon"))]
gene_list <- unique(gene_list)

# If your genes are symbols, convert to Entrez IDs (recommended)
gene_ids <- bitr(gene_list, 
                 fromType = "SYMBOL", 
                 toType = c("ENTREZID", "ENSEMBL"), 
                 OrgDb = org.Dr.eg.db)

head(gene_ids)

# Using gene symbols directly
ego_bp <- enrichGO(gene = gene_list,
                   OrgDb = org.Dr.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   readable = TRUE)

# View results
head(ego_bp@result)

signif_go <- ego_bp@result[ego_bp@result$qvalue <= 0.05, ]
signif_go_filtered_skeletal_fin <- signif_go[signif_go$Description %in% c("fin morphogenesis", "fin development", "skeletal system development", "muscle structure development"), ]
signif_go_filtered_fin <- signif_go[signif_go$Description %in% c("fin morphogenesis", "fin development", "skeletal system development", "muscle structure development"), ]

filtered_genes <- unique(unlist(strsplit(signif_go_filtered_fin$geneID, "/")))

peak_anno_list$Teleost_CNE@anno[ peak_anno_list$Teleost_CNE@anno$geneId %in% filtered_genes, ]
teleost_filtered <- as.data.frame(peak_anno_list$Teleost_CNE@anno[ 
  (peak_anno_list$Teleost_CNE@anno$geneId %in% filtered_genes) & 
    !(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, 'Exon')),  ])


# Clean annotation and count
df <- teleost_filtered %>%
  mutate(annotation = str_remove(annotation, " \\(.*$")) %>%
  count(annotation) %>%
  mutate(perc = n / sum(n) * 100,
         label = paste0(round(perc, 1), "%"))

# Donut chart
ggplot(df, aes(x = 2, y = n, fill = annotation)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # creates the "hole" in the middle
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), size = 4, color = 'white') +
  theme_void() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Annotation")) +
  scale_fill_brewer(palette = 'Dark2')
ggsave('../output/GO-filtered_cne_annotations.pdf')


##############
### EXPORT ###
##############

export <- peak_anno_list$Teleost_CNE@anno[!(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, 'Exon')), ]
signif_go_filtered_skeletal_fin <- signif_go[signif_go$Description %in% c("fin morphogenesis", "fin development", "skeletal system development", "muscle structure development"), ]
signif_go_filtered_fin <- signif_go[signif_go$Description %in% c("fin morphogenesis", "fin development"), ]

filtered_genes_skeletal_fin <- unique(unlist(strsplit(signif_go_filtered_skeletal_fin$geneID, "/")))
filtered_genes_fin <- unique(unlist(strsplit(signif_go_filtered_fin$geneID, "/")))

teleost_filtered_skeletal_fin <- as.data.frame(peak_anno_list$Teleost_CNE@anno[ 
  (peak_anno_list$Teleost_CNE@anno$geneId %in% filtered_genes_skeletal_fin) & 
    !(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, 'Exon')),  ])
teleost_filtered_fin <- as.data.frame(peak_anno_list$Teleost_CNE@anno[ 
  (peak_anno_list$Teleost_CNE@anno$geneId %in% filtered_genes_fin) & 
    !(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, 'Exon')),  ])

write.table(teleost_filtered_fin, '../output/teleost_GO-fin_specific_cne.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
write.table(teleost_filtered_skeletal_fin, '../output/teleost_GO-fin-skeletal_specific_cne.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
write.table(export, '../output/teleost_specific_cne.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

go_signif <- ego_bp@result[ego_bp@result$qvalue <= 0.05, ]

write.table(go_signif, '../output/GO_enrichment_results.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

### ACTIVE GENE FILTERING ###

# Export CNEs near active genes (all annotations, non-exonic)
export_active <- as.data.frame(peak_anno_list$Teleost_CNE@anno[
  (peak_anno_list$Teleost_CNE@anno$nearest_gene_active == TRUE) &
    !(str_detect(peak_anno_list$Teleost_CNE@anno$annotation, 'Exon')),  ])

if(nrow(export_active) > 0) {
  write.table(export_active, '../output/teleost_cne_near_active_genes.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
  cat("\nCNEs near active genes (non-exonic):", nrow(export_active), "\n")
} else {
  cat("\nNo CNEs near active genes found.\n")
}

# Export GO-filtered CNEs that are near active genes (fin-specific)
teleost_filtered_fin_active <- teleost_filtered_fin %>%
  filter(flank_is_active == TRUE)

if(nrow(teleost_filtered_fin_active) > 0) {
  write.table(teleost_filtered_fin_active, '../output/teleost_GO-fin_specific_cne_active_genes.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
  cat("Fin-specific CNEs near active genes:", nrow(teleost_filtered_fin_active), "\n")
} else {
  cat("No fin-specific CNEs near active genes found.\n")
}

# Export GO-filtered CNEs that are near active genes (fin+skeletal-specific)
teleost_filtered_skeletal_fin_active <- teleost_filtered_skeletal_fin %>%
  filter(flank_is_active == TRUE)

if(nrow(teleost_filtered_skeletal_fin_active) > 0) {
  write.table(teleost_filtered_skeletal_fin_active, '../output/teleost_GO-fin-skeletal_specific_cne_active_genes.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
  cat("Fin+skeletal-specific CNEs near active genes:", nrow(teleost_filtered_skeletal_fin_active), "\n")
} else {
  cat("No fin+skeletal-specific CNEs near active genes found.\n")
}

############################
### ATACseq Peak Overlap ###
############################

# Read in ATACseq consensus peaks BED file
atac_peaks <- read_tsv("../input/consensus_peaks.mLb.clN.bed", col_names = FALSE)
colnames(atac_peaks) <- c('chromosome', 'start', 'end', 'peak_name', 'score', 'strand')

# Convert ATACseq peaks to GRanges object
atac_peaks_gr <- GRanges(
  seqnames = atac_peaks$chromosome,
  ranges = IRanges(start = atac_peaks$start, end = atac_peaks$end),
  peak_name = atac_peaks$peak_name,
  score = atac_peaks$score
)

# Find overlaps between teleost CNEs and ATACseq peaks
cne_atac_overlaps <- findOverlaps(teleost_cne_gr, atac_peaks_gr)

# Extract overlapping CNEs and peaks
overlapping_cnes <- teleost_cne_gr[queryHits(cne_atac_overlaps)]
overlapping_peaks <- atac_peaks_gr[subjectHits(cne_atac_overlaps)]

# Create a summary table of overlaps
overlap_summary <- data.frame(
  cne_chr = seqnames(overlapping_cnes),
  cne_start = start(overlapping_cnes),
  cne_end = end(overlapping_cnes),
  cne_phastcons = overlapping_cnes$phastcons,
  atac_chr = seqnames(overlapping_peaks),
  atac_start = start(overlapping_peaks),
  atac_end = end(overlapping_peaks),
  atac_peak_name = overlapping_peaks$peak_name,
  atac_score = overlapping_peaks$score
)

# Calculate overlap statistics
n_cnes_total <- length(teleost_cne_gr)
n_atac_total <- length(atac_peaks_gr)
n_cnes_overlapping <- length(unique(queryHits(cne_atac_overlaps)))
n_atac_overlapping <- length(unique(subjectHits(cne_atac_overlaps)))
pct_cnes_overlapping <- (n_cnes_overlapping / n_cnes_total) * 100
pct_atac_overlapping <- (n_atac_overlapping / n_atac_total) * 100

# Print overlap statistics
cat("\n=== ATACseq Peak vs Teleost CNE Overlap Statistics ===\n")
cat("Total CNEs:", n_cnes_total, "\n")
cat("Total ATACseq peaks:", n_atac_total, "\n")
cat("CNEs with overlap:", n_cnes_overlapping, paste0("(", round(pct_cnes_overlapping, 2), "%)\n"))
cat("ATACseq peaks with overlap:", n_atac_overlapping, paste0("(", round(pct_atac_overlapping, 2), "%)\n"))
cat("Total overlap instances:", nrow(overlap_summary), "\n")

# Export overlap results
write.table(overlap_summary, '../output/cne_atac_overlaps.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

# Export overlapping CNEs with annotation (non-exonic) and gene expression info using proper overlap logic
anno_gr <- peak_anno_list$Teleost_CNE@anno
anno_atac_overlaps <- findOverlaps(anno_gr, atac_peaks_gr)

# Subset annotated CNEs that overlap ATACseq peaks and are non-exonic
cne_atac_annotated <- anno_gr[queryHits(anno_atac_overlaps)]
cne_atac_annotated <- cne_atac_annotated[!(str_detect(cne_atac_annotated$annotation, 'Exon'))]

cne_atac_annotated_df <- as.data.frame(cne_atac_annotated)
write.table(cne_atac_annotated_df, '../output/teleost_cne_atac_overlaps_annotated.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

# Export ATACseq-overlapping CNEs from ACTIVE genes only
cne_atac_annotated_active <- cne_atac_annotated_df %>%
  filter(flank_is_active == TRUE)

if(nrow(cne_atac_annotated_active) > 0) {
  write.table(cne_atac_annotated_active, '../output/teleost_cne_atac_overlaps_active_genes.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
  cat("\nATACseq-overlapping CNEs from active genes:", nrow(cne_atac_annotated_active), "\n")
} else {
  cat("\nNo ATACseq-overlapping CNEs from active genes found.\n")
}

# Create visualization of overlap distribution
overlap_by_annotation <- cne_atac_annotated_df %>%
  as_tibble() %>%
  mutate(annotation = str_remove(annotation, " \\(.*$")) %>%
  dplyr::count(annotation) %>%
  mutate(perc = n / sum(n) * 100) %>%
  arrange(desc(n))

overlap_plot <- ggplot(overlap_by_annotation, aes(x = reorder(annotation, -n), y = n, fill = annotation)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(perc, 1), "%)")), vjust = -0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Annotation Type", y = "Count", title = "Distribution of ATACseq Peak-overlapping CNEs by Annotation")

ggsave("../output/cne_atac_overlaps_by_annotation.png", overlap_plot, width = 8, height = 6, units = 'in', dpi = 300)

# Filter ATACseq-overlapping CNEs to fin-specific genes only
cne_atac_fin_specific <- cne_atac_annotated_df %>%
  filter(geneId %in% filtered_genes_fin)

# Export fin-specific ATACseq overlap results
write.table(cne_atac_fin_specific, '../output/teleost_cne_atac_overlaps_fin_specific.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

# Filter to fin-specific AND active genes
cne_atac_fin_specific_active <- cne_atac_fin_specific %>%
  filter(flank_is_active == TRUE)

if(nrow(cne_atac_fin_specific_active) > 0) {
  write.table(cne_atac_fin_specific_active, '../output/teleost_cne_atac_overlaps_fin_specific_active_genes.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
  cat("ATACseq-overlapping fin-specific CNEs from active genes:", nrow(cne_atac_fin_specific_active), "\n")
} else {
  cat("No ATACseq-overlapping fin-specific CNEs from active genes found.\n")
}

# Calculate fin-specific overlap statistics
n_atac_fin_specific <- nrow(cne_atac_fin_specific)
pct_atac_fin_specific <- (n_atac_fin_specific / nrow(cne_atac_annotated_df)) * 100 

cat("\n=== Fin-Specific ATACseq-overlapping CNEs ===\n")
cat("Fin-specific CNEs with ATACseq overlap:", n_atac_fin_specific, paste0("(", round(pct_atac_fin_specific, 2), "% of all ATACseq-overlapping CNEs)\n"))

# Create visualization of fin-specific overlap distribution by annotation
if(nrow(cne_atac_fin_specific) > 0) {
  fin_overlap_by_annotation <- cne_atac_fin_specific %>%
    as_tibble() %>%
    mutate(annotation = str_remove(annotation, " \\(.*$")) %>%
    count(annotation) %>%
    mutate(perc = n / sum(n) * 100) %>%
    arrange(desc(n))
  
  fin_overlap_plot <- ggplot(fin_overlap_by_annotation, aes(x = reorder(annotation, -n), y = n, fill = annotation)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", round(perc, 1), "%)")), vjust = -0.5, size = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = "Annotation Type", y = "Count", title = "ATACseq-overlapping Fin-Specific CNEs by Annotation")
  
  ggsave("../output/cne_atac_overlaps_fin_specific_by_annotation.png", fin_overlap_plot, width = 8, height = 6, units = 'in', dpi = 300)
} else {
  cat("No fin-specific CNEs found with ATACseq overlap.\n")
}

#############################
### Teleost-unique + ATAC ###
#############################

# Select CNEs that are unique to teleosts (no overlap with vertebrata CNEs)
if (exists("teleost_cne_gr") && exists("vertebrata_cne_gr")) {
  ov_tel_verte <- findOverlaps(teleost_cne_gr, vertebrata_cne_gr)
  teleost_only_idx <- setdiff(seq_along(teleost_cne_gr), unique(queryHits(ov_tel_verte)))
  teleost_only_gr <- teleost_cne_gr[teleost_only_idx]

  if (length(teleost_only_gr) == 0) {
    message("No teleost-unique CNEs found (all teleost CNEs overlap vertebrata set).")
  } else if (!exists("atac_peaks_gr") || length(atac_peaks_gr) == 0) {
    message("ATAC peaks not available; cannot subset teleost-unique CNEs by ATAC.")
  } else {
    # Find teleost-unique CNEs that overlap ATAC peaks
    ov_tel_atac <- findOverlaps(teleost_only_gr, atac_peaks_gr)
    if (length(ov_tel_atac) == 0) {
      message("No teleost-unique CNEs overlap ATAC peaks.")
    } else {
      tel_only_with_atac <- teleost_only_gr[unique(queryHits(ov_tel_atac))]

      # Map annotated entries (from peak_anno_list$Teleost_CNE@anno) to these ranges
      anno_gr <- peak_anno_list$Teleost_CNE@anno
      ov_anno_telonly <- findOverlaps(anno_gr, tel_only_with_atac)
      if (length(ov_anno_telonly) == 0) {
        message("Annotated teleost peaks do not match teleost-unique ATAC-overlapping ranges.")
      } else {
        final_anno_tel <- anno_gr[unique(queryHits(ov_anno_telonly))]
        final_df_tel <- as.data.frame(final_anno_tel)
        out_path <- '../output/teleost_unique_cne_with_atac_annotated.tsv'
        write.table(final_df_tel, out_path, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
        message(sprintf("Wrote %d teleost-unique, ATAC-overlapping annotated CNEs to %s", nrow(final_df_tel), out_path))
        
        # Filter to only those near active genes
        final_df_tel_active <- final_df_tel %>%
          filter(flank_is_active == TRUE)
        
        if(nrow(final_df_tel_active) > 0) {
          out_path_active <- '../output/teleost_unique_cne_with_atac_annotated_active_genes.tsv'
          write.table(final_df_tel_active, out_path_active, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
          message(sprintf("Wrote %d teleost-unique, ATAC-overlapping CNEs from ACTIVE genes to %s", nrow(final_df_tel_active), out_path_active))
        } else {
          message("No teleost-unique ATAC-overlapping CNEs from active genes found.")
        }
      }
    }
  }
}

### GENE ACTIVITY SUMMARY ###
cat("\n=== GENE ACTIVITY CLASSIFICATION SUMMARY ===\n")
cat("Total genes analyzed:", nrow(gene_activity), "\n")
cat("Active genes (z-FPKM > 0):", sum(gene_activity$is_active), "\n")
cat("Inactive genes (z-FPKM <= 0):", sum(!gene_activity$is_active), "\n")
cat("Mean TPM range:", round(min(gene_activity$mean_tpm), 4), "-", round(max(gene_activity$mean_tpm), 4), "\n")
cat("Mean z-FPKM range:", round(min(gene_activity$z_fpkm), 4), "-", round(max(gene_activity$z_fpkm), 4), "\n\n")

###################################
### UpSet Plot - CNE Categories ###
###################################
library(GenomicRanges)
library(dplyr)
library(stringr)
library(ComplexHeatmap)

# ---------- helper: stable IDs ----------
gr_id <- function(gr, include_strand = FALSE) {
  if (!include_strand) {
    paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  } else {
    paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr), ":", as.character(strand(gr)))
  }
}

# ---------- Universe ----------
anno_gr <- peak_anno_list$Teleost_CNE@anno

# optional filter
anno_gr2 <- anno_gr[!str_detect(anno_gr$annotation, "Exon")]

U <- unique(gr_id(anno_gr2))  # universe IDs

# ---------- Set 1: Teleost CNEs ----------
S_teleost <- U

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

anno_gr_active <- subset_gr_true(anno_gr2, "flank_is_active")
S_active <- unique(gr_id(anno_gr_active))

# ---------- Set 4: fin developmental genes nearby ----------
S_fin <- unique(gr_id(anno_gr2[anno_gr2$geneId %in% filtered_genes_fin]))

# ---------- Set 5: YueSong Sheet7 overlap (teleost_cne_ov) ----------
teleost_cne_ov <- teleost_cne_ov %>% mutate(chromosome = paste0(chromosome, ".1"))
S_sheet7 <- paste0(as.character(teleost_cne_ov$chromosome), ":", 
                   as.character(teleost_cne_ov$start), "-", 
                   as.character(teleost_cne_ov$end))

# ---------- Set 5: Chan el al.  overlap (teleost_cne_ov) ----------
ep <- read_tsv('../input/Chan_et_al_EP-loops.tsv')
en <- read_tsv('../ancilliary_files/enhancer.grcz12.bed', col_names = F)
pr <- read_tsv('../ancilliary_files/promoter.grcz12.bed', col_names = F)

enh_gr <- GRanges(
  seqnames = en$X1,
  ranges   = IRanges(start = en$X2 + 1, end = en$X3),  # BED -> GRanges
  ep_id    = en$X4
)

hits_enh <- findOverlaps(anno_gr2, enh_gr, ignore.strand = TRUE)

S_chan_enh <- unique(gr_id(anno_gr2[unique(queryHits(hits_enh))]))


teleost_cne_ov <- teleost_cne_ov %>% mutate(chromosome = paste0(chromosome, ".1"))
S_sheet7 <- paste0(as.character(teleost_cne_ov$chromosome), ":", 
                   as.character(teleost_cne_ov$start), "-", 
                   as.character(teleost_cne_ov$end))

# ---------- assemble sets, and (optionally) clamp to universe ----------
set_list_ids <- list(
  `Teleost CNEs` = S_teleost,
  `ATAC Peaks` = S_atac,
  `Active genes nearby` = S_active,
  `Fin-dev genes nearby` = S_fin,
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

