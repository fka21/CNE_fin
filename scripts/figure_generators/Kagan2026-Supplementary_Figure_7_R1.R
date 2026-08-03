# ============================================================================
# TOP_N flanking-gene CNE hotspot summary, with GO annotation (BP / MF / CC)
#
# For each of the two CNE datasets (actinopterygii-specific and gnathostome-
# conserved):
#   1. assign every CNE to its gene via the geneId column of the input table;
#   2. OPTIONALLY drop genes whose name starts with "loc" or "si"
#      (uncharacterised / clone-based identifiers) - see DROP_UNCHARACTERISED;
#   3. rank the remaining genes by number of CNEs, take the top 25;
#   4. annotate each top gene with its key GO terms from ALL THREE ontologies
#      - Biological Process, Molecular Function, Cellular Component;
#   5. DROP any gene for which no GO term is found in any ontology;
#   6. write a summary TSV per dataset;
#   7. draw a summary lollipop visualization, faceted by dataset.
#
# Inputs (working directory):
#   actinopterygii_cne_final_table.tsv
#   gnathostomata_cne_final_table.tsv
#   GO_enrichment_results.tsv            (optional - for enrichment-term match)
#
# Outputs:
#   topN_flank_hotspots_actinopterygii.tsv
#   topN_flank_hotspots_gnathostome.tsv
#   topN_flank_hotspots_summary.png
#
# Namespace note: dplyr verbs are fully qualified (dplyr::) so the script is
# robust even when AnnotationDbi / GO.db are attached and mask select()/etc.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# ----------------------------------------------------------------------------
# USER OPTIONS
# ----------------------------------------------------------------------------
# Drop genes whose symbol starts with "loc" or "si" (case-insensitive):
#   loc... = NCBI uncharacterised genes
#   si...  = zebrafish clone-based provisional names (e.g. si:dkey-...)
DROP_UNCHARACTERISED <- TRUE

# How many top genes to keep AFTER all filtering.
TOP_N <- 25

# Max GO terms to list per ontology, per gene.
MAX_TERMS_PER_ONT <- 3
# ----------------------------------------------------------------------------

# GO lookup needs the zebrafish OrgDb. Because unannotated genes are dropped,
# the output would be empty without it - so stop early with a clear message.
have_orgdb <- requireNamespace("org.Dr.eg.db", quietly = TRUE) &&
              requireNamespace("AnnotationDbi", quietly = TRUE) &&
              requireNamespace("GO.db", quietly = TRUE)
if (!have_orgdb) {
  stop("Packages org.Dr.eg.db / AnnotationDbi / GO.db are required for GO ",
       "annotation. Install with:\n",
       "  BiocManager::install(c('org.Dr.eg.db','AnnotationDbi','GO.db'))")
}

# ---- palette (colour-blind safe, Wong 2011) --------------------------------
ACTI  <- "#CC79A7"
GNATH <- "#009E73"

# ============================================================================
# 1. helper: load a CNE table and reassign genes via column Q
# ============================================================================


load_cne_table <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE,
                  guess_max = 100000, na = c("", "NA")) |>
    dplyr::mutate(
      gene_clean = tolower(trimws(sub("^gene-", "", as.character(geneId)))),
      gene_clean = dplyr::if_else(gene_clean == "", NA_character_, gene_clean)
    )
}

# ============================================================================
# 2. helper: rank genes, apply the loc/si filter, take the top N
# ============================================================================
top_hotspots <- function(cne_tbl, top_n = TOP_N,
                          drop_uncharacterised = DROP_UNCHARACTERISED) {
  out <- cne_tbl |>
    dplyr::filter(!is.na(gene_clean), gene_clean != "")

  if (drop_uncharacterised) {
    n_before <- dplyr::n_distinct(out$gene_clean)
    out <- out |>
      dplyr::filter(!stringr::str_detect(gene_clean,
                                         stringr::regex("^(loc|si)",
                                                        ignore_case = TRUE)))
    n_after <- dplyr::n_distinct(out$gene_clean)
    cat(sprintf("  loc/si filter: %d -> %d genes (dropped %d)\n",
                n_before, n_after, n_before - n_after))
  }

  out |>
    dplyr::count(gene_clean, name = "n_cnes") |>
    dplyr::arrange(dplyr::desc(n_cnes)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(rank = dplyr::row_number())
}

# ============================================================================
# 3. GO terms per gene - ALL THREE ontologies (BP, MF, CC)
#    Returns a data frame: gene, go_bp, go_mf, go_cc, any_go (logical)
# ============================================================================
go_terms_all_ontologies <- function(genes, max_terms = MAX_TERMS_PER_ONT) {

  res <- data.frame(gene = genes, go_bp = NA_character_,
                    go_mf = NA_character_, go_cc = NA_character_,
                    stringsAsFactors = FALSE)

  # gene symbol -> Entrez ID
sym2eg <- tryCatch(
    AnnotationDbi::mapIds(org.Dr.eg.db::org.Dr.eg.db,
                          keys = genes, keytype = "ALIAS",
                          column = "ENTREZID", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(genes)), genes))

  for (i in seq_along(genes)) {
    eg <- sym2eg[genes[i]]          # single-bracket: returns NA if name absent
    eg <- unname(eg)
    if (is.null(eg) || length(eg) == 0 || is.na(eg)) next

    go_tab <- tryCatch(
      AnnotationDbi::select(org.Dr.eg.db::org.Dr.eg.db, keys = eg,
                            keytype = "ENTREZID",
                            columns = c("GO", "ONTOLOGY")),
      error = function(e) NULL)
    if (is.null(go_tab) || nrow(go_tab) == 0) next

    for (ont in c("BP", "MF", "CC")) {
      ids <- unique(go_tab$GO[go_tab$ONTOLOGY == ont & !is.na(go_tab$GO)])
      if (length(ids) == 0) next
      nm <- AnnotationDbi::select(GO.db::GO.db, keys = ids,
                                  columns = "TERM", keytype = "GOID")
      nm <- unique(nm$TERM[!is.na(nm$TERM)])
      if (length(nm) == 0) next
      col <- paste0("go_", tolower(ont))
      res[[col]][i] <- paste(utils::head(nm, max_terms), collapse = "; ")
    }
  }

  res$any_go <- !(is.na(res$go_bp) & is.na(res$go_mf) & is.na(res$go_cc))
  res
}

# ============================================================================
# 4. optional: enrichment-result terms (best effort; BP context only)
# ============================================================================
enrichment_terms_for_genes <- function(genes, go_file, cluster_label) {
  if (!file.exists(go_file))
    return(setNames(rep(NA_character_, length(genes)), genes))
  go <- readr::read_tsv(go_file, show_col_types = FALSE)
  gene_col <- intersect(c("geneID", "geneId", "core_enrichment", "Genes"),
                        colnames(go))
  if (length(gene_col) == 0)
    return(setNames(rep(NA_character_, length(genes)), genes))
  gene_col <- gene_col[1]
  go <- go |> dplyr::filter(Cluster == cluster_label)
  out <- setNames(rep(NA_character_, length(genes)), genes)
  for (i in seq_along(genes)) {
    hit <- go[stringr::str_detect(tolower(go[[gene_col]]),
                                  stringr::fixed(tolower(genes[i]))), ]
    if (nrow(hit) == 0) next
    hit <- hit |> dplyr::arrange(p.adjust) |> dplyr::slice_head(n = 3)
    out[genes[i]] <- paste(hit$Description, collapse = "; ")
  }
  out
}

# ============================================================================
# 5. build the summary for one dataset
# ============================================================================
build_summary <- function(cne_path, go_cluster_label, set_name) {
  cat(sprintf("\n=== %s ===\n", set_name))
  cne <- load_cne_table(cne_path)

  # rank more than TOP_N initially, because some genes will be dropped for
  # lacking GO annotation; re-trim to TOP_N afterwards.
  ranked <- top_hotspots(cne, top_n = TOP_N * 2)

  go_tab    <- go_terms_all_ontologies(ranked$gene_clean)
  go_enrich <- enrichment_terms_for_genes(ranked$gene_clean,
                                          "GO_enrichment_results_R1.tsv",
                                          go_cluster_label)

  res <- ranked |>
    dplyr::left_join(go_tab, by = c("gene_clean" = "gene")) |>
    dplyr::mutate(
      dataset         = set_name,
      go_terms_enrich = go_enrich[gene_clean]
    ) |>
    # DROP genes with no GO term in any ontology
    dplyr::filter(any_go) |>
    # re-rank and keep the final TOP_N
    dplyr::arrange(dplyr::desc(n_cnes)) |>
    dplyr::slice_head(n = TOP_N) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::select(dataset, rank, gene = gene_clean, n_cnes,
                  go_bp, go_mf, go_cc, go_terms_enrich)

  cat(sprintf("  kept %d annotated genes (top %d requested)\n",
              nrow(res), TOP_N))
  res
}

summary_acti  <- build_summary("actinopterygii_cne_final_table_R1_20260731.tsv",
                               "actinopterygii", "Actinopterygii-specific")
summary_gnath <- build_summary("gnathostomata_cne_final_table_R1_20260731.tsv",
                               "gnathostomata", "Gnathostome-conserved")

readr::write_tsv(
  summary_acti,
  sprintf("top%d_flank_hotspots_actinopterygii.tsv", TOP_N))
readr::write_tsv(
  summary_gnath,
  sprintf("top%d_flank_hotspots_gnathostome.tsv", TOP_N))
cat(sprintf("\nWrote top%d_flank_hotspots_*.tsv\n", TOP_N))

# ============================================================================
# 6. visualization - lollipop of CNE counts, faceted by dataset
#    Tip label = leading term, BP preferred, then MF, then CC.
# ============================================================================
plot_df <- dplyr::bind_rows(summary_acti, summary_gnath) |>
  dplyr::mutate(
    dataset = factor(dataset, levels = c("Actinopterygii-specific",
                                         "Gnathostome-conserved")),
    lead_term = dplyr::coalesce(
      stringr::str_split_i(go_bp, ";", 1),
      stringr::str_split_i(go_mf, ";", 1),
      stringr::str_split_i(go_cc, ";", 1)),
    go_short = stringr::str_trunc(lead_term, 40)
  )

plot_df <- plot_df |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(gene = factor(gene, levels = gene[order(n_cnes)])) |>
  dplyr::ungroup()

p <- ggplot(plot_df, aes(x = n_cnes, y = gene, colour = dataset, fill = dataset)) +
   geom_col(width = 0.7) +
  geom_text(aes(label = go_short), hjust = 0, nudge_x = 0.6,
            size = 2.5, colour = "grey25") +
  facet_wrap(~ dataset, scales = "free", ncol = 1) +
  scale_colour_manual(values = c("Actinopterygii-specific" = ACTI,
                                 "Gnathostome-conserved"   = GNATH),
                      guide = "none") +
  scale_fill_manual(values = c("Actinopterygii-specific" = ACTI,
                                 "Gnathostome-conserved"   = GNATH),
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.45))) +
  labs(
    # title    = "Top-25 CNE hotspot genes (flank-based assignment)",
    # subtitle = paste0("loc/si genes excluded; unannotated genes dropped. ",
    #                  "Label = leading GO term (BP > MF > CC)."),
    x = "Number of CNEs", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    panel.grid.major.y = element_blank(),
    strip.text    = element_text(face = "bold", size = 11),
    axis.text.y   = element_text(size = 8, face = "italic"))

ggsave(sprintf("top%d_flank_hotspots_summary_R1_20260731.pdf", TOP_N), p,
       width = 7, height = 8,  bg = "white")
cat(sprintf("\nWrote top%d_flank_hotspots_summary_R1_20260731.pdf\n", TOP_N))
