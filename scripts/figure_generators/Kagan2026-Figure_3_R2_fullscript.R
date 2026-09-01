# ==========================================================================
# Kagan2026 Figure 3 - revised version
#
# Changes relative to the original Figure 3:
#   * Panel A redesigned to improve readability.
#   * Panel B retains the original DAG-style layout, but labels are no longer
#     abbreviated with "..."; full GO-term labels are wrapped instead.
#   * Panel C unchanged in concept.
#
# Input:
#   GO_enrichment_results.tsv
#
# Output:
#   Kagan2026-Fig3_R2.pdf
# ===========================================================================

suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

PADJ_BINOM <- 0.01
PADJ_HYPER <- 0.05
FOLD_MIN   <- 2
INPUT_FILE <- "GO_enrichment_results.tsv"

COL_GNATH <- "#009E73"
COL_MID   <- "grey88"
COL_ACTI  <- "#CC79A7"
COL_LINE  <- "grey72"

wrap_label <- function(x, width = 24) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"),
         character(1))
}

# --------------------------------------------------------------------------
# Load data and derive enriched-term summaries
# --------------------------------------------------------------------------
go_raw <- readr::read_tsv(INPUT_FILE, show_col_types = FALSE) |>
  dplyr::group_by(Cluster, ID) |>
  dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    significant = p.adjust < PADJ_BINOM &
      p.adjust_hyper < PADJ_HYPER &
      FoldEnrichment >= FOLD_MIN
  )

go_sig <- go_raw |> dplyr::filter(significant)

acti <- go_sig |>
  dplyr::filter(Cluster == "actinopterygii") |>
  dplyr::select(ID, Description, FoldEnrichment, p.adjust, GeneHits)

gnath <- go_sig |>
  dplyr::filter(Cluster == "gnathostomata") |>
  dplyr::select(ID, Description, FoldEnrichment, p.adjust, GeneHits)

terms <- dplyr::full_join(
  acti |>
    dplyr::rename(fold_a = FoldEnrichment,
                  padj_a = p.adjust,
                  count_a = GeneHits),
  gnath |>
    dplyr::rename(fold_g = FoldEnrichment,
                  padj_g = p.adjust,
                  count_g = GeneHits),
  by = c("ID", "Description")
) |>
  dplyr::mutate(
    fold_a = dplyr::coalesce(fold_a, 1),
    fold_g = dplyr::coalesce(fold_g, 1),
    padj_a = dplyr::coalesce(padj_a, 1),
    padj_g = dplyr::coalesce(padj_g, 1),
    count_a = dplyr::coalesce(count_a, 0L),
    count_g = dplyr::coalesce(count_g, 0L),
    contrast = log2(fold_a / fold_g),
    gene_count = pmax(count_a, count_g)
  )

valid_keys   <- AnnotationDbi::keys(GO.db::GOTERM)
enriched_ids <- terms$ID[terms$ID %in% valid_keys]

anc_full <- AnnotationDbi::as.list(GO.db::GOBPANCESTOR)[enriched_ids]
anc_full <- lapply(anc_full, function(x) {
  if (is.null(x)) character(0) else x[!is.na(x) & x != "all"]
})
names(anc_full) <- enriched_ids

# --------------------------------------------------------------------------
# Panel A - revised developmental-process presentation
# --------------------------------------------------------------------------
dev_terms <- c(
  "GO:0030900", # forebrain development
  "GO:0021555", # midbrain-hindbrain boundary morphogenesis
  "GO:0021984", # adenohypophysis development
  "GO:0060042", # retina morphogenesis in camera-type eye
  "GO:0071599", # otic vesicle development
  "GO:0060872", # semicircular canal development
  "GO:0007517", # muscle organ development
  "GO:0042692", # muscle cell differentiation
  "GO:0035118", # embryonic pectoral fin morphogenesis
  "GO:0001501", # skeletal system development
  "GO:0035775", # pronephric glomerulus morphogenesis
  "GO:0048793", # pronephros development
  "GO:0014032", # neural crest cell development
  "GO:0001755", # neural crest cell migration
  "GO:1904888", # cranial skeletal system development
  "GO:0048066", # developmental pigmentation
  "GO:0050931"  # pigment cell differentiation
)

panelA_long <- tidyr::expand_grid(
  ID = dev_terms,
  Cluster = c("gnathostomata", "actinopterygii")
) |>
  dplyr::left_join(
    go_raw |>
      dplyr::select(ID, Description, Cluster, FoldEnrichment,
                    p.adjust, p.adjust_hyper, GeneHits, significant),
    by = c("ID", "Cluster")
  ) |>
  dplyr::group_by(ID) |>
  dplyr::mutate(Description = dplyr::coalesce(
    Description,
    dplyr::first(Description[!is.na(Description)])
  )) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(FoldEnrichment)) |>
  dplyr::mutate(
    ID = factor(ID, levels = rev(dev_terms)),
    Cluster = factor(
      Cluster,
      levels = c("gnathostomata", "actinopterygii"),
      labels = c("Gnath-cons", "Acti-spec")
    ),
    label = wrap_label(Description, 38)
  )

panelA_wide <- panelA_long |>
  dplyr::select(ID, Cluster, FoldEnrichment, label) |>
  tidyr::pivot_wider(names_from = Cluster, values_from = FoldEnrichment)

panelA <- ggplot() +
  geom_segment(
    data = panelA_wide,
    aes(x = `Gnath-cons`, xend = `Acti-spec`, y = ID, yend = ID),
    colour = COL_LINE, linewidth = 0.55
  ) +
  geom_point(
    data = panelA_long,
    aes(x = FoldEnrichment, y = ID,
        colour = Cluster, fill = Cluster, shape = significant),
    size = 2.8, stroke = 0.8
  ) +
  scale_colour_manual(
    values = c("Gnath-cons" = COL_GNATH, "Acti-spec" = COL_ACTI),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Gnath-cons" = COL_GNATH, "Acti-spec" = COL_ACTI),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 1),
    labels = c(`TRUE` = "passes enrichment threshold",
               `FALSE` = "does not pass threshold"),
    name = NULL
  ) +
  scale_y_discrete(
    labels = setNames(panelA_wide$label, as.character(panelA_wide$ID))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Developmental processes",
    x = "Fold enrichment",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8.2, lineheight = 0.92),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

# --------------------------------------------------------------------------
# Panel B - signaling-process DAG, retaining original layout logic but using
#           wrapped full labels instead of "..." truncation
# --------------------------------------------------------------------------
signaling_anchors <- c(
  "GO:0008543", # fibroblast growth factor receptor signaling pathway
  "GO:0016055", # Wnt signaling pathway
  "GO:0060070", # canonical Wnt signaling pathway
  "GO:0030111", # regulation of Wnt signaling pathway
  "GO:0060828", # regulation of canonical Wnt signaling pathway
  "GO:0030178", # negative regulation of Wnt signaling pathway
  "GO:0090090", # negative regulation of canonical Wnt signaling pathway
  "GO:0007224", # smoothened signaling pathway
  "GO:0032924", # activin receptor signaling pathway
  "GO:0032925", # regulation of activin receptor signaling pathway
  "GO:0032926", # negative regulation of activin receptor signaling pathway
  "GO:0048384", # retinoic acid receptor signaling pathway
  "GO:0048385"  # regulation of retinoic acid receptor signaling pathway
)

sig_offspring <- unlist(lapply(signaling_anchors, function(a) {
  tryCatch(AnnotationDbi::get(a, GO.db::GOBPOFFSPRING),
           error = function(e) character(0))
}))
sig_subtree <- unique(c(signaling_anchors, sig_offspring))
sig_subtree <- sig_subtree[sig_subtree %in% valid_keys]

sig_enriched <- enriched_ids[enriched_ids %in% sig_subtree]

sig_anc <- lapply(anc_full[sig_enriched], function(a) intersect(a, sig_subtree))
sig_counts <- table(unlist(sig_anc))
sig_connectors <- names(sig_counts)[sig_counts >= 2]

sig_nodes <- unique(c(sig_enriched, sig_connectors))
sig_nodes <- sig_nodes[sig_nodes %in% valid_keys & sig_nodes != "all"]

sig_par <- AnnotationDbi::as.list(GO.db::GOBPPARENTS)[sig_nodes]
sig_edges <- do.call(rbind, lapply(names(sig_par), function(ch) {
  pr <- intersect(sig_par[[ch]], sig_nodes)
  if (length(pr) == 0) NULL else data.frame(
    from = pr, to = ch, stringsAsFactors = FALSE
  )
}))
sig_edges <- unique(sig_edges)

sig_names <- AnnotationDbi::select(
  GO.db::GO.db, keys = sig_nodes, columns = "TERM", keytype = "GOID"
)
sig_names <- sig_names[!duplicated(sig_names$GOID), ]
sig_name_map <- setNames(sig_names$TERM, sig_names$GOID)

sig_nd <- data.frame(name = sig_nodes, stringsAsFactors = FALSE) |>
  dplyr::mutate(
    term_name = unname(sig_name_map[name]),
    enriched = name %in% sig_enriched
  ) |>
  dplyr::left_join(
    terms |> dplyr::select(ID, contrast, gene_count),
    by = c("name" = "ID")
  ) |>
  dplyr::mutate(
    contrast = dplyr::coalesce(contrast, 0),
    gene_count = dplyr::coalesce(gene_count, 0L),
    term_name = dplyr::coalesce(term_name, name),
    # full wrapped labels; no "..." abbreviations
    term_label = ifelse(enriched, wrap_label(term_name, 24), NA_character_)
  )

sig_edges <- sig_edges |>
  dplyr::filter(from %in% sig_nd$name, to %in% sig_nd$name)

sig_graph <- tidygraph::tbl_graph(
  nodes = sig_nd, edges = sig_edges, directed = TRUE
)

panelB <- ggraph(sig_graph, layout = "sugiyama") +
  geom_edge_link(
    edge_colour = "grey80", edge_width = 0.3,
    arrow = arrow(length = unit(1, "mm"), type = "closed"),
    end_cap = circle(1.4, "mm"),
    start_cap = circle(1.4, "mm")
  ) +
  geom_node_point(
    aes(size = gene_count, fill = contrast, shape = enriched),
    colour = "grey30", stroke = 0.3
  ) +
  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = term_label),
    size = 2.3,
    lineheight = 0.92,
    segment.size = 0.2,
    segment.colour = "grey60",
    box.padding = 0.18,
    point.padding = 0.08,
    min.segment.length = 0,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = COL_GNATH, mid = COL_MID, high = COL_ACTI, midpoint = 0,
    name = "lineage contrast\n(log2 fold ratio:\nacti / gnath)",
    limits = range(terms$contrast, na.rm = TRUE)
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 24),
    labels = c(`TRUE` = "enriched term", `FALSE` = "connecting ancestor"),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.5, 7), name = "gene count",
    limits = c(0, max(terms$gene_count, na.rm = TRUE))
  ) +
  labs(
    title = "Signaling process",
    subtitle = sprintf("%d enriched + %d connectors",
                       length(sig_enriched), length(sig_connectors))
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey30"),
    legend.position = "right"
  )

# --------------------------------------------------------------------------
# Panel C - GO-slim summary
# --------------------------------------------------------------------------
slim_ids <- c(
  "GO:0008283", "GO:0030154", "GO:0007049", "GO:0006950", "GO:0007154",
  "GO:0008219", "GO:0040007", "GO:0048856", "GO:0009790", "GO:0007399",
  "GO:0048513", "GO:0030198", "GO:0006259", "GO:0016070", "GO:0006412",
  "GO:0006629", "GO:0005975", "GO:0006464", "GO:0006810", "GO:0007165",
  "GO:0065007", "GO:0002376", "GO:0032502", "GO:0040011", "GO:0000003"
)
slim_ids <- slim_ids[slim_ids %in% valid_keys]

term2slim <- lapply(seq_along(enriched_ids), function(i) {
  intersect(c(enriched_ids[i], anc_full[[i]]), slim_ids)
})
names(term2slim) <- enriched_ids

slim_map <- tibble::tibble(
  ID = rep(enriched_ids, lengths(term2slim)),
  slim = unlist(term2slim)
) |>
  dplyr::left_join(terms |> dplyr::select(ID, contrast), by = "ID") |>
  dplyr::mutate(slim_name = AnnotationDbi::Term(GO.db::GOTERM[slim]))

n_min <- max(3, min(10, floor(nrow(terms) / 50)))

slim_summary <- slim_map |>
  dplyr::group_by(slim, slim_name) |>
  dplyr::summarise(
    n_terms = dplyr::n(),
    mean_contr = mean(contrast, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n_terms >= n_min) |>
  dplyr::arrange(mean_contr) |>
  dplyr::mutate(slim_name = factor(slim_name, levels = slim_name))

panelC <- ggplot(
  slim_summary,
  aes(x = mean_contr, y = slim_name, fill = mean_contr)
) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_text(
    aes(label = paste0("n=", n_terms),
        hjust = ifelse(mean_contr >= 0, -0.2, 1.2)),
    size = 2.6, colour = "grey25"
  ) +
  scale_fill_gradient2(
    low = COL_GNATH, mid = COL_MID, high = COL_ACTI,
    midpoint = 0, guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    subtitle = paste0(
      "Enriched terms collapsed to generic GO slim; bar = mean ",
      "log2(fold_acti / fold_gnath)."
    ),
    x = "<- stronger in Gnath-cons   |   stronger in Acti-spec ->",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9)
  )

# --------------------------------------------------------------------------
# Compose and export
# --------------------------------------------------------------------------
fig <- panelA / panelB / panelC +
  patchwork::plot_layout(heights = c(2.25, 1.25, 1.25), guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(
  "Kagan2026-Fig3_R2.pdf",
  fig,
  width = 9.2,
  height = 14.8,
  units = "in",
  bg = "white"
)

cat("Wrote Kagan2026-Fig3_R2.pdf\n")
