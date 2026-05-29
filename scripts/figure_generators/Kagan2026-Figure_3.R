# ============================================================================
# Three-row GO comparison figure (actinopterygii-specific vs gnathostome)
# rGREAT-compatible version. Replaces the original two-panel design with a
# domain-split layout:
#
#   Row 1 : Developmental DAG (full width)
#   Row 2 : Signaling DAG  |  Cellular / metabolic DAG
#   Row 3 : GO-slim diverging bar chart (full width, unsplit summary)
#
# All three GO DAGs share fill / size / shape scales (single collected
# legend) so they are directly comparable. Node fill = log2(fold_acti /
# fold_gnath), size = gene count, shape = enriched vs connecting ancestor.
#
# Domain assignment: each enriched term is classified into exactly one
# domain in the order Signaling -> Cellular/metabolic -> Developmental.
# This priority keeps the narrow, hand-curated categories (the first two)
# from being absorbed by the broad "developmental process" root.
#
# Label visibility is controlled per-domain via LABEL_CONTRAST (named list).
#
# Input (working directory):
#   GO_enrichment_results.tsv   rGREAT output
#
# Output:
#   go_two_panel_figure.pdf
#   go_two_panel_panelA_nodes.tsv   (includes a domain column)
# ============================================================================

suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(GSEABase)
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
})

# ---- user-tunable thresholds ----------------------------------------------
PADJ_BINOM <- 0.01    # rGREAT binomial adjusted p
PADJ_HYPER <- 0.05    # rGREAT hypergeometric adjusted p
FOLD_MIN   <- 2       # fold enrichment floor

# Per-domain label cutoff on |log2(fold_acti / fold_gnath)|. A term gets a
# text label only if it is enriched AND its absolute contrast meets this
# threshold. Bump the value for a domain whose sub-panel looks cluttered.
# Names MUST match the domain names in category_anchors below.
LABEL_CONTRAST <- list(
  "Developmental"        = 1.3,
  "Signaling"            = 1,
  "Cellular / metabolic" = 1
)
INCLUDE_SIGNALING <- FALSE   # FALSE drops the Signaling sub-panel

# ---- palette --------------------------------------------------------------
COL_GNATH <- "#009E73"
COL_MID   <- "grey88"
COL_ACTI  <- "#CC79A7"

# ---- domain definitions (anchors used to build BP sub-trees) --------------
# Lifted from the go_subtrees script. Signaling and cellular/metabolic are
# narrowly defined; developmental uses the broad root so it picks up
# everything else.
category_anchors <- list(
  "Signaling" = c(
    "GO:0008543",  # fibroblast growth factor receptor signaling pathway
    "GO:0016055",  # Wnt signaling pathway
    "GO:0007219",  # Notch signaling pathway
    "GO:0007224",  # smoothened signaling pathway
    "GO:0032925",  # regulation of activin receptor signaling pathway
    "GO:0034199",  # activation of cAMP-mediated signaling
    "GO:0007188",  # adenylate cyclase-modulating GPCR signaling
    "GO:0007264",  # small GTPase mediated signal transduction
    "GO:0048545"   # response to steroid hormone
  ),
  "Cellular / metabolic" = c(
    "GO:0006260",  # DNA replication
    "GO:0006325",  # chromatin organization
    "GO:0007049",  # cell cycle
    "GO:0006914",  # autophagy
    "GO:0006869",  # lipid transport
    "GO:0008203",  # cholesterol metabolic process
    "GO:0006811",  # monoatomic ion transport
    "GO:0001508"   # action potential
  ),
  "Developmental" = c(
    "GO:0032502"   # developmental process (catches all dev offspring)
  )
)

if (!INCLUDE_SIGNALING) category_anchors$Signaling <- NULL
# Priority order for mutually-exclusive assignment.
category_priority <- names(category_anchors)

# ============================================================================
# 1. Load, filter, build per-term contrast table
# ============================================================================
go_raw <- readr::read_tsv("GO_enrichment_results.tsv", show_col_types = FALSE)

cat(sprintf("Loaded %d rows from rGREAT result table.\n", nrow(go_raw)))

go <- go_raw |>
  dplyr::group_by(Cluster, ID) |>
  dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

go_sig <- go |>
  dplyr::filter(p.adjust       < PADJ_BINOM,
                p.adjust_hyper < PADJ_HYPER,
                FoldEnrichment >= FOLD_MIN)

cat(sprintf("After filter (p.adjust<%g, p.adjust_hyper<%g, fold>=%g): %d rows.\n",
            PADJ_BINOM, PADJ_HYPER, FOLD_MIN, nrow(go_sig)))

acti <- go_sig |>
  dplyr::filter(Cluster == "actinopterygii") |>
  dplyr::select(ID, Description, FoldEnrichment, p.adjust, GeneHits)
gnath <- go_sig |>
  dplyr::filter(Cluster == "gnathostomata") |>
  dplyr::select(ID, Description, FoldEnrichment, p.adjust, GeneHits)

terms <- dplyr::full_join(
  acti  |> dplyr::rename(fold_a = FoldEnrichment, padj_a = p.adjust,
                         count_a = GeneHits),
  gnath |> dplyr::rename(fold_g = FoldEnrichment, padj_g = p.adjust,
                         count_g = GeneHits),
  by = c("ID", "Description")
) |>
  dplyr::mutate(
    fold_a  = dplyr::coalesce(fold_a, 1),
    fold_g  = dplyr::coalesce(fold_g, 1),
    padj_a  = dplyr::coalesce(padj_a, 1),
    padj_g  = dplyr::coalesce(padj_g, 1),
    count_a = dplyr::coalesce(count_a, 0L),
    count_g = dplyr::coalesce(count_g, 0L),
    contrast   = log2(fold_a / fold_g),
    gene_count = pmax(count_a, count_g)
  )

cat(sprintf("Unique enriched GO BP terms: %d\n", nrow(terms)))

if (nrow(terms) == 0) {
  stop("No terms pass the significance filter. Loosen PADJ_BINOM / ",
       "PADJ_HYPER / FOLD_MIN at the top of the script.")
}

valid_keys   <- AnnotationDbi::keys(GO.db::GOTERM)
enriched_ids <- terms$ID[terms$ID %in% valid_keys]

# ============================================================================
# 2. Precompute domain sub-trees and classify each enriched term
# ============================================================================
category_subtrees <- lapply(category_anchors, function(anchors) {
  off <- unlist(lapply(anchors, function(a)
    tryCatch(AnnotationDbi::get(a, GO.db::GOBPOFFSPRING),
             error = function(e) character(0))))
  subtree <- unique(c(anchors, off))
  subtree[subtree %in% valid_keys]
})

# ancestor lookup (term itself + all BP ancestors), used both for
# classification and for finding connectors per domain.
anc_full <- AnnotationDbi::as.list(GO.db::GOBPANCESTOR)[enriched_ids]
anc_full <- lapply(anc_full, function(x) {
  if (is.null(x)) character(0) else x[!is.na(x) & x != "all"]
})
names(anc_full) <- enriched_ids
ancestry <- mapply(function(id, anc) c(id, anc),
                   enriched_ids, anc_full, SIMPLIFY = FALSE)

classify_one <- function(anc_set) {
  # mutually exclusive: first matching category wins
  for (cat in category_priority) {
    if (length(intersect(anc_set, category_subtrees[[cat]])) > 0) return(cat)
  }
  NA_character_
}
term_domain <- vapply(ancestry, classify_one, character(1))
names(term_domain) <- enriched_ids

cat("Domain assignment (mutually exclusive, by priority):\n")
print(table(term_domain, useNA = "ifany"))

# ============================================================================
# 3. Per-domain sub-panel builder
# ============================================================================
build_domain_panel <- function(domain_name) {

  domain_subtree <- category_subtrees[[domain_name]]
  in_domain      <- names(term_domain)[term_domain %in% domain_name]
  label_cut      <- LABEL_CONTRAST[[domain_name]]
  if (is.null(label_cut)) {
    stop(sprintf("No LABEL_CONTRAST entry for domain '%s'.", domain_name))
  }

  if (length(in_domain) < 2) {
    return(ggplot() +
             annotate("text", x = 0, y = 0,
                      label = sprintf("%s\n(< 2 enriched terms)", domain_name),
                      size = 3.5) +
             theme_void())
  }

  # connectors within this domain subtree, sitting above >=2 enriched terms
  anc_dom <- lapply(anc_full[in_domain], function(a) intersect(a, domain_subtree))
  cnt <- table(unlist(anc_dom))
  connectors <- names(cnt)[cnt >= 2]

  node_ids <- unique(c(in_domain, connectors))
  node_ids <- node_ids[node_ids %in% valid_keys & node_ids != "all"]

  # edges: parent->child, both endpoints in node_ids
  par <- AnnotationDbi::as.list(GO.db::GOBPPARENTS)[node_ids]
  ed <- do.call(rbind, lapply(names(par), function(ch) {
    pr <- intersect(par[[ch]], node_ids)
    if (length(pr) == 0) NULL else data.frame(from = pr, to = ch,
                                              stringsAsFactors = FALSE)
  }))
  ed <- unique(ed)

  # node attributes
  nm <- AnnotationDbi::select(GO.db::GO.db, keys = node_ids,
                              columns = "TERM", keytype = "GOID")
  nm <- nm[!duplicated(nm$GOID), ]
  go_name <- setNames(nm$TERM, nm$GOID)

  nd <- data.frame(name = node_ids, stringsAsFactors = FALSE) |>
    dplyr::mutate(term_name = unname(go_name[name]),
                  enriched  = name %in% in_domain) |>
    dplyr::left_join(terms |> dplyr::select(ID, contrast, gene_count),
                     by = c("name" = "ID")) |>
    dplyr::mutate(contrast   = dplyr::coalesce(contrast, 0),
                  gene_count = dplyr::coalesce(gene_count, 0L),
                  term_name  = dplyr::coalesce(term_name, name),
                  show_label = enriched & abs(contrast) >= label_cut,
                  short_lab  = ifelse(show_label,
                                      ifelse(nchar(term_name) > 32,
                                             paste0(substr(term_name,1,29),"..."),
                                             term_name),
                                      NA_character_))
  ed <- ed |> dplyr::filter(from %in% nd$name, to %in% nd$name)

  g <- tidygraph::tbl_graph(nodes = nd, edges = ed, directed = TRUE)

  ggraph(g, layout = "sugiyama") +
    geom_edge_link(edge_colour = "grey80", edge_width = 0.3,
                   arrow = arrow(length = unit(1, "mm"), type = "closed"),
                   end_cap = circle(1.4, "mm"),
                   start_cap = circle(1.4, "mm")) +
    geom_node_point(aes(size = gene_count, fill = contrast,
                        shape = enriched),
                    colour = "grey30", stroke = 0.3) +
    ggrepel::geom_text_repel(
      aes(x = x, y = y, label = short_lab),
      size = 2.3, segment.size = 0.2, segment.colour = "grey60",
      min.segment.length = 0, max.overlaps = 30, na.rm = TRUE) +
    scale_fill_gradient2(
      low = COL_GNATH, mid = COL_MID, high = COL_ACTI, midpoint = 0,
      name = "lineage contrast\n(log2 fold ratio:\nacti / gnath)",
      limits = range(terms$contrast, na.rm = TRUE)) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24),
                       labels = c(`TRUE` = "enriched term",
                                  `FALSE` = "connecting ancestor"),
                       name = NULL) +
    scale_size_continuous(range = c(1.5, 7), name = "gene count",
                          limits = c(0, max(terms$gene_count, na.rm = TRUE))) +
    labs(title = domain_name,
         subtitle = sprintf("%d enriched + %d connectors",
                            length(in_domain), length(connectors))) +
    theme_void(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}

# ---- export the union node table across the three sub-panels --------------
extract_nodes <- function(domain_name) {
  domain_subtree <- category_subtrees[[domain_name]]
  in_domain      <- names(term_domain)[term_domain %in% domain_name]
  label_cut      <- LABEL_CONTRAST[[domain_name]]
  if (length(in_domain) < 2) return(NULL)

  anc_dom <- lapply(anc_full[in_domain], function(a) intersect(a, domain_subtree))
  cnt <- table(unlist(anc_dom))
  connectors <- names(cnt)[cnt >= 2]

  node_ids <- unique(c(in_domain, connectors))
  node_ids <- node_ids[node_ids %in% valid_keys & node_ids != "all"]

  nm <- AnnotationDbi::select(GO.db::GO.db, keys = node_ids,
                              columns = "TERM", keytype = "GOID")
  nm <- nm[!duplicated(nm$GOID), ]
  go_name <- setNames(nm$TERM, nm$GOID)

  tibble::tibble(GO_ID = node_ids) |>
    dplyr::mutate(term_name = unname(go_name[GO_ID]),
                  node_type = ifelse(GO_ID %in% in_domain,
                                     "enriched term", "connecting ancestor"),
                  domain    = domain_name) |>
    dplyr::left_join(terms |> dplyr::select(ID, contrast, gene_count),
                     by = c("GO_ID" = "ID")) |>
    dplyr::mutate(contrast   = dplyr::coalesce(contrast, 0),
                  gene_count = dplyr::coalesce(gene_count, 0L),
                  term_name  = dplyr::coalesce(term_name, GO_ID),
                  labelled   = node_type == "enriched term" &
                               abs(contrast) >= label_cut)
}

panelA_table <- dplyr::bind_rows(lapply(category_priority, extract_nodes)) |>
  dplyr::arrange(domain, dplyr::desc(labelled), dplyr::desc(abs(contrast)))

readr::write_tsv(panelA_table, "go_two_panel_panelA_nodes.tsv")
cat(sprintf("Wrote go_two_panel_panelA_nodes.tsv (%d nodes across %d domains)\n",
            nrow(panelA_table), dplyr::n_distinct(panelA_table$domain)))

# ---- build the three sub-panels -------------------------------------------
# Each sub-panel is a free-standing ggraph; they get composed at the end
# (Developmental on row 1 full-width, Signaling | Cellular/metabolic on row 2,
# Panel B on row 3).
panelA_sub <- lapply(category_priority, build_domain_panel)
names(panelA_sub) <- category_priority

# ============================================================================
# 4. PANEL B - GO-slim diverging bar chart
# ============================================================================
slim_ids <- c(
  "GO:0008283", "GO:0030154", "GO:0007049", "GO:0006950", "GO:0007154",
  "GO:0008219", "GO:0040007", "GO:0048856", "GO:0009790", "GO:0007399",
  "GO:0048513", "GO:0030198", "GO:0006259", "GO:0016070", "GO:0006412",
  "GO:0006629", "GO:0005975", "GO:0006464", "GO:0006810", "GO:0007165",
  "GO:0065007", "GO:0002376", "GO:0032502", "GO:0040011", "GO:0000003"
)
slim_ids <- slim_ids[slim_ids %in% valid_keys]
slim_ok  <- length(slim_ids) > 0

if (slim_ok) {
  term2slim <- lapply(seq_along(enriched_ids), function(i) {
    intersect(ancestry[[i]], slim_ids)
  })
  names(term2slim) <- enriched_ids

  slim_map <- tibble::tibble(
    ID   = rep(enriched_ids, lengths(term2slim)),
    slim = unlist(term2slim)
  ) |>
    dplyr::left_join(terms |> dplyr::select(ID, contrast, padj_a, padj_g),
                     by = "ID") |>
    dplyr::mutate(slim_name = AnnotationDbi::Term(GO.db::GOTERM[slim]))

  n_min <- max(3, min(10, floor(nrow(terms) / 50)))

  slim_summary <- slim_map |>
    dplyr::group_by(slim, slim_name) |>
    dplyr::summarise(
      n_terms      = dplyr::n(),
      mean_contr   = mean(contrast, na.rm = TRUE),
      n_acti_sig   = sum(padj_a < 0.05, na.rm = TRUE),
      n_gnath_sig  = sum(padj_g < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(n_terms >= n_min) |>
    dplyr::arrange(mean_contr) |>
    dplyr::mutate(slim_name = factor(slim_name, levels = slim_name))

  cat(sprintf("Panel B: %d GO-slim categories (n_terms >= %d).\n",
              nrow(slim_summary), n_min))

  panelB <- ggplot(slim_summary,
                   aes(x = mean_contr, y = slim_name, fill = mean_contr)) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_text(aes(label = paste0("n=", n_terms),
                  hjust = ifelse(mean_contr >= 0, -0.2, 1.2)),
              size = 2.6, colour = "grey25") +
    scale_fill_gradient2(low = COL_GNATH, mid = COL_MID, high = COL_ACTI,
                         midpoint = 0, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(title = "GO-slim summary of the lineage contrast",
         subtitle = paste0("Enriched terms collapsed to generic GO slim; ",
                           "bar = mean log2(fold_acti / fold_gnath)."),
         x = "<- stronger in gnathostome      |      stronger in actinopterygii ->",
         y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          panel.grid.major.y = element_blank(),
          axis.text.y   = element_text(size = 9))
} else {
  panelB <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "GO-slim collection unavailable") +
    theme_void()
}

# ============================================================================
# 5. Compose
# Row 1: Developmental (full width)
# Row 2: Signaling | Cellular / metabolic
# Row 3: GO-slim summary (Panel B, full width)
# ============================================================================
# Parens around (Signaling | Cellular/metabolic) are required: in R, `/`
# binds tighter than `|`, so without them patchwork would parse the
# expression as (Dev/Sig) | (Cell/PanelB).


row2 <- if ("Signaling" %in% names(panelA_sub)) {
  panelA_sub[["Signaling"]] | panelA_sub[["Cellular / metabolic"]]
} else {
  panelA_sub[["Cellular / metabolic"]]
}

fig <- panelA_sub[["Developmental"]] / row2 / panelB +
  patchwork::plot_layout(heights = c(1.8, 1.2, 1), guides = "collect")

# Taller canvas because the figure is now three rows. Tune to taste; if the
# Developmental panel still looks cramped, bump its row in `heights` above.
ggsave("go_two_panel_figure.pdf", fig,
       width = 12, height = 16, bg = "white")

cat("Wrote go_two_panel_figure.pdf\n")