# ============================================================================
# Process-specific GO sub-tree panels (rGREAT-compatible)
# Zooms into biological branches and tests whether the modulator-vs-core
# (depth-vs-lineage-contrast) pattern holds within each.
#
# Input:  GO_enrichment_results.tsv  (rGREAT output)
# Output: go_subtree_<process>.pdf  x N, plus go_subtrees_combined.png
#
# Changes vs the original (clusterProfiler-style) script:
#   * Cluster label spelling: "actinopteriigy" -> "actinopterygii"
#   * Column rename:           Count -> GeneHits
#   * Adds rGREAT-appropriate significance filtering (binomial AND
#     hypergeometric q-values plus a fold-enrichment floor); without it
#     "enriched_ids" would include every term in the BP universe.
#   * Legend label corrected to log2 fold ratio (matches the formula).
# ============================================================================

suppressPackageStartupMessages({
  library(GO.db); library(AnnotationDbi)
  library(igraph); library(ggraph); library(tidygraph)
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

COL_GNATH <- "#009E73"; COL_MID <- "grey88"; COL_ACTI <- "#CC79A7"

# ---- user-tunable thresholds ----------------------------------------------
# rGREAT returns the full BP universe; must filter before treating the
# result table as "enriched terms". Same defaults as the two-panel script.
PADJ_BINOM <- 0.01    # adjusted binomial p-value (rGREAT's `p.adjust`)
PADJ_HYPER <- 0.05    # adjusted hypergeometric p-value (`p.adjust_hyper`)
FOLD_MIN   <- 2       # fold enrichment floor

# ---- per-term contrast table (same as the two-panel script) ---------------
go_raw <- readr::read_tsv("GO_enrichment_results_R1_20260731.tsv", show_col_types = FALSE)

# collapse duplicate (Cluster, ID) rows to the best p.adjust, then filter
go <- go_raw |>
  dplyr::group_by(Cluster, ID) |>
  dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::filter(p.adjust       < PADJ_BINOM,
                p.adjust_hyper < PADJ_HYPER,
                FoldEnrichment >= FOLD_MIN)

cat(sprintf("After filter (p.adjust<%g, p.adjust_hyper<%g, fold>=%g): %d rows.\n",
            PADJ_BINOM, PADJ_HYPER, FOLD_MIN, nrow(go)))

mk <- function(cl) go |> dplyr::filter(Cluster == cl) |>
  dplyr::select(ID, Description, FoldEnrichment, p.adjust, GeneHits)

terms <- dplyr::full_join(
  mk("actinopterygii") |>
    dplyr::rename(fa = FoldEnrichment, pa = p.adjust, ca = GeneHits),
  mk("gnathostomata") |>
    dplyr::rename(fg = FoldEnrichment, pg = p.adjust, cg = GeneHits),
  by = c("ID", "Description")) |>
  dplyr::mutate(
    fa = dplyr::coalesce(fa, 1), fg = dplyr::coalesce(fg, 1),
    pa = dplyr::coalesce(pa, 1), pg = dplyr::coalesce(pg, 1),
    ca = dplyr::coalesce(ca, 0L), cg = dplyr::coalesce(cg, 0L),
    contrast   = log2(fa / fg),
    gene_count = pmax(ca, cg))

cat(sprintf("Unique enriched GO BP terms after join: %d\n", nrow(terms)))

if (nrow(terms) == 0) {
  stop("No terms pass the significance filter. Loosen PADJ_BINOM / ",
       "PADJ_HYPER / FOLD_MIN at the top of the script.")
}

enriched_ids <- terms$ID
valid_keys   <- AnnotationDbi::keys(GO.db::GOTERM)
enriched_ids <- enriched_ids[enriched_ids %in% valid_keys]

# ---- the processes, each defined by one or more anchor GO terms -----------
processes <- list(
  "Nervous system" = c(
    "GO:0060322", # head development
    # "GO:0021551", # central nervous system morphogenesis
    "GO:0021536", # diencephalon development
    "GO:0021953", # central nervous system neuron differentiation
    "GO:0021761", # limbic system development
    # "GO:0048665", # neuron fate specification
    # "GO:0021675", # nerve development
    # "GO:0031103", # axon regeneration
    # "GO:0048678", # response to axon injury
    # "GO:0001764", # neuron migration
    # "GO:0045685", # regulation of glial cell differentiation
    # "GO:0022010", # central nervous system myelination
    # "GO:0008088",	# axo-dendritic transport
    # "GO:0032291"  # axon ensheathment in central nervous system
    "GO:0030900" # forebrain development
    ),
  "Eye / visual system" = c(
  "GO:0048592",  # eye morphogenesis
  "GO:0048593",  # camera-type eye morphogenesis
  "GO:0060041",  # retina development in camera-type eye
  "GO:0060042",  # retina morphogenesis in camera-type eye
  "GO:0002089",  # lens morphogenesis in camera-type eye
  "GO:0031076",  # embryonic camera-type eye development
  "GO:0002072",  # optic cup morphogenesis involved in camera-type eye development
  "GO:0061386",  # closure of optic fissure
  "GO:0001754"   # eye photoreceptor cell differentiation
),
  "Inner ear / otolith development" = c(
    "GO:0032474",	# otolith morphogenesis
    "GO:0032475",	# otolith formation
    "GO:0060872",	# semicircular canal development
    "GO:0042471",	# ear morphogenesis
    "GO:0042472", #	inner ear morphogenesis
    "GO:0043583",  # ear development
    "GO:0048839",	# inner ear development
    "GO:0048840",	# otolith development
    "GO:0071599", # otic vesicle development
    "GO:0071600", # otic vesicle morphogenesis
    "GO:0030916", # otic vesicle formation
    "GO:0043049", # otic placode formation
    "GO:0060788", # ectodermal placode formation
    "GO:0071696", # ectodermal placode development/morphogenesis
    "GO:0071697" # ectodermal placode development/morphogenesis
  ),
  "Muscle development"            = c(
    "GO:0061061", # muscle structure development
    "GO:0060537", # muscle tissue development
    "GO:0055013" # cardiac muscle cell development
    # "GO:0006942", # regulation of striated muscle contraction
    #"GO:0060297",	# regulation of sarcomere organization
    # "GO:0014904",	# myotube cell development
    # "GO:0048339",	# paraxial mesoderm development
    # "GO:0001756"  # somitogenesis
    # "GO:0007379",  # segment specification
    # "GO:0061053"  # somite development    
  ),
    "Limb / fin / skeletal development"        = c(
    "GO:0060173",  # limb development
    "GO:0048706",  # embryonic skeletal system development (padj 7.7e-38 gnath)
    "GO:0001501",  # skeletal system development
    "GO:0032330",  # regulation of chondrocyte differentiation (FE 6.68 gnath)
    "GO:0030282",   # bone mineralization
    "GO:0048736", # appendage development 
    "GO:0033334", # fin morphogenesis
    "GO:0033333" # fin development
  ),
  "Kidney / pronephros" = c(
    "GO:0072001", # renal system development
    "GO:0048793", # pronephros development
    "GO:0039019", # pronephric nephron development
    "GO:0039020", # pronephric nephron tubule development
    # "GO:0003014", # renal system process
    "GO:0072088",	# nephron epithelium morphogenesis
    "GO:0072028",	# nephron morphogenesis
    "GO:0072080", #	nephron tubule development
    "GO:0072102"	# glomerulus morphogenesis
    ),
  "Neural crest / head"   = c(
    "GO:0014032", # neural crest cell development
    "GO:0014033", # neural crest cell differentiation
    "GO:0001755", # neural crest cell migration
    "GO:1904888",  # cranial skeletal system development
    "GO:0060323",	# head morphogenesis
    "GO:0060324"	# face development
    ),
  "Pigmentation"                  = c(
    "GO:0043473", # pigmentation
    "GO:0048066", # developmental pigmentation
    "GO:0030318", # melanocyte differentiation
    "GO:0050931"  # pigment cell differentiation
  )

)

# ---- builder: induced subgraph panel for one process ----------------------
build_panel <- function(proc_name, anchors, label_contrast = 2) {

  # sub-tree = anchors + all their BP offspring
  off <- unlist(lapply(anchors, function(a)
    tryCatch(AnnotationDbi::get(a, GO.db::GOBPOFFSPRING),
             error = function(e) character(0))))
  subtree <- unique(c(anchors, off))
  subtree <- subtree[subtree %in% valid_keys]

  # enriched terms that fall in this sub-tree
  in_sub <- intersect(enriched_ids, subtree)
  if (length(in_sub) < 2) {
    message(sprintf("[%s] only %d enriched term(s) - skipped.",
                    proc_name, length(in_sub)))
    return(NULL)
  }

  # connectors: ancestors (within the sub-tree) shared by >=2 enriched terms,
  # so the induced graph is connected
  anc <- AnnotationDbi::as.list(GO.db::GOBPANCESTOR)[in_sub]
  anc <- lapply(anc, function(x) intersect(x, subtree))
  cnt <- table(unlist(anc))
  connectors <- names(cnt)[cnt >= 2]
  node_ids <- unique(c(in_sub, connectors, anchors))
  node_ids <- node_ids[node_ids %in% valid_keys]

  # edges from GOBPPARENTS, both ends in node_ids
  par <- AnnotationDbi::as.list(GO.db::GOBPPARENTS)[node_ids]
  ed <- do.call(rbind, lapply(names(par), function(ch) {
    pr <- intersect(par[[ch]], node_ids)
    if (length(pr) == 0) NULL else data.frame(from = pr, to = ch)
  }))
  ed <- unique(ed)

  # node attributes
  nm <- AnnotationDbi::select(GO.db::GO.db, keys = node_ids,
                              columns = "TERM", keytype = "GOID")
  nm <- nm[!duplicated(nm$GOID), ]
  go_name <- setNames(nm$TERM, nm$GOID)

  nd <- data.frame(name = node_ids) |>
    dplyr::mutate(term_name = unname(go_name[name]),
                  enriched  = name %in% in_sub) |>
    dplyr::left_join(terms |> dplyr::select(ID, contrast, gene_count),
                     by = c("name" = "ID")) |>
    dplyr::mutate(contrast   = dplyr::coalesce(contrast, 0),
                  gene_count = dplyr::coalesce(gene_count, 0L),
                  term_name  = dplyr::coalesce(term_name, name),
                  # show_label = enriched & abs(contrast) >= label_contrast,
                  show_label = TRUE,
                  short_lab  = ifelse(show_label,
                                      ifelse(nchar(term_name) > 32,
                                             paste0(substr(term_name,1,29),"..."),
                                             term_name), NA))
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
    geom_text_repel(aes(x = x, y = y, label = short_lab),
                    size = 2.3, segment.size = 0.2, max.overlaps = 30,
                    min.segment.length = 0, na.rm = TRUE) +
    scale_fill_gradient2(low = COL_GNATH, mid = COL_MID, high = COL_ACTI,
                         midpoint = 0,
                         name = "contrast\n(log2 fold ratio:\nacti / gnath)") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), guide = "none") +
    scale_size_continuous(range = c(1.5, 7), name = "gene count") +
    labs(title = proc_name,
         subtitle = sprintf("%d enriched terms in sub-tree", length(in_sub))) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}

# ---- build all panels, save individually + combined -----------------------
panels <- list()
for (pn in names(processes)) {
  p <- build_panel(pn, processes[[pn]])
  if (is.null(p)) next
  panels[[pn]] <- p
  fn <- sprintf("go_subtree_%s.pdf",
                gsub("[^A-Za-z0-9]+", "_", tolower(pn)))
  ggsave(fn, p, width = 8, height = 6.5, bg = "white")
  cat("Wrote", fn, "\n")
}

if (length(panels) > 0) {
  combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "GO sub-tree zoom: does the modulator-vs-core pattern hold per process?",
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave("2026Kagan-Supplementary_Figure9_R1.pdf", combined,
         width = 10, height = 15, bg = "white")
  cat("Wrote 2026Kagan-Supplementary_Figure9_R1.pdf\n")
}
