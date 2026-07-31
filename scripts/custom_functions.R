# Base R gained %||% in 4.4.0 and rlang exports it; define it only if neither
# is in scope so the scripts run on older R.
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# ---------- helper: stable IDs ----------
gr_id <- function(gr, include_strand = FALSE) {
  if (!include_strand) {
    paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  } else {
    paste0(
      as.character(seqnames(gr)),
      ":",
      start(gr),
      "-",
      end(gr),
      ":",
      as.character(strand(gr))
    )
  }
}


bin_granges <- function(gr_obj, bins, autosomes) {
  gr_obj <- gr_obj[seqnames(gr_obj) %in% autosomes]
  gr_obj <- keepSeqlevels(gr_obj, autosomes, pruning.mode = "coarse")

  hits <- countOverlaps(bins, gr_obj)

  df <- cbind(
    as.data.frame(bins)[, c("seqnames", "start", "end")],
    hit_count = as.numeric(hits)
  ) |>
    dplyr::rename(`Chromosome name` = seqnames) |>
    mutate(
      `Chromosome name` = factor(
        `Chromosome name`,
        levels = autosomes,
        ordered = TRUE
      )
    )

  df
}


plot_cne_width_powerlaw <- function(
  gr,
  pdf_file = NULL,
  main = "Power-law like distribution of CNE widths",
  xlab = "CNE width",
  ylab = "CDF",
  text_y = 0.7,
  line_col = 2,
  line_lwd = 2,
  pch = 1
) {
  widths <- width(gr)

  fits <- poweRlaw::displ$new(widths)

  xmin_est <- poweRlaw::estimate_xmin(fits)
  fits$setXmin(xmin_est)

  pars_est <- poweRlaw::estimate_pars(fits)

  txt <- paste(
    "xmin:",
    fits$xmin,
    "\n",
    "alpha:",
    format(fits$pars, digits = 3)
  )

  if (!is.null(pdf_file)) {
    pdf(pdf_file)
    on.exit(dev.off(), add = TRUE)
  }

  poweRlaw::plot(
    fits,
    xlab = xlab,
    ylab = ylab,
    main = main,
    pch = pch
  )

  poweRlaw::lines(fits, col = line_col, lwd = line_lwd)

  text(
    x = max(fits$dat),
    y = text_y,
    labels = txt,
    adj = c(1, 1)
  )

  invisible(list(
    fit = fits,
    xmin = fits$xmin,
    alpha = fits$pars
  ))
}

plot_gviz_zoom <- function(
  gr,
  chr,
  start,
  end,
  filepath = NULL,
  width = 8,
  height = 8,
  genome = NA,
  gene_id_col = "gene_id",
  gene_activity = NULL,
  active_col = "royalblue",
  inactive_col = "skyblue",
  unknown_col = "grey80",
  teleost_fill = NULL,
  verte_fill = NULL,
  atac_fill = NULL,
  background_title = "brown",
  collapseTranscripts = TRUE,
  shape = "arrow",
  transcriptAnnotation = "symbol",
  ...
) {
  stopifnot(is.character(chr), length(chr) == 1)
  stopifnot(
    is.numeric(start),
    is.numeric(end),
    length(start) == 1,
    length(end) == 1
  )
  if (start > end) {
    stop("`start` must be <= `end`.")
  }

  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(Gviz)
  })

  roi <- GRanges(chr, IRanges(start, end))

  subset_chr_overlap <- function(x, roi) {
    x_chr <- x[as.character(seqnames(x)) == as.character(seqnames(roi))[1]]
    subsetByOverlaps(x_chr, roi, ignore.strand = TRUE)
  }

  # --- subset inputs
  gene_of_interest <- subset_chr_overlap(gr$genome_genes_gr, roi)
  teleost_of_interest <- subset_chr_overlap(gr$actinopteriigy_cne_gr, roi)
  verte_of_interest <- subset_chr_overlap(gr$gnathostomata_cne_gr, roi)
  peaks_of_interest <- subset_chr_overlap(
    gr$atac_peaks_overlapping_actinopteriigy_cne_gr,
    roi
  )

  axisTrack <- GenomeAxisTrack()

  # gene_id grouping / labels
  if (!gene_id_col %in% names(mcols(gene_of_interest))) {
    warning(sprintf(
      "Column '%s' not found in genes metadata; falling back to showId=FALSE.",
      gene_id_col
    ))
    group_vec <- NULL
    show_id <- FALSE
  } else {
    group_vec <- mcols(gene_of_interest)[[gene_id_col]]
    show_id <- TRUE
  }

  # --- map gene_id -> active/inactive -> colors
  gene_fill <- NULL
  gene_col <- NULL

  if (!is.null(gene_activity) && !is.null(group_vec)) {
    if (!all(c("gene_id", "is_active") %in% colnames(gene_activity))) {
      stop("`gene_activity` must contain columns: gene_id, is_active")
    }

    ga <- gene_activity[, c("gene_id", "is_active")]
    ga <- ga[!is.na(ga$gene_id), ]
    ga <- ga[!duplicated(ga$gene_id), ]
    active_lookup <- setNames(as.logical(ga$is_active), ga$gene_id)

    is_active_vec <- unname(active_lookup[group_vec])
    gene_fill <- ifelse(
      is.na(is_active_vec),
      unknown_col,
      ifelse(is_active_vec, active_col, inactive_col)
    )
    gene_col <- gene_fill
  }

  atrack <- AnnotationTrack(
    gene_of_interest,
    chromosome = chr,
    name = "Genes",
    group = group_vec,
    showId = show_id,
    fill = gene_fill,
    col = gene_col,
    background.title = background_title
  )

  ctrack <- AnnotationTrack(
    teleost_of_interest,
    chromosome = chr,
    name = "Actinopteriigy\nCNEs",
    fill = teleost_fill,
    background.title = background_title
  )

  vtrack <- AnnotationTrack(
    verte_of_interest,
    chromosome = chr,
    name = "Gnathostomata\nCNEs",
    fill = verte_fill,
    background.title = background_title
  )

  ptrack <- AnnotationTrack(
    peaks_of_interest,
    chromosome = chr,
    name = "ATACseq\npeaks",
    fill = atac_fill,
    background.title = background_title
  )

  do_plot <- function() {
    plotTracks(
      list(axisTrack, atrack, ctrack, vtrack, ptrack),
      from = start,
      to = end,
      collapseTranscripts = collapseTranscripts,
      shape = shape,
      transcriptAnnotation = transcriptAnnotation,
      ...
    )
  }

  if (!is.null(filepath)) {
    pdf(filepath, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    do_plot()
  } else {
    do_plot()
  }

  invisible(list(
    axisTrack = axisTrack,
    atrack = atrack,
    ctrack = ctrack,
    vtrack = vtrack,
    ptrack = ptrack
  ))
}


# Tolerance-aware overlap: returns indices of `query` overlapping `subject`.
# `slack = 0` is strict overlap; > 0 allows that many bp of coordinate drift
# (useful for lifted-over coordinates). All call sites pass `slack` explicitly.
overlapping_idx <- function(query, subject, slack = 0L) {
  hits <- findOverlaps(query, subject, maxgap = slack, ignore.strand = TRUE)
  unique(queryHits(hits))
}


# ── Master per-universe analysis ─────────────────────────────────────────────

analyse_cne_universe <- function(
  anno_gr,
  label,
  atac_peaks_gr,
  enh_gr = NULL,
  yuesong_gr = NULL,
  fin_geneIds = NULL,
  slack = 0L,
  out_dir = '../output',
  upset_sets = NULL
) {
  # 1. Universe
  universe_gr <- unique(non_exon(anno_gr))
  n_uni <- length(universe_gr)

  # 2. Membership indices
  in_atac <- overlapping_idx(universe_gr, atac_peaks_gr, slack = slack)

  in_chan <- if (!is.null(enh_gr)) {
    overlapping_idx(universe_gr, enh_gr, slack = slack)
  } else {
    integer(0)
  }

  in_yue <- if (!is.null(yuesong_gr)) {
    overlapping_idx(universe_gr, yuesong_gr, slack = slack)
  } else {
    integer(0)
  }

  in_active <- which(
    !is.na(mcols(universe_gr)$is_active) &
      mcols(universe_gr)$is_active
  )

  in_fin <- if (length(fin_geneIds)) {
    which(universe_gr$geneId %in% fin_geneIds)
  } else {
    integer(0)
  }

  # 3. Membership table
  universe_col <- paste(label, "CNEs")

  mem <- data.frame(
    placeholder = TRUE,
    `ATAC peaks` = seq_len(n_uni) %in% in_atac,
    `Active genes nearby` = seq_len(n_uni) %in% in_active,
    check.names = FALSE
  )

  names(mem)[1] <- universe_col

  if (!is.null(enh_gr)) {
    mem$`Chan enhancers` <- seq_len(n_uni) %in% in_chan
  }

  if (!is.null(yuesong_gr)) {
    mem$`Yuesong CNEs` <- seq_len(n_uni) %in% in_yue
  }

  mem <- mem[, vapply(mem, any, logical(1)), drop = FALSE]

  # `upset_sets` limits which membership columns are drawn. The exported table
  # below still carries every column, so restricting the plot does not change
  # what downstream scripts and the Shiny app read.
  mem_plot <- if (is.null(upset_sets)) {
    mem
  } else {
    wanted <- c(universe_col, upset_sets)
    missing <- setdiff(wanted, names(mem))
    if (length(missing)) {
      warning(
        "upset_sets not present for ",
        label,
        ": ",
        paste(missing, collapse = ", ")
      )
    }
    mem[, intersect(wanted, names(mem)), drop = FALSE]
  }

  # 4. UpSet
  comb <- make_comb_mat(as.matrix(mem_plot + 0L))

  pdf(
    file.path(out_dir, sprintf("upset_overlaps_%s.pdf", label)),
    width = 10,
    height = 6
  )

  draw(
    UpSet(
      comb,
      set_order = colnames(mem_plot),
      top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
      right_annotation = upset_right_annotation(comb, add_numbers = TRUE)
    )
  )

  dev.off()

  # 6. Export
  final_df <- as.data.frame(universe_gr)

  final_df$in_atac_peak <-
    as.integer(seq_len(n_uni) %in% in_atac)

  final_df$nearby_gene_active <-
    as.integer(seq_len(n_uni) %in% in_active)

  final_df$nearby_gene_fin_dev <-
    as.integer(seq_len(n_uni) %in% in_fin)

  if (!is.null(enh_gr)) {
    final_df$overlaps_chan_enhancer <-
      as.integer(seq_len(n_uni) %in% in_chan)
  }

  if (!is.null(yuesong_gr)) {
    final_df$overlaps_yuesong_cne <-
      as.integer(seq_len(n_uni) %in% in_yue)
  }

  stopifnot(
    sum(final_df$in_atac_peak) == length(in_atac),
    sum(final_df$nearby_gene_active) == length(in_active),
    sum(final_df$nearby_gene_fin_dev) == length(in_fin)
  )

  if (!is.null(enh_gr)) {
    stopifnot(
      sum(final_df$overlaps_chan_enhancer) == length(in_chan)
    )
  }

  if (!is.null(yuesong_gr)) {
    stopifnot(
      sum(final_df$overlaps_yuesong_cne) == length(in_yue)
    )
  }

  write.table(
    final_df,
    file.path(out_dir, sprintf("%s_cne_final_table.tsv", label)),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(
    list(
      universe = universe_gr,
      mem = mem,
      final = final_df
    )
  )
}


# Export the non-exonic annotated CNEs as TSVs (kept here because these are
# pure descriptive products of the preprocessing step).
non_exon <- function(gr) {
  gr[!str_detect(gr$annotation, "Exon") & !(is.na(gr$is_active))]
}


combine_great <- function(tbl, cluster_name) {
  as_tibble(tbl) %>%
    transmute(
      ID = id,
      Description = description,
      Cluster = cluster_name,
      FoldEnrichment = fold_enrichment,
      pvalue = p_value,
      p.adjust = p_adjust,
      p.adjust_hyper = p_adjust_hyper,
      GeneHits = observed_gene_hits,
      GeneSetSize = gene_set_size,
      RegionHits = observed_region_hits
    )
}

# Dual-test rule (McLean et al. 2010): require binomial AND hypergeometric
# adj p <= 0.05 for the headline call.
strict_call <- function(df) {
  df %>% filter(p.adjust <= 0.05 & p.adjust_hyper <= 0.05)
}


run_simplify <- function(
  tbl,
  label,
  out_dir = great_dir,
  padj_cut = 0.05,
  hyper_cut = 0.05
) {
  sig_ids <- as_tibble(tbl) %>%
    filter(p_adjust <= padj_cut, p_adjust_hyper <= hyper_cut) %>%
    pull(id) %>%
    unique()

  if (length(sig_ids) < 5) {
    message(
      "Too few sig terms for simplifyGO in ",
      label,
      " (n = ",
      length(sig_ids),
      ")"
    )
    return(invisible(NULL))
  }

  sim_mat <- GO_similarity(sig_ids, ont = "BP", db = "org.Dr.eg.db")

  pdf(
    file.path(out_dir, paste0("simplifyGO_", label, ".pdf")),
    width = 10,
    height = 8
  )
  go_clusters <- simplifyGO(sim_mat)
  dev.off()

  write.table(
    go_clusters,
    file.path(out_dir, paste0("simplifyGO_clusters_", label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  invisible(go_clusters)
}


extract_genes_for_terms <- function(res, tbl, term_descriptions) {
  term_ids <- tbl$id[tbl$description %in% term_descriptions]
  if (!length(term_ids)) {
    return(character(0))
  }
  unique(unlist(lapply(term_ids, function(tid) {
    assoc <- getRegionGeneAssociations(res, term_id = tid)
    unique(unlist(assoc$annotated_genes))
  })))
}


# Relabel RefSeq accessions -> chrN, but only where the relabel succeeds.
# (Previously this could silently NA out an already-renamed GRanges.)
relabel_seqlevels <- function(x, genome_df) {
  current <- seqlevels(x)
  new_levels <- genome_df$`Sequence name`[
    match(current, genome_df$`RefSeq seq accession`)
  ]
  if (all(is.na(new_levels))) {
    # already in chr* form — nothing to do
    return(x)
  }
  if (any(is.na(new_levels))) {
    keep <- !is.na(new_levels)
    seqlevels(x, pruning.mode = "coarse") <- current[keep]
    new_levels <- new_levels[keep]
  }
  seqlevels(x) <- new_levels
  seqlengths(x) <- genome_df$`Seq length`[
    match(seqlevels(x), genome_df$`Sequence name`)
  ]
  x
}


# ── CNE BED input ────────────────────────────────────────────────────────────

# phastCons BED output as written by the Snakemake workflow. BED is 0-based
# half-open, so `start` is shifted to 1-based inclusive on read. Seqnames are
# taken from column 4 (the phastCons element name carries the sequence
# accession) with any version suffix stripped, matching the naming used by the
# TxDb and the chrom-size tables.
read_cne_bed <- function(path, min_width = 25L) {
  bed <- readr::read_tsv(
    path,
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

  gr <- GRanges(
    seqnames = sub("\\.[0-9]+$", "", bed$cne_name),
    ranges = IRanges(bed$start + 1L, bed$end),
    phastcons = bed$phastcons_score
  )

  gr[width(gr) > min_width]
}


# ── Reviewer response: SYMBOL mapping rate ───────────────────────────────────

# Fraction of the custom GFF annotation that reaches org.*.eg.db, plus the
# fraction of GO gene-set members that receive a regulatory domain. The second
# rate is what actually sets the universe of the GREAT binomial test, so both
# are reported.
symbol_mapping_table <- function(gene_gr, extended_tss, gene_sets, orgdb) {
  sym_anno <- unique(as.character(mcols(gene_gr)$gene_id))
  sym_anno <- sym_anno[!is.na(sym_anno) & nzchar(sym_anno)]

  sym_org <- AnnotationDbi::keys(orgdb, "SYMBOL")
  sym_domain <- unique(as.character(mcols(extended_tss)$gene_id))
  sym_sets <- unique(unlist(gene_sets, use.names = FALSE))

  n_matched <- sum(sym_anno %in% sym_org)
  n_sets_dom <- sum(sym_sets %in% sym_domain)

  tibble::tibble(
    metric = c(
      "Genes in custom annotation",
      "Matched to org.Dr.eg.db SYMBOL",
      "Unmatched (dropped)",
      "Regulatory domains built",
      "GO gene-set members",
      "GO members with a regulatory domain"
    ),
    n = c(
      length(sym_anno),
      n_matched,
      length(sym_anno) - n_matched,
      length(extended_tss),
      length(sym_sets),
      n_sets_dom
    ),
    percent = c(
      100,
      round(100 * n_matched / length(sym_anno), 1),
      round(100 * (length(sym_anno) - n_matched) / length(sym_anno), 1),
      NA_real_,
      100,
      round(100 * n_sets_dom / length(sym_sets), 1)
    )
  )
}


# ── Reviewer response: domain-normalised CNE hotspots ────────────────────────

# Per-gene CNE density scored against the same GREAT regulatory domains used
# for the GO analysis. A CNE is counted once per overlapping domain, so the
# total number of CNE-domain assignments is the binomial n and each domain's
# success probability is its share of the total domain space. Domain length is
# therefore part of the null expectation rather than an unmodelled confounder,
# which is what raw nearest-gene counts cannot do.
hotspot_domains <- function(cne_gr, extended_tss, label = NA_character_) {
  shared <- intersect(seqlevels(cne_gr), seqlevels(extended_tss))
  if (!length(shared)) {
    stop("`cne_gr` and `extended_tss` share no seqlevels.")
  }

  cne <- keepSeqlevels(cne_gr, shared, pruning.mode = "coarse")
  dom <- keepSeqlevels(extended_tss, shared, pruning.mode = "coarse")

  hits <- findOverlaps(cne, dom, ignore.strand = TRUE)
  n_assign <- length(hits)
  total_space <- sum(as.numeric(width(dom)))

  obs_tbl <- tibble::tibble(idx = subjectHits(hits)) |>
    dplyr::count(idx, name = "observed")

  tibble::tibble(
    idx = seq_along(dom),
    gene = as.character(mcols(dom)$gene_id),
    domain_width = as.numeric(width(dom))
  ) |>
    dplyr::left_join(obs_tbl, by = "idx") |>
    dplyr::mutate(
      set = label,
      observed = dplyr::coalesce(observed, 0L),
      expected = n_assign * (domain_width / total_space),
      obs_exp = observed / expected,
      per_100kb = observed / (domain_width / 1e5),
      p_binom = stats::pbinom(
        observed - 1L,
        size = n_assign,
        prob = domain_width / total_space,
        lower.tail = FALSE
      ),
      fdr = stats::p.adjust(p_binom, method = "BH")
    ) |>
    dplyr::select(-idx) |>
    dplyr::relocate(set) |>
    dplyr::arrange(fdr, dplyr::desc(obs_exp))
}


# Teleost-specific duplicates arising from the 3R whole-genome duplication are
# conventionally named with `a`/`b` suffixes in ZFIN. Flagging them keeps
# ohnologue pairs (ryr1a/ryr1b) visible as separate loci in the hotspot table
# rather than being read as a single inflated locus.
flag_ohnologue_pairs <- function(hotspots) {
  stem <- sub("([ab])$", "", hotspots$gene)
  suffixed <- grepl("[ab]$", hotspots$gene)
  stem_counts <- table(stem[suffixed])

  hotspots |>
    dplyr::mutate(
      ohnologue_stem = ifelse(
        suffixed & stem %in% names(stem_counts)[stem_counts >= 2],
        stem,
        NA_character_
      ),
      is_ohnologue = !is.na(ohnologue_stem)
    )
}


# One panel of the hotspot figure. Kept as a function so the two CNE sets are
# drawn by identical code and differ only in fill colour; the ohnologue flag is
# carried on the bar outline rather than the fill so that a single colour legend
# is shared when the panels are combined.
plot_hotspot_panel <- function(
  hotspots,
  set_name,
  fill_colour,
  top_n = 20,
  title = NULL,
  label_size = 2.4
) {
  df <- hotspots |>
    dplyr::filter(set == set_name, expected > 0, observed > 0) |>
    dplyr::slice_min(fdr, n = top_n, with_ties = FALSE) |>
    dplyr::mutate(
      gene = forcats::fct_reorder(gene, obs_exp),
      duplicate_status = dplyr::if_else(
        is_ohnologue,
        "3R ohnologue pair member",
        "single-copy locus"
      )
    )

  ggplot2::ggplot(df, ggplot2::aes(gene, obs_exp)) +
    ggplot2::geom_col(
      ggplot2::aes(colour = duplicate_status),
      fill = fill_colour,
      alpha = 0.85,
      width = 0.75,
      linewidth = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(observed, " / ", round(domain_width / 1e3), " kb")
      ),
      hjust = -0.08,
      size = label_size
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_manual(
      values = c(
        "3R ohnologue pair member" = "black",
        "single-copy locus" = "grey55"
      )
    ) +
    ggplot2::expand_limits(y = max(df$obs_exp) * 1.35) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Observed / expected CNEs per regulatory domain",
      colour = NULL,
      title = title %||% set_name
    )
}


# ── Human-coordinate CNE locus views ─────────────────────────────────────────

# Locate the densest windows of a CNE set on one chromosome. `region` optionally
# restricts the search, which is how the telomeric windows are selected without
# hard-coding coordinates that would change with the element set.
find_cne_spikes <- function(
  gr,
  chrom,
  bin_size = 1e6,
  region = NULL,
  n = 3,
  chrom_length = NULL
) {
  x <- gr[as.character(seqnames(gr)) == chrom]
  if (!length(x)) {
    return(tibble::tibble())
  }

  sl <- seqlengths(gr)
  end_pos <- chrom_length
  if (is.null(end_pos) && chrom %in% names(sl)) {
    end_pos <- unname(sl[[chrom]])
  }
  if (is.null(end_pos) || is.na(end_pos)) {
    end_pos <- max(end(x))
  }

  starts <- seq(1, end_pos, by = bin_size)
  bins <- GRanges(
    chrom,
    IRanges(starts, pmin(starts + bin_size - 1L, end_pos))
  )

  if (!is.null(region)) {
    bins <- bins[start(bins) >= region[1] & end(bins) <= region[2]]
  }
  if (!length(bins)) {
    return(tibble::tibble())
  }

  tibble::tibble(
    chrom = chrom,
    start = start(bins),
    end = end(bins),
    n_cne = countOverlaps(bins, x, ignore.strand = TRUE)
  ) |>
    dplyr::arrange(dplyr::desc(n_cne)) |>
    dplyr::slice_head(n = n)
}


# Gviz panel for a human locus: binned density of each CNE set above the
# individual elements, with gene models underneath for context. Density is what
# makes a spike legible at megabase scale; the element track shows that the
# spike is many elements rather than a few wide ones.
plot_hsap_cne_locus <- function(
  cne_list,
  chrom,
  start,
  end,
  filepath,
  track_colours,
  txdb = NULL,
  orgdb = NULL,
  density_bin = 50e3,
  width = 10,
  height = 7
) {
  suppressPackageStartupMessages(library(Gviz))

  roi <- GRanges(chrom, IRanges(start, end))
  n_bins <- max(1L, ceiling((end - start + 1) / density_bin))
  bin_starts <- seq(start, end, length.out = n_bins + 1)
  bins <- GRanges(
    chrom,
    IRanges(
      round(head(bin_starts, -1)),
      round(tail(bin_starts, -1)) - 1L
    )
  )

  tracks <- list(GenomeAxisTrack())

  for (nm in names(cne_list)) {
    x <- cne_list[[nm]]
    x <- x[as.character(seqnames(x)) == chrom]

    tracks <- c(
      tracks,
      DataTrack(
        range = bins,
        data = countOverlaps(bins, x, ignore.strand = TRUE),
        type = "histogram",
        name = paste0(nm, "\nper ", round(density_bin / 1e3), " kb"),
        col.histogram = track_colours[[nm]],
        fill.histogram = track_colours[[nm]],
        background.title = "grey30"
      ),
      AnnotationTrack(
        subsetByOverlaps(x, roi, ignore.strand = TRUE),
        chromosome = chrom,
        name = nm,
        fill = track_colours[[nm]],
        col = track_colours[[nm]],
        stacking = "dense",
        background.title = "grey30"
      )
    )
  }

  if (!is.null(txdb)) {
    grt <- tryCatch(
      {
        g <- GeneRegionTrack(
          txdb,
          chromosome = chrom,
          start = start,
          end = end,
          name = "Genes",
          background.title = "grey30"
        )
        if (!is.null(orgdb) && length(gene(g))) {
          sym <- suppressMessages(AnnotationDbi::mapIds(
            orgdb,
            keys = sub("\\..*$", "", gene(g)),
            keytype = "ENTREZID",
            column = "SYMBOL"
          ))
          symbol(g) <- ifelse(is.na(sym), gene(g), sym)
        }
        g
      },
      error = function(e) {
        message("Gene track unavailable for ", chrom, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(grt)) {
      tracks <- c(tracks, grt)
    }
  }

  pdf(filepath, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  plotTracks(
    tracks,
    from = start,
    to = end,
    chromosome = chrom,
    transcriptAnnotation = "symbol",
    collapseTranscripts = "meta",
    shape = "arrow",
    main = sprintf(
      "%s:%.2f-%.2f Mb",
      chrom,
      start / 1e6,
      end / 1e6
    ),
    cex.main = 0.9
  )

  invisible(roi)
}


# ── Circos element density ───────────────────────────────────────────────────

# Reusable version of the element-density circos so that the zebrafish and
# human panels are drawn by identical code. `tracks` is an ordered, named list
# of GRanges already relabelled to the chr* scheme; the first element is the
# outermost track. `chrom_sizes` is a named numeric vector giving the plotting
# order.
plot_cne_circos <- function(
  tracks,
  chrom_sizes,
  filepath,
  track_colours,
  legend_labels = names(tracks),
  bin_size = 2e6,
  track_height = 0.12,
  axis_step = NULL,
  label_gap = 0.30,
  width = 9,
  height = 9,
  legend_title = NULL
) {
  stopifnot(length(tracks) == length(track_colours))

  chroms <- names(chrom_sizes)
  bins <- tileGenome(
    chrom_sizes,
    tilewidth = bin_size,
    cut.last.tile.in.chrom = TRUE
  )

  # A track with no elements on some chromosome would otherwise drop that
  # seqlevel and fail the keepSeqlevels() call inside bin_granges().
  tracks <- lapply(tracks, function(x) {
    x <- x[as.character(seqnames(x)) %in% chroms]
    seqlevels(x) <- chroms
    x
  })

  binned <- lapply(tracks, bin_granges, bins = bins, autosomes = chroms)

  # Axis ticks every `axis_step` bp, chosen from the longest chromosome so the
  # same function works for a 78 Mb zebrafish chromosome and a 248 Mb human one.
  if (is.null(axis_step)) {
    axis_step <- if (max(chrom_sizes) > 150e6) 50e6 else 20e6
  }
  ticks <- seq(0, max(chrom_sizes), by = axis_step)
  tick_labels <- paste0(ticks / 1e6, " Mb")

  xlim <- cbind(start = rep(1, length(chroms)), end = as.numeric(chrom_sizes))

  pdf(filepath, width = width, height = height)
  on.exit(
    {
      circos.clear()
      dev.off()
    },
    add = TRUE
  )

  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = c(rep(1, length(chroms) - 1), 4),
    points.overflow.warning = FALSE
  )
  circos.initialize(factors = chroms, xlim = xlim)

  for (i in seq_along(binned)) {
    df <- binned[[i]]
    col_i <- track_colours[[i]]
    draw_axis <- i == 1L

    circos.trackPlotRegion(
      factors = df$`Chromosome name`,
      x = df$start,
      y = df$hit_count,
      track.height = track_height,
      panel.fun = function(x, y) {
        circos.lines(x, y, col = col_i, area = TRUE, lwd = 1, type = "s")
        circos.segments(
          x0 = x,
          y0 = 0,
          x1 = x,
          y1 = y,
          col = adjustcolor("grey70", alpha.f = 0.3)
        )
        if (draw_axis) {
          circos.xaxis(
            labels.facing = "clockwise",
            labels.niceFacing = TRUE,
            major.at = ticks,
            labels = tick_labels,
            labels.cex = 0.5
          )
        }
      }
    )
  }

  # Chromosome labels placed relative to the first track's own y-range, so the
  # offset does not have to be retuned when the count scale changes.
  # Place chromosome labels inside the innermost track
  innermost_track <- length(binned)

  for (chr in chroms) {
    ylim_inner <- get.cell.meta.data(
      "ylim",
      sector.index = chr,
      track.index = innermost_track
    )

    circos.text(
      sector.index = chr,
      track.index = innermost_track,
      x = as.numeric(chrom_sizes[[chr]]) / 2,
      y = ylim_inner[1] - diff(ylim_inner) * label_gap, # Offsets inwards past the bottom of the inner track
      labels = chr,
      facing = "clockwise", # Keeps short labels legible without bending distortion
      niceFacing = TRUE,
      cex = 0.55,
      adj = c(0, 0.5) # Smooth center alignment
    )
  }

  for (i in seq_along(binned)) {
    circos.yaxis(
      sector.index = chroms[1],
      track.index = i,
      labels.cex = 0.45
    )
  }

  par(fig = c(0, 1, 0, 1), new = TRUE)
  plot.new()
  legend(
    "center",
    legend = legend_labels,
    fill = unlist(track_colours),
    border = unlist(track_colours),
    bty = "n",
    cex = 0.9,
    title = legend_title
  )

  invisible(binned)
}
