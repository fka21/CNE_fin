#' Calculate A/T Nucleotide Frequencies in Binned Windows
#'
#' Calculates the frequency of A or T nucleotides within bins of a user-specified width across each sequence in a \code{DNAStringSet} object.
#'
#' @param dnastringset DNAStringSet. Input set of DNA sequences (e.g., from Biostrings).
#' @param bin_width Integer (default 10). Width of the bins (in bases) to split each sequence for frequency calculation.
#'
#' @return Matrix of A/T frequencies. Rows correspond to individual sequences; columns correspond to bins along the sequence.
#'
#' @examples
#' # Example usage:
#' # calculate_at_frequencies(dnastringset, bin_width = 50)
#'
#' @import Biostrings IRanges stringr
#' @export
calculate_at_frequencies <- function(dnastringset, bin_width = 10) {
  max_seq_length <- max(width(dnastringset))
  bins <- IRanges(
    start = seq(1, max_seq_length - bin_width + 1, by = bin_width),
    width = bin_width
  )
  at_frequencies_matrix <- matrix(
    0,
    nrow = length(dnastringset),
    ncol = length(bins)
  )

  for (i in 1:length(dnastringset)) {
    cat(
      "\r",
      paste0("Calculating: ", round((i / length(dnastringset)) * 100, 2), "%"),
      "\b\b",
      file = stdout()
    )

    for (j in 1:length(bins)) {
      bin_seq <- subseq(dnastringset[i], bins[j])
      at_freq <- sum(str_count(as.character(bin_seq), "A|T")) / bin_width
      at_frequencies_matrix[i, j] <- at_freq
    }
  }
  return(at_frequencies_matrix)
}


#' Generate Random DNA Sequences of Fixed Window Size
#'
#' Randomly samples a defined number of windows of fixed size from a set of DNA sequences.
#'
#' All input sequences are concatenated into one super-sequence before sampling.
#'
#' @param dna_string_set DNAStringSet. Input set of DNA sequences to sample from.
#' @param window_size Integer. Length (in bases) of each random window to sample.
#' @param num_sequences Integer. Number of random sequences (windows) to return.
#'
#' @return DNAStringSet with randomly sampled sequences.
#'
#' @examples
#' # generate_random_sequences(dnastringset, 100, 10)
#'
#' @import Biostrings
#' @export
generate_random_sequences <- function(
  dna_string_set,
  window_size,
  num_sequences
) {
  # Collapse all sequences into a single sequence
  collapsed_sequence <- paste(dna_string_set, collapse = "")

  # Initialize a DNAStringSet to store the sampled sequences
  sampled_sequences <- DNAStringSet()

  while (length(sampled_sequences) < num_sequences) {
    max_start_pos <- width(collapsed_sequence) - window_size
    start_pos <- sample(1:max_start_pos, 1)

    temp <- subseq(
      collapsed_sequence,
      start = start_pos,
      end = start_pos + window_size - 1
    )
    sampled_sequences <- c(sampled_sequences, temp)
  }

  return(sampled_sequences)
}


#' Extract Genomic Windows and Perform Motif Enrichment Analysis
#'
#' Extracts DNA sequences defined by genomic intervals in a \code{GRanges} object
#' and performs motif enrichment analysis using AME (from the MEME Suite).
#'
#' @param GRangeObj GRanges object. Genomic intervals to extract.
#' @param SequenceObj Named DNAStringSet object. Genome sequences, indexed by chromosome name.
#'
#' @details
#' Each interval from \code{GRangeObj} is used to extract the corresponding sequence from \code{SequenceObj}.
#' Motif enrichment analysis is conducted on these sequences via the \code{runAme} function.
#'
#' @return A runAme enrichment result object.
#'
#' @examples
#' # enrich <- subtract_and_enrich(granges, dnasset)
#'
#' @import Biostrings GenomicRanges
#' @export
subtract_and_enrich <- function(GRangeObj, SequenceObj) {
  targets <- DNAStringSet()

  for (i in 1:length(start(GRangeObj))) {
    tryCatch(
      {
        start_temp <- start(GRangeObj)[i]
        end_temp <- end(GRangeObj)[i]
        chr_temp <- as.character(GRangeObj@seqnames[i])

        temp_seq <- subseq(SequenceObj[chr_temp], start_temp, end_temp)

        targets <- c(targets, as.character(temp_seq))
      },
      error = function(e) {
        cat("Error occurred in iteration", i, ": ", conditionMessage(e), "\n")
      }
    )
  }

  enrich <- runAme(
    targets,
    meme_path = "~/Documents/Tools/meme/bin/",
    database = "~/Documents/Tools/meme/Databases/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
    silent = TRUE
  )

  return(enrich)
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
  gene_activity = NULL, # <- NEW (tibble/data.frame with gene_id + is_active)
  active_col = "royalblue", # <- NEW
  inactive_col = "skyblue", # <- NEW
  unknown_col = "grey80", # <- NEW (if a gene_id isn't in the table)
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

    # robust lookup: named logical vector
    # (if duplicated gene_ids exist, keep the first; you can change this if needed)
    ga <- gene_activity[, c("gene_id", "is_active")]
    ga <- ga[!is.na(ga$gene_id), ]
    ga <- ga[!duplicated(ga$gene_id), ]
    active_lookup <- setNames(as.logical(ga$is_active), ga$gene_id)

    is_active_vec <- unname(active_lookup[group_vec]) # aligns to gene_of_interest rows
    # genes missing from lookup become NA -> unknown_col
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
    fill = gene_fill, # <- the key line
    col = gene_col, # <- outlines match fill
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
