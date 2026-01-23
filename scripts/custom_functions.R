

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
  bins <- IRanges(start = seq(1, max_seq_length - bin_width + 1, by = bin_width), width = bin_width)
  at_frequencies_matrix <- matrix(0, nrow = length(dnastringset), ncol = length(bins))
  
  for (i in 1:length(dnastringset)) {
    cat("\r", paste0("Calculating: ", round((i / length(dnastringset)) * 100, 2), "%"), "\b\b", file = stdout())
    
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
generate_random_sequences <- function(dna_string_set, window_size, num_sequences) {
  # Collapse all sequences into a single sequence
  collapsed_sequence <- paste(dna_string_set, collapse = "")
  
  # Initialize a DNAStringSet to store the sampled sequences
  sampled_sequences <- DNAStringSet()
  
  while (length(sampled_sequences) < num_sequences) {
    max_start_pos <- width(collapsed_sequence) - window_size
    start_pos <- sample(1:max_start_pos, 1)
    
    temp <- subseq(collapsed_sequence, start = start_pos, end = start_pos + window_size - 1)
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
    tryCatch({
      start_temp <- start(GRangeObj)[i]
      end_temp <- end(GRangeObj)[i]
      chr_temp <- as.character(GRangeObj@seqnames[i])
      
      temp_seq <- subseq(SequenceObj[chr_temp], start_temp, end_temp)
      
      targets <- c(targets, as.character(temp_seq))
    }, error = function(e) {
      cat("Error occurred in iteration", i, ": ", conditionMessage(e), "\n")
    })
  }
  
  enrich <- runAme(targets,
                   meme_path = "~/Documents/Tools/meme/bin/",
                   database = "~/Documents/Tools/meme/Databases/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
                   silent = TRUE)
  
  return(enrich)
}
