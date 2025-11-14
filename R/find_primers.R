#' Analyze genome for target amplicons using primer pairs (optimized - FINAL FIXED)
#'
#' Searches a genome file in FASTA format for potential amplicons using sets of
#' forward and reverse primers. Processes each contig, identifies primer matches,
#' and generates amplicon sequences with quality metrics. Optimized for speed.
#' Duplicate amplicons with identical coordinates and sequences are counted only once.
#'
#' @param genome_file Character string specifying path to genome FASTA file.
#' @param forward_primers Character vector of forward primer sequences.
#' @param reverse_primers Character vector of reverse primer sequences.
#' @param max_mismatch Integer specifying maximum mismatches for primer binding.
#'   Defaults to 1.
#' @param max_amplicon_length Integer specifying maximum amplicon length including
#'   primers. Defaults to 5000.
#' @param min_amplicon_length Integer specifying minimum amplicon length including
#'   primers. Defaults to 100.
#' @param output_dir Character string specifying output directory for results.
#'   Defaults to "contigs_results". Directory is created if it does not exist.
#' @param all Logical. If TRUE processes all contigs in genome. If FALSE stops after
#'   first contig with successful amplicon match. Defaults to FALSE.
#' @param min_contig_length Integer specifying minimum contig length to process.
#'   Contigs shorter than this are skipped. Defaults to 1000.
#'
#' @return List of lists containing details for all unique found amplicons. Each amplicon
#'   entry includes contig name, primer sequences, amplicon coordinates, length,
#'   GC content, and sequence data. Identical amplicons at same position are deduplicated.
#'
#' @examples
#' \dontrun{
#' forward_primers <- c("CGGTTACCAAAACCGGCGTAGAAAGGACCGCCGT")
#' reverse_primers <- c("ACACCGTCCCCTTCTATGCCGGTTTTGGCAACCG")
#' results <- analyze_trichoderma_genome(
#'   genome_file = "genome.fasta",
#'   forward_primers = forward_primers,
#'   reverse_primers = reverse_primers,
#'   all = TRUE
#' )
#' }
#'
#' @export
analyze_trichoderma_genome <- function(genome_file,
                                       forward_primers,
                                       reverse_primers,
                                       max_mismatch = 1,
                                       max_amplicon_length = 5000,
                                       min_amplicon_length = 100,
                                       output_dir = "contigs_results",
                                       all = FALSE,
                                       min_contig_length = 1000) {

  cat("=== GENOME ANALYSIS FOR SPECIFIED PRIMERS (OPTIMIZED - FINAL FIXED) ===\n")
  cat("Genome file:", genome_file, "\n")

  cat("\nLoading all contigs...\n")
  genome_seqs <- readDNAStringSet(genome_file)
  
  genome_seqs <- DNAStringSet(lapply(genome_seqs, function(x) DNAString(toupper(as.character(x)))))

  cat("Total number of contigs/sequences:", length(genome_seqs), "\n")

  contig_lengths <- width(genome_seqs)

  cat("Contig lengths:\n")
  cat("  Longest contig:", max(contig_lengths), "bp\n")
  cat("  Shortest contig:", min(contig_lengths), "bp\n")
  cat("  Average length:", round(mean(contig_lengths), 0), "bp\n")
  cat("  Total genome length:", sum(contig_lengths), "bp\n")

  contigs_to_skip <- sum(contig_lengths < min_contig_length)
  if (contigs_to_skip > 0) {
    cat("  Contigs to skip (<", min_contig_length, "bp):", contigs_to_skip, "\n")
  }

  cat("\nFirst 10 contig names:\n")
  first_names <- names(genome_seqs)[1:min(10, length(genome_seqs))]
  for (i in seq_along(first_names)) {
    cat("  ", i, ":", first_names[i], "(", contig_lengths[i], "bp)\n")
  }
  if (length(genome_seqs) > 10) {
    cat("  ... and", length(genome_seqs) - 10, "more contigs\n")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  all_results <- list()
  total_amplicons <- 0
  processed_contigs <- 0
  contigs_with_hits <- 0
  skipped_contigs <- 0

  cat("\n=== STARTING CONTIG ANALYSIS ===\n")

  start_time <- Sys.time()

  for (i in seq_along(genome_seqs)) {
    contig_name <- names(genome_seqs)[i]
    if (is.null(contig_name) || contig_name == "") {
      contig_name <- paste0("contig_", i)
    }

    contig_length <- contig_lengths[i]

    if (contig_length < min_contig_length) {
      skipped_contigs <- skipped_contigs + 1
      next
    }

    processed_contigs <- processed_contigs + 1

    if (processed_contigs %% 100 == 0 || processed_contigs == sum(contig_lengths >= min_contig_length)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat("Processed", processed_contigs, "contigs (", round(elapsed, 1), "s)...\n")
    }

    contig_results <- analyze_single_contig(
      sequence = genome_seqs[[i]],
      contig_name = contig_name,
      forward_primers = forward_primers,
      reverse_primers = reverse_primers,
      max_mismatch = max_mismatch,
      max_amplicon_length = max_amplicon_length,
      min_amplicon_length = min_amplicon_length,
      verbose = FALSE
    )

    if (length(contig_results) > 0) {
      all_results[[contig_name]] <- contig_results
      total_amplicons <- total_amplicons + length(contig_results)
      contigs_with_hits <- contigs_with_hits + 1

      cat("*** FOUND", length(contig_results), "amplicons in contig:", contig_name, "***\n")

      if (all == FALSE) {
        break
      }
    }
  }

  elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  cat("\n=== GENOME ANALYSIS SUMMARY ===\n")
  cat("Contigs processed:", processed_contigs, "\n")
  cat("Contigs skipped (<", min_contig_length, "bp):", skipped_contigs, "\n")
  cat("Contigs with amplicons:", contigs_with_hits, "\n")
  cat("Total number of amplicons found:", total_amplicons, "\n")
  cat("Total time elapsed:", round(elapsed_total, 2), "seconds\n")

  if (total_amplicons > 0) {
    cat("\nContigs with amplicons:\n")
    for (contig_name in names(all_results)) {
      cat("  ", contig_name, ":", length(all_results[[contig_name]]), "amplicons\n")
    }

    save_contig_results(all_results, output_dir, genome_file, all = all)
  } else {
    cat("\nNo amplicons found in the entire genome.\n")
    cat("Please check:\n")
    cat("  - If the primers are correct\n")
    cat("  - If the max_mismatch parameter is not too restrictive\n")
    cat("  - If the gene is actually present in this genome\n")
  }

  return(all_results)
}




#' Analyze single contig for amplicons (optimized - FINAL FIXED)
#'
#' Helper function that searches a single DNA sequence for all possible amplicons
#' formed by provided primer pairs. Handles all primer orientations on both strands.
#' Duplicate amplicons with identical coordinates and sequences are removed.
#'
#' @param sequence DNAString object representing contig to analyze.
#' @param contig_name Character string with contig identifier.
#' @param forward_primers Character vector of forward primer sequences.
#' @param reverse_primers Character vector of reverse primer sequences.
#' @param max_mismatch Integer for maximum allowed mismatches. Defaults to 1.
#' @param max_amplicon_length Integer for maximum amplicon length. Defaults to 5000.
#' @param min_amplicon_length Integer for minimum amplicon length. Defaults to 100.
#' @param verbose Logical. If TRUE prints detailed analysis information. Defaults to TRUE.
#'
#' @return List of lists where each entry contains detailed amplicon information.
#'   Identical amplicons are deduplicated.
#'
#' @keywords internal
analyze_single_contig <- function(sequence, contig_name,
                                  forward_primers, reverse_primers,
                                  max_mismatch = 1,
                                  max_amplicon_length = 5000,
                                  min_amplicon_length = 100,
                                  verbose = TRUE) {

  sequence <- DNAString(toupper(as.character(sequence)))
  forward_primers <- toupper(forward_primers)
  reverse_primers <- toupper(reverse_primers)

  if (verbose) {
    cat("\n--- Analyzing contig:", contig_name, "---\n")
    cat("Length:", length(sequence), "bp\n")
  }

  amplicons <- list()
  amplicon_signatures <- character()
  amplicon_count <- 0

  for (i in seq_along(forward_primers)) {
    fwd_primer <- forward_primers[i]

    fwd_matches <- matchPattern(fwd_primer, sequence, max.mismatch = max_mismatch)
    fwd_primer_rc <- as.character(reverseComplement(DNAString(fwd_primer)))
    fwd_matches_rc <- matchPattern(fwd_primer_rc, sequence, max.mismatch = max_mismatch)

    if (length(fwd_matches) == 0 && length(fwd_matches_rc) == 0) next

    for (j in seq_along(reverse_primers)) {
      rev_primer <- reverse_primers[j]

      rev_primer_rc <- as.character(reverseComplement(DNAString(rev_primer)))
      rev_matches_rc <- matchPattern(rev_primer_rc, sequence, max.mismatch = max_mismatch)
      
      rev_matches_direct <- matchPattern(rev_primer, sequence, max.mismatch = max_mismatch)

      if (length(rev_matches_rc) == 0 && length(rev_matches_direct) == 0) next

      if (verbose && (length(fwd_matches) > 0 || length(fwd_matches_rc) > 0) && 
          (length(rev_matches_rc) > 0 || length(rev_matches_direct) > 0)) {
        cat("Primer pair", i, "x", j, "- FWD:", length(fwd_matches) + length(fwd_matches_rc),
            "| REV:", length(rev_matches_rc) + length(rev_matches_direct), "\n")
      }

      if (length(fwd_matches) > 0 && length(rev_matches_rc) > 0) {
        for (f in seq_along(fwd_matches)) {
          fwd_start <- start(fwd_matches[f])
          fwd_end <- end(fwd_matches[f])

          for (r in seq_along(rev_matches_rc)) {
            rev_start <- start(rev_matches_rc[r])
            rev_end <- end(rev_matches_rc[r])

            if (rev_start <= fwd_end) next

            amplicon_length <- rev_end - fwd_start + 1
            if (amplicon_length > max_amplicon_length) next
            if (amplicon_length < min_amplicon_length) next

            gene_length <- rev_start - fwd_end - 1
            if (gene_length < 50) next

            amplicon_with_primers <- subseq(sequence, fwd_start, rev_end)
            amplicon_signature <- paste(fwd_start, rev_end, as.character(amplicon_with_primers), sep = "___")

            if (amplicon_signature %in% amplicon_signatures) next

            amplicon_count <- amplicon_count + 1
            amplicon_without_primers <- subseq(sequence, fwd_end + 1, rev_start - 1)

            gene_seq_str <- as.character(amplicon_without_primers)
            gc_count <- sum(strsplit(gene_seq_str, "")[[1]] %in% c("G", "C"))
            gc_content <- (gc_count / nchar(gene_seq_str)) * 100

            amplicon_info <- list(
              contig_name = contig_name,
              amplicon_id = paste0(gsub("[^A-Za-z0-9]", "_", contig_name), "_amp", amplicon_count),
              forward_primer = fwd_primer,
              reverse_primer = rev_primer,
              fwd_start = fwd_start,
              fwd_end = fwd_end,
              rev_start = rev_start,
              rev_end = rev_end,
              amplicon_length = amplicon_length,
              gene_length = gene_length,
              gc_content = gc_content,
              amplicon_with_primers = amplicon_with_primers,
              amplicon_without_primers = amplicon_without_primers
            )

            amplicons[[amplicon_count]] <- amplicon_info
            amplicon_signatures[amplicon_count] <- amplicon_signature

            if (verbose) {
              cat("  *** AMPLICON FOUND (fwd+rev_rc) ***\n")
              cat("  ID:", amplicon_info$amplicon_id, "\n")
              cat("  Amplicon:", amplicon_length, "bp | Gene:", gene_length, "bp\n")
              cat("  GC:", round(gc_content, 1), "% | Position:", fwd_start, "-", rev_end, "\n")
            }
          }
        }
      }

      if (length(fwd_matches) > 0 && length(rev_matches_direct) > 0) {
        for (f in seq_along(fwd_matches)) {
          fwd_start <- start(fwd_matches[f])
          fwd_end <- end(fwd_matches[f])

          for (r in seq_along(rev_matches_direct)) {
            rev_start <- start(rev_matches_direct[r])
            rev_end <- end(rev_matches_direct[r])

            if (rev_start <= fwd_end) next

            amplicon_length <- rev_end - fwd_start + 1
            if (amplicon_length > max_amplicon_length) next
            if (amplicon_length < min_amplicon_length) next

            gene_length <- rev_start - fwd_end - 1
            if (gene_length < 50) next

            amplicon_with_primers <- subseq(sequence, fwd_start, rev_end)
            amplicon_signature <- paste(fwd_start, rev_end, as.character(amplicon_with_primers), sep = "___")

            if (amplicon_signature %in% amplicon_signatures) next

            amplicon_count <- amplicon_count + 1
            amplicon_without_primers <- subseq(sequence, fwd_end + 1, rev_start - 1)

            gene_seq_str <- as.character(amplicon_without_primers)
            gc_count <- sum(strsplit(gene_seq_str, "")[[1]] %in% c("G", "C"))
            gc_content <- (gc_count / nchar(gene_seq_str)) * 100

            amplicon_info <- list(
              contig_name = contig_name,
              amplicon_id = paste0(gsub("[^A-Za-z0-9]", "_", contig_name), "_amp", amplicon_count),
              forward_primer = fwd_primer,
              reverse_primer = rev_primer,
              fwd_start = fwd_start,
              fwd_end = fwd_end,
              rev_start = rev_start,
              rev_end = rev_end,
              amplicon_length = amplicon_length,
              gene_length = gene_length,
              gc_content = gc_content,
              amplicon_with_primers = amplicon_with_primers,
              amplicon_without_primers = amplicon_without_primers
            )

            amplicons[[amplicon_count]] <- amplicon_info
            amplicon_signatures[amplicon_count] <- amplicon_signature

            if (verbose) {
              cat("  *** AMPLICON FOUND (fwd+rev_direct) ***\n")
              cat("  ID:", amplicon_info$amplicon_id, "\n")
              cat("  Amplicon:", amplicon_length, "bp | Gene:", gene_length, "bp\n")
              cat("  GC:", round(gc_content, 1), "% | Position:", fwd_start, "-", rev_end, "\n")
            }
          }
        }
      }

      if (length(fwd_matches_rc) > 0 && length(rev_matches_rc) > 0) {
        for (f in seq_along(fwd_matches_rc)) {
          fwd_start <- start(fwd_matches_rc[f])
          fwd_end <- end(fwd_matches_rc[f])

          for (r in seq_along(rev_matches_rc)) {
            rev_start <- start(rev_matches_rc[r])
            rev_end <- end(rev_matches_rc[r])

            if (rev_start <= fwd_end) next

            amplicon_length <- rev_end - fwd_start + 1
            if (amplicon_length > max_amplicon_length) next
            if (amplicon_length < min_amplicon_length) next

            gene_length <- rev_start - fwd_end - 1
            if (gene_length < 50) next

            amplicon_with_primers <- subseq(sequence, fwd_start, rev_end)
            amplicon_signature <- paste(fwd_start, rev_end, as.character(amplicon_with_primers), sep = "___")

            if (amplicon_signature %in% amplicon_signatures) next

            amplicon_count <- amplicon_count + 1
            amplicon_without_primers <- subseq(sequence, fwd_end + 1, rev_start - 1)

            gene_seq_str <- as.character(amplicon_without_primers)
            gc_count <- sum(strsplit(gene_seq_str, "")[[1]] %in% c("G", "C"))
            gc_content <- (gc_count / nchar(gene_seq_str)) * 100

            amplicon_info <- list(
              contig_name = contig_name,
              amplicon_id = paste0(gsub("[^A-Za-z0-9]", "_", contig_name), "_amp", amplicon_count),
              forward_primer = fwd_primer_rc,
              reverse_primer = rev_primer,
              fwd_start = fwd_start,
              fwd_end = fwd_end,
              rev_start = rev_start,
              rev_end = rev_end,
              amplicon_length = amplicon_length,
              gene_length = gene_length,
              gc_content = gc_content,
              amplicon_with_primers = amplicon_with_primers,
              amplicon_without_primers = amplicon_without_primers
            )

            amplicons[[amplicon_count]] <- amplicon_info
            amplicon_signatures[amplicon_count] <- amplicon_signature

            if (verbose) {
              cat("  *** AMPLICON FOUND (fwd_rc+rev_rc) ***\n")
              cat("  ID:", amplicon_info$amplicon_id, "\n")
              cat("  Amplicon:", amplicon_length, "bp | Gene:", gene_length, "bp\n")
              cat("  GC:", round(gc_content, 1), "% | Position:", fwd_start, "-", rev_end, "\n")
            }
          }
        }
      }

      if (length(fwd_matches_rc) > 0 && length(rev_matches_direct) > 0) {
        for (f in seq_along(fwd_matches_rc)) {
          fwd_start <- start(fwd_matches_rc[f])
          fwd_end <- end(fwd_matches_rc[f])

          for (r in seq_along(rev_matches_direct)) {
            rev_start <- start(rev_matches_direct[r])
            rev_end <- end(rev_matches_direct[r])

            if (rev_start <= fwd_end) next

            amplicon_length <- rev_end - fwd_start + 1
            if (amplicon_length > max_amplicon_length) next
            if (amplicon_length < min_amplicon_length) next

            gene_length <- rev_start - fwd_end - 1
            if (gene_length < 50) next

            amplicon_with_primers <- subseq(sequence, fwd_start, rev_end)
            amplicon_signature <- paste(fwd_start, rev_end, as.character(amplicon_with_primers), sep = "___")

            if (amplicon_signature %in% amplicon_signatures) next

            amplicon_count <- amplicon_count + 1
            amplicon_without_primers <- subseq(sequence, fwd_end + 1, rev_start - 1)

            gene_seq_str <- as.character(amplicon_without_primers)
            gc_count <- sum(strsplit(gene_seq_str, "")[[1]] %in% c("G", "C"))
            gc_content <- (gc_count / nchar(gene_seq_str)) * 100

            amplicon_info <- list(
              contig_name = contig_name,
              amplicon_id = paste0(gsub("[^A-Za-z0-9]", "_", contig_name), "_amp", amplicon_count),
              forward_primer = fwd_primer_rc,
              reverse_primer = rev_primer,
              fwd_start = fwd_start,
              fwd_end = fwd_end,
              rev_start = rev_start,
              rev_end = rev_end,
              amplicon_length = amplicon_length,
              gene_length = gene_length,
              gc_content = gc_content,
              amplicon_with_primers = amplicon_with_primers,
              amplicon_without_primers = amplicon_without_primers
            )

            amplicons[[amplicon_count]] <- amplicon_info
            amplicon_signatures[amplicon_count] <- amplicon_signature

            if (verbose) {
              cat("  *** AMPLICON FOUND (fwd_rc+rev_direct) ***\n")
              cat("  ID:", amplicon_info$amplicon_id, "\n")
              cat("  Amplicon:", amplicon_length, "bp | Gene:", gene_length, "bp\n")
              cat("  GC:", round(gc_content, 1), "% | Position:", fwd_start, "-", rev_end, "\n")
            }
          }
        }
      }
    }
  }

  return(amplicons)
}




#' Save amplicon analysis results to files
#'
#' Processes amplicon analysis results and writes them to FASTA format files.
#' Generates two output files: one containing complete amplicons with primers,
#' and another containing gene sequences without primers.
#'
#' @param all_results List of lists containing amplicons from genome analysis.
#' @param output_dir Character string specifying directory for output files.
#' @param genome_file Character string with path to original genome file.
#' @param all Logical. If TRUE retains all amplicons. If FALSE removes duplicates
#'   with identical genomic coordinates. Defaults to FALSE.
#'
#' @return Invisible NULL.
#'
#' @keywords internal
save_contig_results <- function(all_results, output_dir, genome_file, all = FALSE) {

  all_amplicons <- unlist(all_results, recursive = FALSE)

  if (all == FALSE) {
    coords <- vapply(all_amplicons,
                     function(x) paste(x$contig_name, x$fwd_start, x$rev_end, sep = "_"),
                     character(1))
    unique_indices <- !duplicated(coords)
    all_amplicons <- all_amplicons[unique_indices]
    cat("Removing duplicates by coordinates (all=FALSE)...\n")
  } else {
    cat("Keeping all amplicons (all=TRUE)...\n")
  }

  cat("\nSaving results...\n")

  seqs_with_primers <- list()
  headers_with_primers <- character()

  seqs_without_primers <- list()
  headers_without_primers <- character()

  for (i in seq_along(all_amplicons)) {
    amp <- all_amplicons[[i]]

    header_with <- paste0(amp$amplicon_id,
                          " | Contig:", amp$contig_name,
                          " | Total:", amp$amplicon_length, "bp",
                          " | Gene:", amp$gene_length, "bp",
                          " | GC:", round(amp$gc_content, 1), "%",
                          " | Pos:", amp$fwd_start, "-", amp$rev_end,
                          " | WITH_PRIMERS")

    header_without <- paste0(amp$amplicon_id,
                             " | Contig:", amp$contig_name,
                             " | Gene:", amp$gene_length, "bp",
                             " | GC:", round(amp$gc_content, 1), "%",
                             " | GENE_ONLY")

    seqs_with_primers[[i]] <- amp$amplicon_with_primers
    headers_with_primers[i] <- header_with

    seqs_without_primers[[i]] <- amp$amplicon_without_primers
    headers_without_primers[i] <- header_without
  }

  sequences_with_primers <- DNAStringSet(seqs_with_primers)
  names(sequences_with_primers) <- headers_with_primers

  sequences_without_primers <- DNAStringSet(seqs_without_primers)
  names(sequences_without_primers) <- headers_without_primers

  fasta_with_file <- file.path(output_dir, "amplicons_with_primers.fasta")
  fasta_without_file <- file.path(output_dir, "amplicons_without_primers.fasta")

  writeXStringSet(sequences_with_primers, fasta_with_file)
  writeXStringSet(sequences_without_primers, fasta_without_file)

  cat("=== RESULTS SAVED ===\n")
  cat("Output folder:", output_dir, "\n")
  cat("Amplicons with primers:", fasta_with_file, "\n")
  cat("Amplicons without primers:", fasta_without_file, "\n")
  cat("Total amplicons saved:", length(all_amplicons), "\n\n")

  if (length(all_amplicons) > 0) {
    cat("=== DETAILS OF FOUND AMPLICONS ===\n")
    for (i in seq_along(all_amplicons)) {
      amp <- all_amplicons[[i]]
      cat("Amplicon", i, ":\n")
      cat("  ID:", amp$amplicon_id, "\n")
      cat("  Contig:", amp$contig_name, "\n")
      cat("  Position:", amp$fwd_start, "-", amp$rev_end, "\n")
      cat("  Total length:", amp$amplicon_length, "bp\n")
      cat("  Gene length:", amp$gene_length, "bp\n")
      cat("  GC content:", round(amp$gc_content, 1), "%\n")
      cat("  Gene sequence (first 60 bp):\n")
      gene_preview <- as.character(subseq(amp$amplicon_without_primers, 1, min(60, length(amp$amplicon_without_primers))))
      cat("   ", gene_preview, "\n\n")
    }
  }

  invisible(NULL)
}

#' Analyze Multiple Trichoderma Loci from Genome
#'
#' @description
#' Performs primer-based amplicon identification for multiple Trichoderma loci
#' (TEF1, RPB2, TEF3, TUB2, LNS2, ACT, PGK, ITS) from a genome sequence file.
#' Results for each locus are saved to separate output directories.
#'
#' @param genome_file Character string specifying the path to the genome FASTA file.
#'
#'
#' Analysis parameters are fixed across all loci:
#' \itemize{
#'   \item Maximum primer mismatches: 1
#'   \item Maximum amplicon length: 5000 bp
#'   \item Minimum amplicon length: 100 bp
#' }
#'
#' @return
#' A named list containing analysis results for each locus.
#' Names correspond to locus identifiers (TEF1, RPB2, etc.).
#' Each element contains the output from \code{\link{analyze_trichoderma_genome}}.
#'
#' @examples
#' \dontrun{
#' # Run with default genome file
#' results <- all_amplicon_identification()
#'
#' # Run with custom genome file
#' results <- all_amplicon_identification(
#'   genome_file = "/path/to/custom/genome.fasta"
#' )
#'
#' # Access results for specific locus
#' tef1_results <- results$TEF1
#' its_results <- results$ITS
#' }
#'
#' @seealso
#' \code{\link{analyze_trichoderma_genome}} for the underlying analysis function
#'
#' @keywords internal
#'
#' @export
#'
all_amplicon_identification <- function(genome_file = "file.fasta") {
  
  loci <- list(
    TEF1 = TEF1,
    RPB2 = RPB2,
    TEF3 = TEF3,
    TUB2 = TUB2,
    LNS2 = LNS2,
    ACT = ACT,
    PGK = PGK,
    ITS = ITS
  )
  
  common_params <- list(
    genome_file = genome_file,
    max_mismatch = 1,
    max_amplicon_length = 5000,
    min_amplicon_length = 100,
    all = TRUE
  )
  
  all_results <- list()
  
  for (locus_name in names(loci)) {
    message("Analyzing: ", locus_name)
    
    results <- analyze_trichoderma_genome(
      forward_primers = loci[[locus_name]],
      reverse_primers = loci[[locus_name]],
      output_dir = paste0(locus_name, "_results"),
      max_mismatch = common_params$max_mismatch,
      max_amplicon_length = common_params$max_amplicon_length,
      min_amplicon_length = common_params$min_amplicon_length,
      genome_file = common_params$genome_file,
      all = common_params$all
    )
    
    all_results[[locus_name]] <- results
  }
  
  message("Analysis of all loci completed.")
  return(all_results)
}

#' Combine Trichoderma Amplicons into Single FASTA
#'
#' @description
#' Combines amplicons from TEF1, RPB2, TEF3, LNS2, ACT, PGK, and ITS loci
#' into a single FASTA file for TrichofindR database comparison.
#' Excludes TUB2 due to duplication in Trichoderma genomes.
#' Creates IMLDS_identification folder containing the combined FASTA.
#'
#' @param genome_file Character string specifying the path to the genome FASTA file.
#'
#' @return
#' Invisibly returns path to the combined FASTA file (ultra.fasta).
#' Creates IMLDS_identification folder with ultra.fasta inside.
#'
#' @examples
#' \dontrun{
#' # Run with default genome file
#' create_imlds_identification()
#'
#' # Run with custom genome file
#' create_imlds_identification(
#'   genome_file = "/path/to/custom/genome.fasta"
#' )
#' }
#'
#' @seealso
#' \code{\link{analyze_trichoderma_genome}} for the underlying analysis function
#'
#' @keywords internal
#'
#' @export
#'

IMLDS_identification <- function(genome_file = "/path/to/your/genome.fasta") {
  
  loci <- list(
    TEF1 = TEF1,
    RPB2 = RPB2,
    TEF3 = TEF3,
    LNS2 = LNS2,
    ACT = ACT,
    PGK = PGK,
    ITS = ITS
  )
  
  common_params <- list(
    genome_file = genome_file,
    max_mismatch = 1,
    max_amplicon_length = 5000,
    min_amplicon_length = 100,
    all = TRUE
  )
  
  all_results <- list()
  locus_folders <- c()
  
  message("Starting amplicon identification for IMLDS loci...")
  
  for (locus_name in names(loci)) {
    message("Analyzing: ", locus_name)
    
    results <- analyze_trichoderma_genome(
      forward_primers = loci[[locus_name]],
      reverse_primers = loci[[locus_name]],
      output_dir = paste0(locus_name, "_results"),
      max_mismatch = common_params$max_mismatch,
      max_amplicon_length = common_params$max_amplicon_length,
      min_amplicon_length = common_params$min_amplicon_length,
      genome_file = common_params$genome_file,
      all = common_params$all
    )
    
    all_results[[locus_name]] <- results
    locus_folders <- c(locus_folders, paste0(locus_name, "_results"))
  }
  
  message("All loci analyzed. Combining sequences...")
  
  output_dir <- "IMLDS_identification"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  combined_sequences <- NULL
  
  for (locus_folder in locus_folders) {
    fasta_files <- list.files(
      locus_folder,
      pattern = "\\.fasta$|\\.fa$",
      full.names = TRUE
    )
    
    if (length(fasta_files) > 0) {
      tryCatch({
        sequences <- Biostrings::readDNAStringSet(fasta_files[1])
        
        names(sequences) <- paste0(
          names(sequences),
          "_",
          sub("_results", "", locus_folder)
        )
        
        combined_sequences <- c(combined_sequences, sequences)
        message("Added ", length(sequences), " sequences from ", locus_folder)
      }, error = function(e) {
        warning("Could not read FASTA from ", locus_folder, ": ", e$message)
      })
    } else {
      warning("No FASTA file found in ", locus_folder)
    }
  }
  
  if (!is.null(combined_sequences)) {
    output_fasta <- file.path(output_dir, "ultra.fasta")
    Biostrings::writeXStringSet(combined_sequences, output_fasta)
    message("Combined FASTA written to: ", output_fasta)
    message("Total sequences: ", length(combined_sequences))
  } else {
    warning("No sequences were combined!")
  }
  
  message("IMLDS identification complete!")
  message("Individual locus folders preserved: ", paste(locus_folders, collapse = ", "))
  invisible(file.path(output_dir, "ultra.fasta"))
}
