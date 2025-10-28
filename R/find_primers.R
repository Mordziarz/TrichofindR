#' Analyzes a genome for target amplicons using primer pairs.
#'
#' This function takes a genome file (e.g., in FASTA format) and a set of
#' forward and reverse primers to search for potential amplicons. It processes
#' each contig, calls `analyze_single_contig` to find matches, and then
#' summarizes and saves the results.
#'
#' @param genome_file A character string specifying the path to the genome file (e.g., a FASTA file).
#' @param forward_primers A character vector of forward primer sequences.
#' @param reverse_primers A character vector of reverse primer sequences.
#' @param max_mismatch An integer specifying the maximum number of mismatches allowed for primer binding. Defaults to 1.
#' @param max_amplicon_length An integer specifying the maximum allowed length of the amplicon (including primers). Defaults to 5000.
#' @param min_amplicon_length An integer specifying the minimum allowed length of the amplicon (including primers). Defaults to 100.
#' @param output_dir A character string specifying the directory to save the output files. Defaults to "contigs_results".
#' @param all A logical value. If `TRUE`, the function processes all contigs. If `FALSE`, it stops after the first contig with a successful amplicon match. Defaults to `FALSE`.
#'
#' @return A list of lists, where each sublist contains details about the found amplicons.
#' @examples

analyze_trichoderma_genome <- function(genome_file,
                                         forward_primers = c(TEF1),
                                         reverse_primers = c(TEF1),
                                         max_mismatch = 1,
                                         max_amplicon_length = 5000,
                                         min_amplicon_length = 100,
                                         output_dir = "contigs_results",
                                         all=FALSE) {

  cat("=== GENOME ANALYSIS FOR SPECIFIED PRIMERS ===\n")
  cat("Genome file:", genome_file, "\n")

  cat("\nLoading all contigs...\n")
  genome_seqs <- readDNAStringSet(genome_file)

  cat("Total number of contigs/sequences:", length(genome_seqs), "\n")

  contig_lengths <- width(genome_seqs)
  cat("Contig lengths:\n")
  cat("  Longest contig:", max(contig_lengths), "bp\n")
  cat("  Shortest contig:", min(contig_lengths), "bp\n")
  cat("  Average length:", round(mean(contig_lengths), 0), "bp\n")
  cat("  Total genome length:", sum(contig_lengths), "bp\n")

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

  cat("\n=== STARTING CONTIG ANALYSIS ===\n")

  for (i in seq_along(genome_seqs)) {
    contig_name <- names(genome_seqs)[i]
    if (is.null(contig_name) || contig_name == "") {
      contig_name <- paste0("contig_", i)
    }

    processed_contigs <- processed_contigs + 1

    if (processed_contigs %% 100 == 0 || processed_contigs == length(genome_seqs)) {
      cat("Processed", processed_contigs, "of", length(genome_seqs), "contigs...\n")
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

      if(all==FALSE){
      break
      }
    }
  }

  cat("\n=== GENOME ANALYSIS SUMMARY ===\n")
  cat("Contigs analyzed:", processed_contigs, "\n")
  cat("Contigs with amplicons:", contigs_with_hits, "\n")
  cat("Total number of amplicons found:", total_amplicons, "\n")

  if (total_amplicons > 0) {
    cat("\nContigs with amplicons:\n")
    for (contig_name in names(all_results)) {
      cat("  ", contig_name, ":", length(all_results[[contig_name]]), "amplicons\n")
    }

    save_contig_results(all_results, output_dir, genome_file)
  } else {
    cat("\nNo amplicons found in the entire genome.\n")
    cat("Please check:\n")
    cat("  - If the primers are correct\n")
    cat("  - If the max_mismatch parameter is not too restrictive\n")
    cat("  - If the gene is actually present in this genome\n")
  }

  return(all_results)
}

#' Analyzes a single contig for amplicons.
#'
#' This is a helper function that takes a single DNA sequence (a contig) and searches
#' for all possible amplicons formed by the provided primer pairs. It considers
#' mismatches and filters amplicons based on length criteria.
#'
#' @param sequence A `DNAString` object representing the contig to be analyzed.
#' @param contig_name A character string with the name of the contig.
#' @param forward_primers A character vector of forward primer sequences.
#' @param reverse_primers A character vector of reverse primer sequences.
#' @param max_mismatch An integer for the maximum number of mismatches. Defaults to 1.
#' @param max_amplicon_length An integer for the maximum amplicon length. Defaults to 5000.
#' @param min_amplicon_length An integer for the minimum amplicon length. Defaults to 100.
#' @param verbose A logical value. If `TRUE`, prints detailed information about the analysis. Defaults to `TRUE`.
#'
#' @return A list of lists, where each sublist contains detailed information for a found amplicon.
analyze_single_contig <- function(sequence, contig_name,
                                  forward_primers, reverse_primers,
                                  max_mismatch = 1,
                                  max_amplicon_length = 5000,
                                  min_amplicon_length = 100,
                                  verbose = TRUE) {

  if (verbose) {
    cat("\n--- Analyzing contig:", contig_name, "---\n")
    cat("Length:", length(sequence), "bp\n")
  }

  amplicons <- list()
  amplicon_count <- 0

  for (i in seq_along(forward_primers)) {
    for (j in seq_along(reverse_primers)) {

      fwd_primer <- forward_primers[i]
      rev_primer <- reverse_primers[j]

      fwd_matches <- matchPattern(fwd_primer, sequence, max.mismatch = max_mismatch)

      rev_primer_rc <- reverseComplement(DNAString(rev_primer))
      rev_matches <- matchPattern(rev_primer_rc, sequence, max.mismatch = max_mismatch)

      if (verbose && (length(fwd_matches) > 0 || length(rev_matches) > 0)) {
        cat("Primer pair", i, "x", j, "- Forward:", length(fwd_matches),
            "| Reverse:", length(rev_matches), "\n")
      }

      if (length(fwd_matches) > 0 && length(rev_matches) > 0) {

        for (f in seq_along(fwd_matches)) {
          for (r in seq_along(rev_matches)) {

            fwd_start <- start(fwd_matches[f])
            fwd_end <- end(fwd_matches[f])
            rev_start <- start(rev_matches[r])
            rev_end <- end(rev_matches[r])

            if (rev_start > fwd_end) {
              amplicon_length <- rev_end - fwd_start + 1
              gene_length <- rev_start - fwd_end - 1

              if (amplicon_length >= min_amplicon_length &&
                  amplicon_length <= max_amplicon_length &&
                  gene_length >= 50) {

                amplicon_count <- amplicon_count + 1

                amplicon_with_primers <- subseq(sequence, fwd_start, rev_end)
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

                if (verbose) {
                  cat("  *** AMPLICON FOUND ***\n")
                  cat("  ID:", amplicon_info$amplicon_id, "\n")
                  cat("  Amplicon:", amplicon_length, "bp | Gene:", gene_length, "bp\n")
                  cat("  GC:", round(gc_content, 1), "% | Position:", fwd_start, "-", rev_end, "\n")
                }
              }
            }
          }
        }
      }
    }
  }

  return(amplicons)
}

#' Saves amplicon analysis results to files.
#'
#' This helper function takes the results from `analyze_trichoderma_genome` and
#' writes them to various output files, including FASTA files for amplicons,
#' a detailed CSV summary, and a text report.
#'
#' @param all_results A list of lists containing all found amplicons from the genome analysis.
#' @param output_dir A character string specifying the directory to save the output files.
#' @param genome_file A character string with the path to the original genome file, used for reporting.
#'
#' @return NULL. The function saves files to the specified directory and prints a summary.
save_contig_results <- function(all_results, output_dir, genome_file) {

  all_amplicons <- unlist(all_results, recursive = FALSE)

  coords <- vapply(all_amplicons,
                   function(x) paste(x$contig_name, x$fwd_start, x$rev_end, sep = "_"),
                   character(1))
  unique_indices <- !duplicated(coords)
  all_amplicons <- all_amplicons[unique_indices]

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
  fasta_without_file <- file.path(output_dir, "genes_without_primers.fasta")

  writeXStringSet(sequences_with_primers, fasta_with_file)
  writeXStringSet(sequences_without_primers, fasta_without_file)

  # CSV Summary
  summary_data <- data.frame(
    Amplicon_ID = sapply(all_amplicons, function(x) x$amplicon_id),
    Contig_Name = sapply(all_amplicons, function(x) x$contig_name),
    Forward_Primer = sapply(all_amplicons, function(x) x$forward_primer),
    Reverse_Primer = sapply(all_amplicons, function(x) x$reverse_primer),
    Forward_Start = sapply(all_amplicons, function(x) x$fwd_start),
    Forward_End = sapply(all_amplicons, function(x) x$fwd_end),
    Reverse_Start = sapply(all_amplicons, function(x) x$rev_start),
    Reverse_End = sapply(all_amplicons, function(x) x$rev_end),
    Amplicon_Length = sapply(all_amplicons, function(x) x$amplicon_length),
    Gene_Length = sapply(all_amplicons, function(x) x$gene_length),
    GC_Content = sapply(all_amplicons, function(x) round(x$gc_content, 2)),
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(output_dir, "all_amplicons_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)

  # Contig statistics
  contig_stats <- data.frame(
    Contig_Name = names(all_results),
    Number_of_Amplicons = sapply(all_results, length),
    stringsAsFactors = FALSE
  )

  contig_stats_file <- file.path(output_dir, "contig_statistics.csv")
  write.csv(contig_stats, contig_stats_file, row.names = FALSE)

  # Final report
  report_lines <- c(
    "=== GENOME ANALYSIS REPORT ===",
    paste("Analysis date:", Sys.time()),
    paste("Genome file:", basename(genome_file)),
    "",
    "FOUND AMPLİCONS:",
    paste("Total number of amplicons:", length(all_amplicons)),
    paste("Number of contigs with amplicons:", length(all_results)),
    "",
    "LENGTH STATISTICS:",
    paste("Average amplicon length:", round(mean(summary_data$Amplicon_Length), 1), "bp"),
    paste("Average gene length:", round(mean(summary_data$Gene_Length), 1), "bp"),
    paste("Average GC content:", round(mean(summary_data$GC_Content), 1), "%"),
    "",
    "CONTIGS WITH AMPLİCONS:"
  )

  for (contig_name in names(all_results)) {
    num_amps <- length(all_results[[contig_name]])
    report_lines <- append(report_lines, paste("  ", contig_name, ":", num_amps, "amplicons"))
  }

  report_lines <- append(report_lines, c(
    "",
    "OUTPUT FILES:",
    paste("  - Amplicons with primers:", basename(fasta_with_file)),
    paste("  - Genes without primers:", basename(fasta_without_file)),
    paste("  - Detailed summary:", basename(summary_file)),
    paste("  - Contig statistics:", basename(contig_stats_file))
  ))

  report_file <- file.path(output_dir, "ANALYSIS_REPORT.txt")
  writeLines(report_lines, report_file)

  cat("=== RESULTS SAVED ===\n")
  cat("Output folder:", output_dir, "\n")
  cat("Amplicons with primers:", fasta_with_file, "\n")
  cat("Genes without primers:", fasta_without_file, "\n")
  cat("Summary:", summary_file, "\n")
  cat("Report:", report_file, "\n\n")

  if (length(all_amplicons) > 0) {
    cat("=== DETAILS OF FOUND AMPLİCONS ===\n")
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
}

# TEF1 primers
TEF1 <- c("CATCGAGAAGTTCGAGAAGG","AACTTGCAGGCAATGTGG")

# RPB2 primers
RPB2 <- c(  "GATGATAGAGATCACTTTGG","GATGATAGAGATCACTTCGG","GATGATAGAGATCACTATGG","GATGATAGAGATCACTATCG",
                          "GATGATAGTGATCACTTTGG","GATGATAGTGATCACTTCGG","GATGATAGTGATCACTATGG","GATGATAGTGATCACTATCG",
                          "GATGATCGAGATCACTTTGG","GATGATCGAGATCACTTCGG","GATGATCGAGATCACTATGG","GATGATCGAGATCACTATCG",
                          "GATGATCGTGATCACTTTGG","GATGATCGTGATCACTTCGG","GATGATCGTGATCACTATGG","GATGATCGTGATCACTATCG",
                          "GATGACAGAGATCACTTTGG","GATGACAGAGATCACTTCGG","GATGACAGAGATCACTATGG","GATGACAGAGATCACTATCG",
                          "GATGACAGTGATCACTTTGG","GATGACAGTGATCACTTCGG","GATGACAGTGATCACTATGG","GATGACAGTGATCACTATCG",
                          "GATGACCGAGATCACTTTGG","GATGACCGAGATCACTTCGG","GATGACCGAGATCACTATGG","GATGACCGAGATCACTATCG",
                          "GATGACCGTGATCACTTTGG","GATGACCGTGATCACTTCGG","GATGACCGTGATCACTATGG","GATGACCGTGATCACTATCG",
                          "GACGATAGAGATCACTTTGG","GACGATAGAGATCACTTCGG","GACGATAGAGATCACTATGG","GACGATAGAGATCACTATCG",
                          "GACGATAGTGATCACTTTGG","GACGATAGTGATCACTTCGG","GACGATAGTGATCACTATGG","GACGATAGTGATCACTATCG",
                          "GACGATCGAGATCACTTTGG","GACGATCGAGATCACTTCGG","GACGATCGAGATCACTATGG","GACGATCGAGATCACTATCG",
                          "GACGATCGTGATCACTTTGG","GACGATCGTGATCACTTCGG","GACGATCGTGATCACTATGG","GACGATCGTGATCACTATCG",
                          "GACGACAGAGATCACTTTGG","GACGACAGAGATCACTTCGG","GACGACAGAGATCACTATGG","GACGACAGAGATCACTATCG",
                          "GACGACAGTGATCACTTTGG","GACGACAGTGATCACTTCGG","GACGACAGTGATCACTATGG","GACGACAGTGATCACTATCG",
                          "GACGACCGAGATCACTTTGG","GACGACCGAGATCACTTCGG","GACGACCGAGATCACTATGG","GACGACCGAGATCACTATCG",
                          "GACGACCGTGATCACTTTGG","GACGACCGTGATCACTTCGG","GACGACCGTGATCACTATGG","GACGACCGTGATCACTATCG",
                          "CCCATAGCTTGTTTACCCAT","CCCATAGCTTGTTTGCCCAT","CCCATAGCTTGCTTACCCAT","CCCATAGCTTGCTTGCCCAT",
                          "CCCATGGCTTGTTTACCCAT","CCCATGGCTTGTTTGCCCAT","CCCATGGCTTGCTTACCCAT","CCCATGGCTTGCTTGCCCAT"
)
# ITS primers
ITS <-  c("GGAAGTAAAAGTCGTAACAAGG","TCCTCCGCTTATTGATATGC")

# TOP1
TOP1 <- c("TGTAAAACGACGGCCAGTACGATACTGCCAAGGTTTTCCGTACATACAACGC", "TGTAAAACGACGGCCAGTACGATACTGCCAAGGTTTTCCGTACCTACAACGC", "TGTAAAACGACGGCCAGTACGATACTGCCAAGGTTTTCCGTACTTACAACGC",
           "CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTAATAGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTAATGGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTGATAGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTGATGGCCCA",
            "CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTAATAGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTAATGGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTGATAGCCCA", "CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTGATGGCCCA")

# PGK
PGK <- c(
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGATTCCACGTTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGCTTCCACGTTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACCTGCGTTTCCACGTTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGATTCCACGTTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGCTTCCACGTTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACACCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACACCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACACTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACACTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACATCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACATCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACATTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACATTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGCCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGCCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGCTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGCTGAGGAGGAGGG",
"TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGTCGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGTCGAGGAGGAGGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGTTGAGGAGGAAGG","TGTAAAACGACGGCCAGTACGATGAGAACTTGCGTTTCCACGTTGAGGAGGAGGG",
"CAGGAAACAGCTATGACTCGAAGACACCAGGGGGACCGTTCCA","CAGGAAACAGCTATGACTCGAAGACACCGGGGGGACCGTTCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTAATAGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTAATGGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTGATAGCCCA",
"CAGGAAACAGCTATGACCCAGTCCTCGTCAACAGACTTGATGGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTAATAGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTAATGGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTGATAGCCCA","CAGGAAACAGCTATGACCCAGTCCTCGTCAACTGACTTGATGGCCCA",
  "CAGGAAACAGCTATGACCTTCTTGAAGGTGAAAGCCAT","CAGGAAACAGCTATGACCTTCTTGAAGGTGAAGGCCAT", "GTCGGTGCCCTGCCAACCATCAA","GTCGGTGCCCTGCCCACCATCAA","GTCGGTGCCTTGCCAACCATCAA","GTCGGTGCCTTGCCCACCATCAA",
"GTCGGTGCTCTGCCAACCATCAA","GTCGGTGCTCTGCCCACCATCAA","GTCGGTGCTTTGCCAACCATCAA","GTCGGTGCTTTGCCCACCATCAA",
"GTCGCTGCCCTGCCAACCATCAA","GTCGCTGCCCTGCCCACCATCAA","GTCGCTGCCTTGCCAACCATCAA","GTCGCTGCCTTGCCCACCATCAA",
"GTCGCTGCTCTGCCAACCATCAA","GTCGCTGCTCTGCCCACCATCAA","GTCGCTGCTTTGCCAACCATCAA","GTCGCTGCTTTGCCCACCATCAA",
"GTTGGTGCCCTGCCAACCATCAA","GTTGGTGCCCTGCCCACCATCAA","GTTGGTGCCTTGCCAACCATCAA","GTTGGTGCCTTGCCCACCATCAA",
"GTTGGTGCTCTGCCAACCATCAA","GTTGGTGCTCTGCCCACCATCAA","GTTGGTGCTTTGCCAACCATCAA","GTTGGTGCTTTGCCCACCATCAA",
"GTTGCTGCCCTGCCAACCATCAA","GTTGCTGCCCTGCCCACCATCAA","GTTGCTGCCTTGCCAACCATCAA","GTTGCTGCCTTGCCCACCATCAA",
"GTTGCTGCTCTGCCAACCATCAA","GTTGCTGCTCTGCCCACCATCAA","GTTGCTGCTTTGCCAACCATCAA","GTTGCTGCTTTGCCCACCATCAA","ATCTTGTCAGAAACCTTAGCACC","ATCTTGTCAGAAACCTTGGCACC","ATCTTGTCAGAGACCTTAGCACC","ATCTTGTCAGAGACCTTGGCACC",
"ATCTTGTCAGCAACCTTAGCACC","ATCTTGTCAGCAACCTTGGCACC","ATCTTGTCAGCGACCTTAGCACC","ATCTTGTCAGCGACCTTGGCACC",
"ATCTTGTCGGAAACCTTAGCACC","ATCTTGTCGGAAACCTTGGCACC","ATCTTGTCGGAGACCTTAGCACC","ATCTTGTCGGAGACCTTGGCACC",
"ATCTTGTCGGCAACCTTAGCACC","ATCTTGTCGGCAACCTTGGCACC","ATCTTGTCGGCGACCTTAGCACC","ATCTTGTCGGCGACCTTGGCACC","GTCGACTTCAACGTCCC","GTCGACTTCAACGTTCC","GTCGACTTCAATGTCCC","GTCGACTTCAATGTTCC",
"GTCGATTTCAACGTCCC","GTCGATTTCAACGTTCC","GTCGATTTCAATGTCCC","GTCGATTTCAATGTTCC",
"GTTGACTTCAACGTCCC","GTTGACTTCAACGTTCC","GTTGACTTCAATGTCCC","GTTGACTTCAATGTTCC",
"GTTGATTTCAACGTCCC","GTTGATTTCAACGTTCC","GTTGATTTCAATGTCCC","GTTGATTTCAATGTTCC","ACACCAGGAGGACCGTTCCA","ACACCAGGAGGGCCGTTCCA","ACACCAGGGGGACCGTTCCA","ACACCAGGGGGGCCGTTCCA",
"ACACCAGGTGGACCGTTCCA","ACACCAGGTGGGCCGTTCCA","ACACCGGGAGGACCGTTCCA","ACACCGGGAGGGCCGTTCCA",
"ACACCGGGGGGACCGTTCCA","ACACCGGGGGGGCCGTTCCA","ACACCGGGTGGACCGTTCCA","ACACCGGGTGGGCCGTTCCA",
"ACACCTGGAGGACCGTTCCA","ACACCTGGAGGGCCGTTCCA","ACACCTGGGGGACCGTTCCA","ACACCTGGGGGGCCGTTCCA",
"ACACCTGGTGGACCGTTCCA","ACACCTGGTGGGCCGTTCCA"
)                     

LSU <- c("ACCCGCTGAACTTAAGC","CCGTGTTTCAAGACGGG","TCCTGAGGGAAACTTC","CGCCGTTACTAGGAAAGTTAA","GTACCCGCTGAACTTAAGC")
