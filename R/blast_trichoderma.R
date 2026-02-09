#' BLAST Search Against Trichoderma Reference Sequences
#'
#' @description
#' Performs a BLAST nucleotide search of query DNA sequences against a reference
#' database of Trichoderma sequences. This function is designed for fungal species
#' identification and phylogenetic classification of Trichoderma isolates based on
#' sequence homology across various genomic regions (ITS, ribosomal genes, protein-coding
#' genes, and other molecular markers).
#'
#' @details
#' The function implements a complete BLAST workflow:
#' \enumerate{
#'   \item Prepares reference sequences by removing whitespace characters
#'   \item Creates a DNAStringSet object from the reference sequences
#'   \item Builds a local BLAST nucleotide database
#'   \item Performs blastn search against the reference database
#'   \item Filters and sorts results by percent identity (pident)
#' }
#'
#' The output provides ranked matches ordered by sequence similarity, facilitating
#' the identification of Trichoderma sequences regardless of the genomic region used.
#' Compatible with any DNA sequences including ITS regions, ribosomal genes (18S, 28S),
#' protein-coding genes, and custom molecular markers.
#'
#' @param query_sequence Character string specifying the path to a FASTA file containing
#'   query DNA sequence(s) to be searched against the reference database. If empty string
#'   (default), the function will fail at the readDNAStringSet step unless a file is provided.
#'   Required: FASTA format file with proper formatting. Can contain single or multiple sequences.
#'
#' @param reference_sequences DNAStringSet or character vector containing reference
#'   Trichoderma sequences. Defaults to \code{trichoderma_reference_sequences}, which should
#'   be a pre-loaded DNAStringSet object in the global environment. Can accept sequences
#'   from any genomic region or molecular marker. Can accept either raw DNA strings or
#'   already-formatted DNAStringSet objects.
#'
#' @return
#' A data frame of class \code{data.frame} containing BLAST results sorted by percent
#' identity (pident) in descending order. Expected columns include:
#' \describe{
#'   \item{query}{Query sequence identifier}
#'   \item{subject}{Subject sequence identifier from reference database}
#'   \item{pident}{Percent identity of the alignment (numeric)}
#'   \item{evalue}{E-value of the match}
#'   \item{bitscore}{Bit score of the alignment}
#'   \item{length}{Length of the alignment}
#'   \item{mismatch}{Number of mismatches}
#'   \item{gapopen}{Number of gap openings}
#' }
#'
#' @examples
#' \dontrun{
#' # Load reference sequences
#' data(trichoderma_reference_sequences)
#'
#' # Perform BLAST search with query sequence (any genomic region)
#' results <- trichoderma_blast(
#'   query_sequence = "path/to/query_sequences.fasta",
#'   reference_sequences = trichoderma_reference_sequences
#' )
#'
#' # Display top matches
#' head(results, 10)
#'
#' # Access best match
#' best_match <- results[1, ]
#'
#' # Filter for high-confidence matches (>95% identity)
#' high_confidence <- results[results$pident >= 95, ]
#' }
#'
#' @seealso
#' \code{\link[Biostrings]{DNAStringSet}} for DNA sequence handling
#' \code{\link[Biostrings]{readDNAStringSet}} for reading FASTA files
#' \code{\link[Biostrings]{writeXStringSet}} for writing sequence files
#' \code{\link[rBLAST]{makeblastdb}} for database creation
#' \code{\link[rBLAST]{blast}} for BLAST object creation
#' \code{\link[rBLAST]{predict}} for sequence similarity searching
#'
#' @importFrom Biostrings DNAStringSet readDNAStringSet writeXStringSet
#' @importFrom rBLAST makeblastdb blast predict
#' @export
#'
#' @author Mateusz Mazdziarz
#' @keywords sequence-matching fungal-identification BLAST nucleotide-search

trichoderma_blast <- function(query_sequence = "",
                               reference_sequences = ITS_reference_sequences) {

  if (!is.character(query_sequence) || nchar(query_sequence) == 0) {
    stop("query_sequence must be a non-empty string pointing to a FASTA file")
  }

  if (!file.exists(query_sequence)) {
    stop(paste("Query sequence file not found:", query_sequence))
  }

  if (inherits(reference_sequences, "XStringSet")) {
    sequences_vector <- as.character(reference_sequences)
  } else {
    sequences_vector <- unlist(reference_sequences)
  }

  sequences_vector <- gsub("\\s+", "", sequences_vector)

  if (length(sequences_vector) == 0) {
    stop("Reference sequences vector is empty")
  }

  dna_sequences <- DNAStringSet(sequences_vector)

  subject_db_path <- "subject_sequences.fasta"
  writeXStringSet(dna_sequences, filepath = subject_db_path)

  makeblastdb(subject_db_path, dbtype = "nucl")

  blast_db <- blast(db = subject_db_path, type = "blastn")

  tryCatch(
    {
      query_sequences <- readDNAStringSet(query_sequence)
    },
    error = function(e) {
      stop(paste("Error reading query sequences from", query_sequence, ":", e$message))
    }
  )

  if (length(query_sequences) == 0) {
    stop("No valid sequences found in query file")
  }

  blast_result <- predict(blast_db, query_sequences)

  results_aggregated <- blast_result %>%
  dplyr::mutate(matching_nts = pident * length / 100) %>%
  dplyr::group_by(qseqid, sseqid) %>%
  dplyr::summarize(
    total_length = sum(length),                 
    total_matching_nts = sum(matching_nts),    
    best_evalue = min(evalue),
    max_bitscore = max(bitscore),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    average_pident_weighted = (total_matching_nts / total_length) * 100
  ) %>%
  dplyr::arrange(desc(average_pident_weighted), desc(total_length))

  return(results_aggregated)
}