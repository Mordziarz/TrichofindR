#' BLAST Search Against Trichoderma Reference Sequences
#' 
#' @description
#' Performs a BLAST nucleotide search of query DNA sequences against a reference
#' database. Optimized for high-throughput genomic identification.
#'
#' @param query_sequence Path to FASTA file.
#' @param reference_sequences DNAStringSet or list of sequences.
#'
#' @return A data frame with aggregated BLAST results.
#' @export

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

  dna_sequences <- Biostrings::DNAStringSet(sequences_vector)

  temp_id <- sample(1000:9999, 1)
  subject_db_path <- paste0("temp_subject_db_", temp_id, ".fasta")
  Biostrings::writeXStringSet(dna_sequences, filepath = subject_db_path)

  on.exit({
    files_to_remove <- list.files(pattern = paste0("temp_subject_db_", temp_id))
    if(length(files_to_remove) > 0) file.remove(files_to_remove)
  })

  rBLAST::makeblastdb(subject_db_path, dbtype = "nucl")
  blast_db <- rBLAST::blast(db = subject_db_path, type = "blastn")

  query_sequences <- Biostrings::readDNAStringSet(query_sequence)
  
  if (length(query_sequences) == 0) {
    stop("No valid sequences found in query file")
  }

  blast_result <- rBLAST::predict(blast_db, query_sequences)

  if (is.null(blast_result) || nrow(blast_result) == 0) {
    message("--- No BLAST hits found for: ", basename(query_sequence), " ---")
    return(data.frame()) 
  }

  results_aggregated <- blast_result %>%
    dplyr::mutate(matching_nts = (pident * length) / 100) %>%
    dplyr::group_by(qseqid, sseqid) %>%
    dplyr::summarise(
      total_length = sum(length, na.rm = TRUE),                 
      total_matching_nts = sum(matching_nts, na.rm = TRUE),    
      best_evalue = min(evalue, na.rm = TRUE),
      max_bitscore = max(bitscore, na.rm = TRUE),
      .groups = 'drop' 
    ) %>%
    dplyr::mutate(
      average_pident_weighted = (total_matching_nts / total_length) * 100
    ) %>%
    dplyr::filter(total_length > 0) %>%
    dplyr::arrange(desc(average_pident_weighted), desc(total_length))

  return(as.data.frame(results_aggregated))
}