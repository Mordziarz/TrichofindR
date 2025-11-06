#' Multi-Marker Identification System for Trichoderma (MIST)
#'
#' @description
#' Performs a sequential multi-gene BLAST identification pipeline for Trichoderma species.
#' The function uses a hierarchical approach analyzing TEF1, RPB2, and ITS genomic regions
#' to progressively narrow down species identification. Results from each marker are used
#' to filter the reference database for the subsequent marker, improving specificity and
#' reducing computational burden.
#'
#' @details
#' The MIST (Multi-marker Identification System for Trichoderma) function implements a
#' three-stage identification workflow:
#' \enumerate{
#'   \item \strong{Stage 1 - TEF1}: Analyzes the translation elongation factor 1-alpha gene.
#'     Filters results for >95% identity. If no matches found, retains top 5 best matches.
#'
#'   \item \strong{Stage 2 - RPB2}: Analyzes the RNA polymerase II largest subunit gene.
#'     Uses TEF1-filtered reference database to narrow down candidates. Applies same
#'     filtering criteria as Stage 1.
#'
#'   \item \strong{Stage 3 - ITS}: Analyzes the Internal Transcribed Spacer region.
#'     Uses RPB2-filtered reference database. Applies adaptive filtering: initially filters
#'     for >95% identity. If more than 2 results are found, applies stricter 99% threshold.
#'     If 2 or fewer results, keeps the 95% matches for broader identification.
#' }
#'
#' This hierarchical approach combines the discriminatory power of multiple genes while
#' maintaining computational efficiency through progressive database reduction.
#'
#' @param query_sequence Character string specifying the path to a genome FASTA file containing
#'   the complete or partial Trichoderma genome to be identified. Required for primer binding
#'   and amplicon extraction. Should be in FASTA format with sequences representing the target
#'   organism's genomic DNA.
#'
#' @param genome_path Character string specifying the path to the genome file. Defaults to
#'   \code{"~/Pobrane/24THUB.fasta"}. Can be modified to analyze different genome files.
#'
#' @param max_mismatch Integer specifying the maximum number of mismatches allowed in primer
#'   binding sites. Defaults to \code{1}. Controls primer specificity.
#'
#' @param min_amplicon_length Integer specifying the minimum amplicon length in base pairs.
#'   Defaults to \code{100}.
#'
#' @param max_amplicon_length Integer specifying the maximum amplicon length in base pairs.
#'   Defaults to \code{5000}.
#'
#' @param tef1_threshold Numeric specifying the percent identity threshold for TEF1 results.
#'   Defaults to \code{95}. Matches above threshold are retained; if none found, top 5 best
#'   matches are kept.
#'
#' @param rpb2_threshold Numeric specifying the percent identity threshold for RPB2 results.
#'   Defaults to \code{95}. Same filtering logic as TEF1.
#'
#' @param its_primary_threshold Numeric specifying the initial percent identity threshold for ITS results.
#'   Defaults to \code{95}. Used for primary filtering in Stage 3.
#'
#' @param its_secondary_threshold Numeric specifying the stricter percent identity threshold applied
#'   when more than 2 results are found at primary threshold. Defaults to \code{99}. Provides
#'   high-confidence identifications when multiple candidates remain.
#'
#' @return
#' A list containing results from all three identification stages:
#' \describe{
#'   \item{stage1_tef1}{Data frame of BLAST results from TEF1 analysis}
#'   \item{stage2_rpb2}{Data frame of BLAST results from RPB2 analysis}
#'   \item{stage3_its}{Data frame of BLAST results from ITS analysis}
#'   \item{final_identification}{Character vector with the most likely species identification}
#'   \item{stage1_reference_ids}{Character vector of reference sequence IDs passing TEF1 filter}
#'   \item{stage2_reference_ids}{Character vector of reference sequence IDs passing RPB2 filter}
#' }
#'
#' @examples
#' \dontrun{
#' # Run MIST identification pipeline on a genome file
#' mist_results <- MIST_identification(
#'   query_sequence = "~/Pobrane/24THUB.fasta"
#' )
#'
#' # View results from each stage
#' head(mist_results$stage1_tef1, 10)
#' head(mist_results$stage2_rpb2, 10)
#' head(mist_results$stage3_its, 10)
#'
#' # Get final species identification
#' mist_results$final_identification
#' }
#'
#' @seealso
#' \code{\link{analyze_trichoderma_genome}} for primer-based amplicon extraction
#' \code{\link{trichoderma_blast}} for BLAST sequence similarity searches
#'
#' @importFrom rBLAST predict
#' @export
#'
#' @author Mateusz Mazdziarz
#' @keywords fungal-identification multi-marker Trichoderma TEF1 RPB2 ITS phylogenetics


MIST_identification <- function(
    genome_path = "path_to_genome",
    max_mismatch = 1,
    min_amplicon_length = 100,
    max_amplicon_length = 5000,
    tef1_threshold = 95,
    rpb2_threshold = 95,
    its_primary_threshold = 95,
    its_secondary_threshold = 99) {
  
  if (!file.exists(genome_path)) {
    stop(paste("Genome file not found:", genome_path))
  }
  
  mist_results <- list(
    stage1_tef1 = NULL,
    stage2_rpb2 = NULL,
    stage3_its = NULL,
    stage1_reference_ids = NULL,
    stage2_reference_ids = NULL,
    final_identification = NA
  )
  
  cat("\n========== MIST Identification Pipeline ==========\n")
  cat("\n[STAGE 1] Analyzing TEF1 (Translation Elongation Factor 1-alpha)...\n")
  
  analyze_trichoderma_genome(
    genome_file = genome_path,
    forward_primers = TEF1,
    reverse_primers = TEF1,
    max_mismatch = max_mismatch,
    max_amplicon_length = max_amplicon_length,
    min_amplicon_length = min_amplicon_length,
    output_dir = "TEF1_results",
    all = FALSE
  )


  stage1_results <- trichoderma_blast(
    query_sequence = "TEF1_results/amplicons_with_primers.fasta",
    reference_sequences = TEF1_reference_sequences
  )
  
  if (nrow(stage1_results) > 0) {
    high_identity_tef1 <- stage1_results[stage1_results$pident > tef1_threshold, ]
    
    if (nrow(high_identity_tef1) > 0) {
      stage1_results <- high_identity_tef1
      cat("✓ Found", nrow(stage1_results), "matches with >", tef1_threshold, "% identity\n")
    } else {
      stage1_results <- stage1_results[1:min(5, nrow(stage1_results)), ]
      cat("⚠ No matches above", tef1_threshold, "% identity. Keeping top 5 best matches.\n")
    }
  }
  
  mist_results$stage1_tef1 <- stage1_results
  
  stage1_reference_ids <- unique(stage1_results$sseqid)
  mist_results$stage1_reference_ids <- stage1_reference_ids
  
  cat("\n--- Stage 1 Results (TEF1) ---\n")
  print(stage1_results[, c("sseqid", "pident", "length")])
  
  cat("\n[STAGE 2] Analyzing RPB2 (RNA Polymerase II Largest Subunit)...\n")
  
  analyze_trichoderma_genome(
    genome_file = genome_path,
    forward_primers = RPB2,
    reverse_primers = RPB2,
    max_mismatch = max_mismatch,
    max_amplicon_length = max_amplicon_length,
    min_amplicon_length = min_amplicon_length,
    output_dir = "RPB2_results/",
    all = FALSE
  )
  
  RPB2_filtered <- RPB2_reference_sequences[names(RPB2_reference_sequences) %in% stage1_reference_ids]
  
  if (length(RPB2_filtered) == 0) {
    cat("⚠ Warning: No RPB2 sequences match TEF1 stage results. Using complete RPB2 database.\n")
    RPB2_filtered <- RPB2_reference_sequences
  }
  
  stage2_results <- trichoderma_blast(
    query_sequence = "RPB2_results/amplicons_with_primers.fasta",
    reference_sequences = RPB2_filtered
  )
  
  if (nrow(stage2_results) > 0) {
    high_identity_rpb2 <- stage2_results[stage2_results$pident > rpb2_threshold, ]
    
    if (nrow(high_identity_rpb2) > 0) {
      stage2_results <- high_identity_rpb2
      cat("✓ Found", nrow(stage2_results), "matches with >", rpb2_threshold, "% identity\n")
    } else {
      stage2_results <- stage2_results[1:min(5, nrow(stage2_results)), ]
      cat("⚠ No matches above", rpb2_threshold, "% identity. Keeping top 5 best matches.\n")
    }
  }
  
  mist_results$stage2_rpb2 <- stage2_results
  
  stage2_reference_ids <- unique(stage2_results$sseqid)
  mist_results$stage2_reference_ids <- stage2_reference_ids
  
  cat("\n--- Stage 2 Results (RPB2) ---\n")
  print(stage2_results[, c("sseqid", "pident", "length")])
  
  cat("\n[STAGE 3] Analyzing ITS (Internal Transcribed Spacer Region)...\n")
  
  analyze_trichoderma_genome(
    genome_file = genome_path,
    forward_primers = ITS,
    reverse_primers = ITS,
    max_mismatch = max_mismatch,
    max_amplicon_length = max_amplicon_length,
    min_amplicon_length = min_amplicon_length,
    output_dir = "ITS_results/",
    all = FALSE
  )
  
  ITS_filtered <- ITS_reference_sequences[names(ITS_reference_sequences) %in% stage2_reference_ids]
  
  if (length(ITS_filtered) == 0) {
    cat("⚠ Warning: No ITS sequences match RPB2 stage results. Using complete ITS database.\n")
    ITS_filtered <- ITS_reference_sequences
  }
  
  stage3_results <- trichoderma_blast(
    query_sequence = "ITS_results/amplicons_with_primers.fasta",
    reference_sequences = ITS_filtered
  )
  
  if (nrow(stage3_results) > 0) {
    high_identity_its <- stage3_results[stage3_results$pident > its_primary_threshold, ]
    
    if (nrow(high_identity_its) > 0) {
      cat("✓ Found", nrow(high_identity_its), "matches with >", its_primary_threshold, "% identity\n")
      
      if (nrow(high_identity_its) > 2) {
        cat("  → More than 2 results found. Applying stricter threshold (>", its_secondary_threshold, "%)...\n")
        stage3_results <- high_identity_its[high_identity_its$pident > its_secondary_threshold, ]
        
        if (nrow(stage3_results) == 0) {
          cat("  → No results at", its_secondary_threshold, "% threshold. Keeping", nrow(high_identity_its), "matches at", its_primary_threshold, "%.\n")
          stage3_results <- high_identity_its
        } else {
          cat("  → Retained", nrow(stage3_results), "matches at", its_secondary_threshold, "% identity\n")
        }
      } else {
        cat("  → ", nrow(high_identity_its), "result(s) found. Keeping at", its_primary_threshold, "% threshold.\n")
        stage3_results <- high_identity_its
      }
    } else {
      cat("⚠ No matches above", its_primary_threshold, "% identity. Keeping top 5 best matches.\n")
      stage3_results <- stage3_results[1:min(5, nrow(stage3_results)), ]
    }
  }
  
  mist_results$stage3_its <- stage3_results

  cat("\n--- Stage 3 Results (ITS) ---\n")
  print(stage3_results[, c("sseqid", "pident", "length")])
  
  cat("\n========== Final Identification ==========\n")

  if (nrow(stage3_results) > 0) {
    high_confidence_results <- stage3_results[stage3_results$pident >= 99, ]
    
    if (nrow(high_confidence_results) > 0) {
      cat("High confidence matches:\n\n")
      
      for (i in seq_len(nrow(high_confidence_results))) {
        match <- high_confidence_results[i, ]
        cat("  ", i, ". ", match$sseqid, " (", match$pident, "% identity)\n")
      }
      
      mist_results$final_identification <- high_confidence_results
      
    } else {
      cat("Probable matches (≥", its_primary_threshold, "% identity):\n\n")
      
      for (i in seq_len(nrow(stage3_results))) {
        match <- stage3_results[i, ]
        cat("  ", i, ". ", match$sseqid, " (", match$pident, "% identity)\n")
      }
      
      mist_results$final_identification <- stage3_results
    }
    
  } else {
    mist_results$final_identification <- "No reliable Trichoderma identification found"
    cat(mist_results$final_identification, "\n")
  }

  cat("\n==================================================\n\n")

  return(mist_results)
}