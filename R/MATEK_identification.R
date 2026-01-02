#' Multi-Marker Identification System for Trichoderma (MATEK)
#'
#' @description
#' Performs a sequential 7-gene BLAST identification pipeline for Trichoderma species.
#' The function uses a hierarchical approach analyzing TEF1, RPB2, TEF3, PGK, ACT, 
#' ACL1, and LNS2 genomic regions to progressively narrow down species identification. 
#' Results from each marker are used to filter the reference database for the 
#' subsequent marker.
#'
#'
#' @param genome_path Character string specifying the path to the genome FASTA file.
#' @param max_mismatch Integer specifying the maximum number of mismatches allowed in primer binding.
#' @param min_amplicon_length Integer specifying the minimum amplicon length in base pairs.
#' @param max_amplicon_length Integer specifying the maximum amplicon length in base pairs.
#' @param tef1_threshold Numeric. Percent identity threshold for Stage 1. Defaults to 95.
#' @param rpb2_threshold Numeric. Percent identity threshold for Stage 2. Defaults to 95.
#' @param tef3_threshold Numeric. Percent identity threshold for Stage 3. Defaults to 95.
#' @param pgk_threshold Numeric. Percent identity threshold for Stage 4. Defaults to 95.
#' @param act_threshold Numeric. Percent identity threshold for Stage 5. Defaults to 95.
#' @param acl1_threshold Numeric. Percent identity threshold for Stage 6. Defaults to 95.
#' @param lns2_primary_threshold Numeric. Initial threshold for LNS2 (Stage 7). Defaults to 95.
#' @param lns2_secondary_threshold Numeric. Strict threshold applied if >1 matches remain in Stage 7. Defaults to 99.
#'
#' @return A list containing results from all 7 identification stages and the final identification.
#'
#' @export

MATEK_identification <- function(
    genome_path = "path_to_genome",
    max_mismatch = 1,
    min_amplicon_length = 100,
    max_amplicon_length = 5000,
    tef1_threshold = 95,
    rpb2_threshold = 95,
    tef3_threshold = 95,
    pgk_threshold = 95,
    act_threshold = 95,
    acl1_threshold = 95,
    lns2_primary_threshold = 95,
    lns2_secondary_threshold = 99) {
  
  if (!file.exists(genome_path)) {
    stop(paste("Genome file not found:", genome_path))
  }
  
  MATEK_results <- list(
    stage1_tef1 = NULL, stage2_rpb2 = NULL, stage3_tef3 = NULL,
    stage4_pgk = NULL, stage5_act = NULL, stage6_acl1 = NULL, 
    stage7_lns2 = NULL, final_identification = NA
  )
  
  run_matek_stage <- function(stage_name, primers, ref_db, filtered_ids, threshold, output_dir) {
    cat(paste0("\n[STAGE] Analyzing ", stage_name, "...\n"))
    
    analyze_trichoderma_genome(
      genome_path = genome_path, forward_primers = primers, reverse_primers = primers,
      max_mismatch = max_mismatch, max_amplicon_length = max_amplicon_length,
      min_amplicon_length = min_amplicon_length, output_dir = output_dir, all = FALSE
    )
    
    current_ref <- if (!is.null(filtered_ids)) ref_db[names(ref_db) %in% filtered_ids] else ref_db
    if (length(current_ref) == 0) {
      cat("⚠ Warning: No sequences match previous stage. Using complete database.\n")
      current_ref <- ref_db
    }
    
    res <- trichoderma_blast(
      query_sequence = paste0(output_dir, "/amplicons_with_primers.fasta"),
      reference_sequences = current_ref
    )
    
    if (nrow(res) > 0) {
      high_id <- res[res$average_pident_weighted > threshold, ]
      if (nrow(high_id) > 0) {
        res <- high_id
        cat("✓ Found", nrow(res), "matches with >", threshold, "% identity\n")
      } else {
        res <- res[1:min(5, nrow(res)), ]
        cat("⚠ No matches above", threshold, "%. Keeping top 5.\n")
      }
    }
    return(res)
  }

  MATEK_results$stage1_tef1 <- run_matek_stage("TEF1", TEF1, TEF1_reference_sequences, NULL, tef1_threshold, "TEF1_results")
  
  MATEK_results$stage2_rpb2 <- run_matek_stage("RPB2", RPB2, RPB2_reference_sequences, unique(MATEK_results$stage1_tef1$sseqid), rpb2_threshold, "RPB2_results")

  MATEK_results$stage3_tef3 <- run_matek_stage("TEF3", TEF3, TEF3_reference_sequences, unique(MATEK_results$stage2_rpb2$sseqid), tef3_threshold, "TEF3_results")

  MATEK_results$stage4_pgk <- run_matek_stage("PGK", PGK, PGK_reference_sequences, unique(MATEK_results$stage3_tef3$sseqid), pgk_threshold, "PGK_results")

  MATEK_results$stage5_act <- run_matek_stage("ACT", ACT, ACT_reference_sequences, unique(MATEK_results$stage4_pgk$sseqid), act_threshold, "ACT_results")

  MATEK_results$stage6_acl1 <- run_matek_stage("ACL1", ACL1, ACL1_reference_sequences, unique(MATEK_results$stage5_act$sseqid), acl1_threshold, "ACL1_results")

  cat("\n[STAGE 7] Analyzing LNS2 ...\n")
  analyze_trichoderma_genome(genome_path = genome_path, forward_primers = LNS2, reverse_primers = LNS2,
                             max_mismatch = max_mismatch, max_amplicon_length = max_amplicon_length,
                             min_amplicon_length = min_amplicon_length, output_dir = "LNS2_results", all = FALSE)
  
  LNS2_filtered <- LNS2_reference_sequences[names(LNS2_reference_sequences) %in% unique(MATEK_results$stage6_acl1$sseqid)]
  if (length(LNS2_filtered) == 0) LNS2_filtered <- LNS2_reference_sequences
  
  stage7_results <- trichoderma_blast(query_sequence = "LNS2_results/amplicons_with_primers.fasta", reference_sequences = LNS2_filtered)
  
  if (nrow(stage7_results) > 0) {
    high_id_lns2 <- stage7_results[stage7_results$average_pident_weighted > lns2_primary_threshold, ]
    
    if (nrow(high_id_lns2) > 0) {
      if (nrow(high_id_lns2) > 1) {
        cat("  → More than 1 genome remain. Applying strict cutoff (>", lns2_secondary_threshold, "%)...\n")
        strict_results <- high_id_lns2[high_id_lns2$average_pident_weighted > lns2_secondary_threshold, ]
        if (nrow(strict_results) > 0) {
          stage7_results <- strict_results
          cat("  → Successfully narrowed down to", nrow(stage7_results), "matches.\n")
        } else {
          cat("  → No matches at 99%. Retaining results from primary threshold.\n")
          stage7_results <- high_id_lns2
        }
      } else {
        stage7_results <- high_id_lns2
      }
    } else {
      cat("⚠ No matches above primary threshold. Keeping top 5.\n")
      stage7_results <- stage7_results[1:min(5, nrow(stage7_results)), ]
    }
  }
  MATEK_results$stage7_lns2 <- stage7_results

  cat("\n========== Final Identification ==========\n")
  if (nrow(MATEK_results$stage7_lns2) > 0) {
    MATEK_results$final_identification <- MATEK_results$stage7_lns2$sseqid
    print(MATEK_results$stage7_lns2[, c("sseqid", "average_pident_weighted")])
  } else {
    MATEK_results$final_identification <- "No reliable identification found"
  }

  return(MATEK_results)
}