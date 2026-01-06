#' Multi-Locus Identification of Trichoderma Genomes (Parallel)
#'
#' @description
#' Performs high-throughput identification of Trichoderma species from genomic FASTA files. 
#' The function utilizes two identification pipelines:
#' \enumerate{
#'   \item \bold{IMLDTS}: Concatenates available oriented amplicons (up to 7 loci) into a supermatrix for BLAST.
#'   \item \bold{MATEK}: A stepwise filtration pipeline (TEF1 -> RPB2 -> ... -> LNS2) that narrows down candidates.
#' }
#' 
#' @param main_dir Character. Path to the directory containing genome FASTA files (.fna or .fasta).
#' @param n_cores Integer. Number of CPU cores to use for parallel processing (default 56).
#'
#' @details 
#' The function enforces strict sequence orientation (Forward-to-Reverse) for all extracted loci.
#' It applies a 95\% identity threshold for the MATEK stages and outputs a cleaned final 
#' identification string without percentage values.
#'
#' @import parallel
#' @import Biostrings
#' @import dplyr
#' @import TrichofindR
#' @export

multi_identification <- function(main_dir="path/to/directory",n_cores = 56) {

  genome_files <- list.files(path = main_dir, pattern = "\\.(fna|fasta)$", full.names = TRUE)
  
  n_genomes <- length(genome_files)
  if (n_genomes == 0) {
    message("⚠️ No genome files found in: ", main_dir)
    return(NULL)
  }
  
  message("🚀 Starting analysis. Genomes found: ", n_genomes)
  message("⚙️ Using ", n_cores, " CPU cores.")
  
  imldts_loci <- c("TEF1", "RPB2", "TEF3", "LNS2", "ACT", "PGK", "ACL1")

for(d in c(paste0(imldts_loci, "_results_all"), "IMLDTS_ultrafastas_all")) {
  dir.create(file.path(main_dir, d), showWarnings = FALSE, recursive = TRUE)
}

format_hit <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("No hit")
  all_hits <- paste0(df$sseqid, " (", round(df$average_pident_weighted, 2), "%)")
  return(paste(all_hits, collapse = ", "))
}

process_genome_full <- function(f_path) {
  g_name <- gsub("\\.(fna|fasta)$", "", basename(f_path))
  g_dir <- file.path(main_dir, paste0("results_", g_name))
  if (!dir.exists(g_dir)) dir.create(g_dir, recursive = TRUE)
  
  old_wd <- getwd()
  setwd(g_dir)
  
  tryCatch({
    all_amplicon_identification(genome_path = f_path)
    found_sequences_oriented <- list()
    for (l in imldts_loci) {
      f_amp <- file.path(g_dir, paste0(l, "_results"), "amplicons_with_primers.fasta")
      if (file.exists(f_amp)) {
        seqs <- readDNAStringSet(f_amp)
        if (length(seqs) > 0) {
          seq_dna <- seqs[[1]]; seq_rc <- reverseComplement(seq_dna)
          fwd_p <- get(paste0(l, "_F")); rev_p <- get(paste0(l, "_R"))
          if (!(find_primer_smart(seq_dna, fwd_p, position = "start") || find_primer_smart(seq_dna, rev_p, position = "end"))) {
            if (find_primer_smart(seq_rc, fwd_p, position = "start")) seq_dna <- seq_rc
          }
          out_s <- DNAStringSet(seq_dna); names(out_s) <- g_name
          writeXStringSet(out_s, file.path(main_dir, paste0(l, "_results_all"), paste0(g_name, ".fasta")))
          p_clean <- file.path(g_dir, paste0(l, "_results"), "amplicons_without_primers.fasta")
          if (file.exists(p_clean)) {
            s_clean <- readDNAStringSet(p_clean)[[1]]
            if (identical(seq_dna, seq_rc)) s_clean <- reverseComplement(s_clean)
            found_sequences_oriented[[l]] <- as.character(s_clean)
          }
        }
      }
    }
    imldts_res <- "No loci found"
    if (length(found_sequences_oriented) > 0) {
      available_loci <- intersect(imldts_loci, names(found_sequences_oriented))
      concat_str <- paste(found_sequences_oriented[available_loci], collapse = "")
      
      ultra_ds <- DNAStringSet(concat_str)
      names(ultra_ds) <- paste0(g_name, "_IMLDTS_loci_", length(available_loci), "_of_7")
      
      writeXStringSet(ultra_ds, file.path(main_dir, "IMLDTS_ultrafastas_all", paste0(g_name, "_IMLDTS.fasta")))
      
      tmp_f <- tempfile(fileext = ".fasta")
      writeXStringSet(ultra_ds, tmp_f)
      imldts_res <- format_hit(trichoderma_blast(tmp_f, IMLDTS_reference_sequences))
    }
    matek_stages <- list(
      list(name="TEF1", db=TEF1_reference_sequences),
      list(name="RPB2", db=RPB2_reference_sequences),
      list(name="TEF3", db=TEF3_reference_sequences),
      list(name="PGK",  db=PGK_reference_sequences),
      list(name="ACT",  db=ACT_reference_sequences),
      list(name="ACL1", db=ACL1_reference_sequences),
      list(name="LNS2", db=LNS2_reference_sequences)
    )
    
    m_results <- list(); current_candidates <- NULL; last_blast_tables <- list() 
    
    for (s in matek_stages) {
      f_amp_m <- file.path(g_dir, paste0(s$name, "_results"), "amplicons_with_primers.fasta")
      if (file.exists(f_amp_m)) {
        ref_db <- s$db
        if (!is.null(current_candidates)) {
          filt_db <- ref_db[names(ref_db) %in% current_candidates]
          if (length(filt_db) > 0) ref_db <- filt_db
        }
        
        res_m <- trichoderma_blast(f_amp_m, ref_db)
        res_m <- res_m[res_m$average_pident_weighted >= 95, ] # Filtr 95%
        
        if (nrow(res_m) > 0) {
          last_blast_tables[[s$name]] <- res_m 
          m_results[[paste0("Stage_", s$name)]] <- format_hit(res_m)
          current_candidates <- unique(res_m$sseqid)
        } else {
          m_results[[paste0("Stage_", s$name)]] <- "No hit above 95%"
          current_candidates <- NULL 
        }
      } else {
        m_results[[paste0("Stage_", s$name)]] <- "Missing"
      }
    }
    final_id <- "No reliable identification"
    l_tab <- last_blast_tables[["LNS2"]]; a_tab <- last_blast_tables[["ACL1"]]
    
    raw_final <- ""
    if (!is.null(l_tab) && nrow(l_tab) > 0) {
      s_lns2 <- l_tab[l_tab$average_pident_weighted > 99, ]
      raw_final <- if(nrow(s_lns2) > 0) format_hit(s_lns2) else format_hit(l_tab)
    } else if (!is.null(a_tab) && nrow(a_tab) > 0) {
      raw_final <- format_hit(a_tab)
    } else {
      v_st <- names(m_results)[sapply(m_results, function(x) !x %in% c("No hit above 95%", "Missing"))]
      if(length(v_st) > 0) raw_final <- m_results[[tail(v_st, 1)]]
    }
    
    if (raw_final != "") { final_id <- gsub(" \\(.*?%\\)", "", raw_final) }
    
    final_row <- data.frame(
      Genome = g_name,
      IMLDTS_results = imldts_res,
      MATEK_final_identification = final_id,
      as.data.frame(m_results),
      Status = "Success",
      stringsAsFactors = FALSE
    )
    setwd(old_wd)
    return(final_row)
    
  }, error = function(e) {
    if (exists("old_wd")) setwd(old_wd)
    return(data.frame(Genome = g_name, Status = paste("Error:", e$message), stringsAsFactors = FALSE))
  })
}
results_list <- mclapply(genome_files, process_genome_full, mc.cores = n_cores)
final_table <- bind_rows(results_list)
write.csv(final_table, file.path(main_dir, "Summary_Report_Parallel.csv"), row.names = FALSE)

for(l in imldts_loci) {
  files <- list.files(file.path(main_dir, paste0(l, "_results_all")), full.names = TRUE)
  if(length(files) > 0) {
    combined <- do.call(c, lapply(files, readDNAStringSet))
    writeXStringSet(combined, file.path(main_dir, paste0("Combined_", l, "_all_genomes.fasta")))
  }
}
message("✨ Everything is ready.")
return(final_table)
}