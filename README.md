# TrichofindR
 Tool for identifying Trichoderma species.

# Installation

```r
install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/TrichofindR')
library(TrichofindR)
```
# Libraries 

```r
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(IRanges)
library(rBLAST)
library(dplyr)
library(parallel)
```
Be sure to install the rBLAST package (https://github.com/mhahsler/rBLAST)﻿

# Example: Identifying the TEF1 gene

This example shows how to search for the TEF1 gene. Note that the EF1-728F and TEF1LLErev primers can be found on both the forward and reverse, so they are included in both the forward_primers and reverse_primers arguments.

`all` a logical value. If `TRUE`, the function processes all contigs. If `FALSE`, it stops after the first contig with a successful amplicon match. Defaults to `FALSE`.

You will find the results of the analysis in the specified output directory. The following files will be generated:

    amplicons_with_primers.fasta: A FASTA file containing all found amplicon sequences, including the primer regions.

    amplicons_without_primers.fasta: A FASTA file containing only the gene sequences, with the primers trimmed off.

```r
results <- analyze_trichoderma_genome(
  genome_path = "path/to/your/genome",
  forward_primers      = c("CATCGAGAAGTTCGAGAAGG","AACTTGCAGGCAATGTGG"),
  reverse_primers      = c("AACTTGCAGGCAATGTGG","CATCGAGAAGTTCGAGAAGG"),  
  max_mismatch         = 1,
  max_amplicon_length  = 5000,
  min_amplicon_length  = 100,  
  output_dir           = "tef1_contigs_results",
  all=FALSE
)
```
# Using pre-defined primer sets

I have pre-coded all of the possible primers mentioned in the article. You can easily switch between them by simply referencing the pre-defined variables: RPB2, ITS, TEF1, LNS2, ACT, TUB2, TEF3, ACL1 or PGK.

For example, to identify the RPB2 gene, you would use this code:

```r
results <- analyze_trichoderma_genome(
  genome_path = "path/to/your/genome",
  forward_primers      = RPB2,
  reverse_primers      = RPB2,  
  max_mismatch         = 1,
  max_amplicon_length  = 5000,
  min_amplicon_length  = 100,  
  output_dir           = "RPB2_contigs_results",
  all=FALSE
)
```
# Identify 9 amplicons at once

I created a function to identify 9 amplicons at once by providing only the path to the Trichoderma genome. The function automatically creates 9 folders in your working directory, each with a different amplicon.

```r
all_amplicon_identification(genome_path = "path/to/your/genome")
```

# BLAST against TrichofindR database

The trichoderma_blast() function allows you to perform a BLAST search against a reference database defined in TrichofindR (the reference_sequences argument). Currently, the following reference databases are available:

1. ITS (ITS_reference_sequences)

2. RPB2 (RPB2_reference_sequences)

3. TEF1 (TEF1_reference_sequences)

4. LNS2 (LNS2_reference_sequences)

5. ACT (ACT_reference_sequences)

6. TUB2 (TUB2_reference_sequences)

7. TEF3 (TEF3_reference_sequences)

8. PGK (PGK_reference_sequences)

9. ACL1 (ACL1_reference_sequences)

10. Ultra_fasta (IMLDTS_reference_sequences)

I've implemented a critical change to how I process BLAST results within the trichoderma_blast function, shifting from a simple percent identity (pident) to a Length-Weighted Average Percent Identity. This modification was necessary because the standard BLAST output often fragments a single sequence match into multiple High-Scoring Segment Pairs (HSPs), meaning a simple arithmetic average of the % identities would be misleading. My new logic works in three steps: I first calculate the absolute number of matching nucleotides for each fragment (matching_nts = pident * length / 100), which assigns a length-based weight. Then, I aggregate these results by species (sseqid), summing both the matching nucleotides and the total covered length for all fragments belonging to that species. Finally, I calculate the true weighted identity (average_pident_weighted = (total_matching_nts / total_length) * 100) by dividing the total matching nucleotides by the total covered length. This approach ensures that longer, more complete, and therefore more significant fragments have a proportionally greater impact on the final score, providing a more robust and biologically meaningful measure of similarity across the entire covered region.

```r
results_blast <- trichoderma_blast(query_sequence = "your_gene_sequence.fasta",
                               reference_sequences = ITS_reference_sequences)
```

# MATEK — Multilocus Analytical Trichoderma Evaluation Key

MATEK - Multilocus Analytical Trichoderma Evaluation Key

I developed an automated Trichoderma identification pipeline based on TEF1, RPB2, TEF3, PGK, ACT, ACL1, and LNS2 identification from FASTA sequence.

```r
MATEK_identification <- MATEK_identification(genome_path = "your_genome_sequence.fasta")

MATEK_identification$stage1_tef1
MATEK_identification$stage2_rpb2
MATEK_identification$stage3_tef3
MATEK_identification$stage4_pgk
MATEK_identification$stage5_act
MATEK_identification$stage6_acl1
MATEK_identification$stage7_lns2

MATEK_identification$final_identification
```
# Automatic Integrative Multilocus Delimitation of Trichoderma Species (IMLDTS)

The goal of this identification is to create a single FASTA file containing all amplicons and compare it against the TrichofindR database. I decided to exclude TUB2 because it always appears in two copies in Trichoderma genomes. ITS was also excluded from this analysis.

```r
IMLDTS_identification <- IMLDTS_identification(genome_path = "your_genome_sequence.fasta",identity_threshold = 95,max_target_seqs = 10)

```

# Batch Barcode Extraction and Species Identification

If you have a folder full of genome files (FASTA/FNA) and you want to analyze them all at once, the TrichofindR::multi_identification() function is your go-to tool.

It automatically scans your folder, finds the barcodes, and runs the identification pipelines (IMLDTS and MATEK) for every single genome it finds. Once it’s finished, you’ll get a clean summary table called Summary_Report_Parallel.csv with all the results.

What does it actually do?

1. Parallel Processing: Uses multiple CPU cores to speed things up (perfect for large datasets).

2. Automatic Cleanup: Organizes results into separate folders for each genome.

3. Orientation Fix: It makes sure every sequence is saved in the correct Forward-to-Reverse orientation.

4. Combined Files: It automatically creates "Combined" FASTA files for each locus (e.g., Combined_ACL1_all_genomes.fasta), so you have all your sequences ready for further analysis like building trees.

```r
multi <- TrichofindR::multi_identification(main_dir = "path/to/your/folder",n_cores = 50)
```