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
```
Be sure to install the rBLAST package (https://github.com/mhahsler/rBLAST)﻿

# Example: Identifying the TEF1 gene

This example shows how to search for the TEF1 gene. Note that the EF1-728F and TEF1LLErev primers can be found on both the forward and reverse, so they are included in both the forward_primers and reverse_primers arguments.

all A logical value. If `TRUE`, the function processes all contigs. If `FALSE`, it stops after the first contig with a successful amplicon match. Defaults to `FALSE`.

You will find the results of the analysis in the specified output directory. The following files will be generated:

    amplicons_with_primers.fasta: A FASTA file containing all found amplicon sequences, including the primer regions.

    amplicons_without_primers.fasta: A FASTA file containing only the gene sequences, with the primers trimmed off.

```r
results <- analyze_trichoderma_genome(
  genome_file = "path/to/your/genome",
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

I have pre-coded all of the possible primers mentioned in the article. You can easily switch between them by simply referencing the pre-defined variables: RPB2, ITS, or TEF1.

For example, to identify the RPB2 gene, you would use this code:

```r
results <- analyze_trichoderma_genome(
  genome_file = "path/to/your/genome",
  forward_primers      = RPB2,
  reverse_primers      = RPB2,  
  max_mismatch         = 1,
  max_amplicon_length  = 5000,
  min_amplicon_length  = 100,  
  output_dir           = "RPB2_contigs_results",
  all=FALSE
)
```
# Identify 8 amplicons at once

I created a function to identify 8 amplicons at once by providing only the path to the Trichoderma genome. The function automatically creates 8 folders in your working directory, each with a different amplicon.

```r
all_amplicon_identification(genome_file = "path/to/your/genome")
```

# BLAST against TrichofindR database

The trichoderma_blast() function allows you to perform a BLAST search against a reference database defined in TrichofindR (the reference_sequences argument). Currently, the following reference databases are available:

1. ITS (ITS_reference_sequences)

2. RPB2 (RPB2_reference_sequences)

3. TEF1 (TEF1_reference_sequences)

4. LNS2 (LNS2_reference_sequences)

5. ACT (ACT_reference_sequences)

6. TUB2 (TUB_reference_sequences)

7. TEF3 (TEF3_reference_sequences)

8. PGK (PGK_reference_sequences)

```r
results_blast <- trichoderma_blast(query_sequence = "your_gene_sequence.fasta",
                               reference_sequences = ITS_reference_sequences)
```

# Automatic TEF1-RPB2-ITS identification (MATEK)

MATEK - Multilocus Analytical Trichoderma Evaluation Key

I developed an automated Trichoderma identification pipeline based on TEF1, RPB2, and ITS identification from FASTA sequence.

```r
MATEK_identification <- MATEK_identification(genome_path = "your_genome_sequence.fasta")

MATEK_identification$stage1_tef1
MATEK_identification$stage2_rpb2
MATEK_identification$stage3_its

```
# Automatic Integrative Multilocus Delimitation of Trichoderma Species (IMLDS)

The goal of this identification is to create a single FASTA file containing all amplicons and compare it against the TrichofindR database. I decided to exclude TUB2 because it always appears in two copies in Trichoderma genomes.

```r
IMLDS_identification <- IMLDS_identification(genome_path = "your_genome_sequence.fasta")

```