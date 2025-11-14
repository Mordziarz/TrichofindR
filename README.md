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

I created a function to identify 8 amplicons at once by providing only the path to the Trichoderma genome. The function automatically creates 8 folders in your working directory, each with a different gene.

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

# Automatic MIST identification

I created an automatic MIST identification system where users don't have to manually click through FASTA sequences—they simply upload their genome.

```r
MIST_identification <- MIST_identification(genome_path = "your_genome_sequence.fasta")

MIST_identification$stage1_tef1
MIST_identification$stage2_rpb2
MIST_identification$stage3_its

```
