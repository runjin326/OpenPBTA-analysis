# Author: Komal S. Rathi
# Function: Filter matrix to protein coding genes

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--input"), type = "character",
              help = "Input expression matrix (.rds)"),
  make_option(c("--gencode_version"), type = "character",
              help = "Gencode version"),
  make_option(c("--output"), type = "character",
              help = "Output file with path (.rds)")
)
opt <- parse_args(OptionParser(option_list = option_list))
input <- opt$input
gencode_version <- opt$gencode_version
output <- opt$output

# read gtf and filter to protein coding 
fname <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_', gencode_version, '/gencode.v', gencode_version, '.primary_assembly.annotation.gtf.gz')
gencode_gtf <- rtracklayer::import(con = fname)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter input matrix 
input <- readRDS(input)
input_pc <- input[rownames(input) %in% gencode_gtf$gene_name,]
saveRDS(input_pc, file = output)

