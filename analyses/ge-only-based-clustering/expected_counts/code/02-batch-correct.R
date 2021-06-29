# Author: Komal S. Rathi
# Function: Batch correction using sva::ComBat and sva::ComBat_seq

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
output_dir <- file.path(root_dir, "analyses", "ge-only-based-clustering", "expected_counts", "input")

# parameters
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "cohort filtered expression matrices path to be used (TPM, FPKM or expected counts) (.rds)"),
  make_option(c("--metadata"), type = "character",
              help = "cohort filtered histologies file path to be used (.tsv)"),
  make_option(c("--output_prefix"), type = "character",
              help = "prefix for output files"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
mat <- opt$mat
metadata <- opt$metadata
output_prefix <- opt$output_prefix

# output files
uncorrected_outfile <- file.path(output_dir, paste0(output_prefix, '_uncorrected.rds'))
corrected_outfile <- file.path(output_dir, paste0(output_prefix, '_corrected.rds'))

# read in cohort filtered expression and histologies file
metadata <- read.delim(metadata, header=T, sep = "\t", stringsAsFactors = F) %>%
  column_to_rownames("Kids_First_Biospecimen_ID")
uncorrected_mat <- readRDS(mat)

# take an intersection if metadata and expression have different number of samples
print(paste0("Samples to use:", nrow(metadata)))

# double-check if everything is lined-up between expression matrix and metadata
if(identical(rownames(metadata), colnames(uncorrected_mat))){
  print("Matching dimensions")
  print("Save uncorrected file...")
  saveRDS(uncorrected_mat, file = uncorrected_outfile)
} else {
  print("Check inputs")
  break
}

# batch correct using ComBat
corrected_mat <- ComBat_seq(counts = as.matrix(uncorrected_mat), batch = metadata$RNA_library)

print("Save corrected file...")
saveRDS(corrected_mat, file = corrected_outfile)