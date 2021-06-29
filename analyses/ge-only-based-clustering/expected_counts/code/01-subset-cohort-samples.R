suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library("optparse"))


option_list <- list(
  make_option(c("-i", "--input_count"),
              type = "character", default = NULL,
              help = "combined expected count matrix as input"),
  make_option(c("-p","--hist_file_path"),
              type = "character", default = NULL,
              help = "Histology file with all samples from all cohorts"),
  make_option(c("-c","--out_cohort_hist_name"),
              type = "character", default = NULL,
              help = "Name of the output histology file with only desired cohort (initial_hgat) files"),
  make_option(c("-o","--out_exp_count_cohort_rds"),
              type = "character", default = NULL,
              help = "Name of of expected count file with only desired cohort (initial_hgat) samples"),
  make_option(c("-t","--out_exp_count_cohort_tsv"),
              type = "character", default = NULL,
              help = "Name of of expected count file with only desired cohort (initial_hgat) samples"))

opt <- parse_args(OptionParser(option_list = option_list))

# reading in batch corrected expected count file with all samples
count_total <- as.data.frame( readr::read_rds(opt$input_count) )

#  Reading in histology file
hist_file <- read.csv(opt$hist_file_path, sep="\t")  %>% as.data.frame()

# Filtering histology file with only RNA-seq and desired cohort (initial_hgat) samples
# Dropping molecular_subtype - without this column you can get unique and remove dup samples
cohort_rnaseq_hist <- hist_file %>%
  filter(short_histology == "HGAT") %>% filter(experimental_strategy == "RNA-Seq") %>%
  filter(tumor_descriptor == "Initial CNS Tumor") %>% 
  dplyr::select(-c(molecular_subtype))
# make it into characteros so the subsequent select would work
cohort_rnaseq_hist$Kids_First_Biospecimen_ID <- as.character(cohort_rnaseq_hist$Kids_First_Biospecimen_ID)

# Creating  a list with desired cohort cohort's BSID
cohort_BSID <- cohort_rnaseq_hist$Kids_First_Biospecimen_ID

# Subsetting expected count with only desired cohort samples and printing table to out RDS file
count_total_cohort <- count_total %>% 
  as.data.frame %>% 
  dplyr::select(cohort_BSID)
saveRDS(count_total_cohort, opt$out_exp_count_cohort_rds)
write.table(count_total_cohort, opt$out_exp_count_cohort_tsv, sep="\t")

# Write initial_hgat histology file
write.table(cohort_rnaseq_hist, opt$out_cohort_hist_name, sep="\t", row.names=FALSE)
