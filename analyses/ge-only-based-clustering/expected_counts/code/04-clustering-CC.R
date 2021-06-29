
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ProteoMM"))
suppressPackageStartupMessages(library("ConsensusClusterPlus"))
suppressPackageStartupMessages(library(edgeR))

source("code/utils/UQ.R")

option_list <- list(
  make_option(c("-c", "--cohort_id"),
    type = "character", default = NULL,
    help = "cohort id "
  ),
  make_option(c("-p", "--patient_id_column"),
    type = "character", default = NULL,
    help = "patient id column id "
  ),
  make_option(c("-o", "--output_folder"),
    type = "character", default = NULL,
    help = " output folder path "
  ),
  make_option(c("-l", "--cohort_label_file"),
    type = "character", default = NULL,
    help = "main manifest for cohort "
  ),
  make_option(c("-e", "--gene_expression_file"),
    type = "character", default = "initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds",
    help = " RDS file with Samples as columns and Genes as Rows"
  ),
  make_option(c("-d", "--distance"),
    type = "character", default = "spearman"
  ),
  make_option(c("-f", "--finalLinkage"),
    type = "character", default = "average"
  ),
  make_option(c("-i", "--innerLinkage"),
    type = "character", default = "average"
  ),
  make_option(c("-a", "--clusterAlg"),
    type = "character", default = "hc"
  ),
  make_option(c("-k", "--maxK"),
    type = "integer", default = "8"
  ),
  make_option(c("-r", "--reps"),
    type = "integer", default = "100"
  ),
  make_option(c("--pItem"),
    type = "double", default = "0.8"
  ),
  make_option(c("--transform"),
              type = "character", default = NULL
  )
)

# set root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


# read opt args
opt <- parse_args(OptionParser(option_list = option_list))
cohort_id <- opt$cohort_id
patient_id_column <- opt$patient_id_column
cohort_label_file <- opt$cohort_label_file
includes_progressive <- opt$includes_progressive
distance <- opt$distance
innerLinkage <- opt$innerLinkage
finalLinkage <- opt$finalLinkage
clusterAlg <- opt$clusterAlg
output_folder <- opt$output_folder
gene_expression_file <- opt$gene_expression_file
transform <- opt$transform

# setting scientific notation to 999
options(scipen=999)

# read cohort specific histology file
cohort_plot_labels <- read_tsv(file.path(
  root_dir,
  "data",
  cohort_label_file
)) %>%
  as.data.frame()

### The input RNAseq dataset
GeneExp <- readRDS(file.path(root_dir, "data", gene_expression_file)) %>%
  as.data.frame() %>%
  dplyr::select(intersect(colnames(.), cohort_plot_labels$Kids_First_Biospecimen_ID))

## if transformation is required using UQ
if (!is.null(transform) && transform == "uq.pgQ2"){
  GeneExp <- uq.pgQ2(as.matrix(GeneExp))
}


## if transformation is required using edgeR calcnorm
if (!is.null(transform) && transform == "edgeR"){
  GeneExp_dgelist <- DGEList(counts=GeneExp)
  GeneExp_dgelist <- estimateCommonDisp(GeneExp_dgelist)
  GeneExp_dgelist <- calcNormFactors(GeneExp_dgelist, method="TMM")
  logCPM <- cpm(GeneExp_dgelist, log=TRUE)
  keep<-filterByExpr(logCPM,group = NULL, design = NULL)
  GeneExp <- logCPM[keep,]
}


## if transformation is required using rlog from DESEq
#if (!is.null(transform) && transform == "rlog"){
#  GeneExp <- GeneExp %>% mutate(across(where(is.numeric), ~ round(., digits = 0)))
#  GeneExp <- rlog(as.matrix(GeneExp), blind = TRUE)
#}

## if transformation is required using vst
if (!is.null(transform) && transform == "vst"){
  # corrected TPM files is a double type matrix (decimal points)
  # However, the VST needs integer values or errors out saying
  # Error in DESeqDataSet(se, design = design, ignoreRank) :
  #    some values in assay are not integers
  GeneExp <- varianceStabilizingTransformation(round(as.matrix(GeneExp)),blind = TRUE,fitType="parametric")
  }

# Log transforming if not vst
if (!is.null(transform) && transform == "log"){
 GeneExp <- log2(GeneExp)
}

# get most variable
means <- rowMeans(GeneExp)
vars <- apply(GeneExp, 1, var)
varorder <- order(vars, decreasing = T)
GeneExp <- GeneExp[varorder, ] %>% head(5000)

# tumorData mean and sd
tumorData_means <- rowMeans(GeneExp, na.rm = TRUE)
tumorData_sd <- apply(GeneExp, 1, sd, na.rm = TRUE)
# subtract mean
GeneExpzscored <- sweep(GeneExp, 1, tumorData_means, FUN = "-")
# divide by SD remove NAs and Inf values from zscore for genes with 0
GeneExpzscored <- sweep(GeneExpzscored, 1, tumorData_sd, FUN = "/") %>%
  na_if(Inf) %>%
  na.omit()

# rename colnames to Participant ID
colnames_index <- unlist(
  lapply(colnames(GeneExpzscored), function(x) {
    grep(x, cohort_plot_labels$Kids_First_Biospecimen_ID)
  })
)
colnames(GeneExpzscored) <- cohort_plot_labels[colnames_index, patient_id_column]

if (distance == "manhattan"){
 manhattan = function(x){ dist(x,method="manhattan")}
}

# moving to results directory per cohort_id
# if missing create folder
if (!file.exists(file.path(output_folder))) {
  dir.create(file.path(output_folder), recursive = T)
}

# set working directory
setwd(file.path(output_folder))

# save params list in each folder
params <- list(
  "cohort_id" = opt$cohort_id,
  "patient_id_column" = opt$patient_id_column,
  "cohort_label_file" = opt$cohort_label_file,
  "includes_progressive" = opt$includes_progressive,
  "distance" = opt$distance,
  "innerLinkage" = opt$innerLinkage,
  "finalLinkage" = opt$finalLinkage,
  "clusterAlg" = opt$clusterAlg,
  "output_folder" = opt$output_folder
)
list.save(params, "params_list.yaml")

## Consensus cluster plot per dataset
result_CC <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(GeneExpzscored),
  finalLinkage = finalLinkage,
  distance = distance,
  clusterAlg = clusterAlg,
  plot = "pdf",
  reps = 100, maxK = 10, pItem = 0.8,
  title = paste0(cohort_id),
  seed = 123
)

save(result_CC, file = file.path(
  paste0(cohort_id),
  "CC.Rdata"
))

saveRDS(GeneExpzscored, file = file.path(
  "GeneExpzscored.rds"
))

# rename colnames to Participant ID
colnames_index <- unlist(
  lapply(colnames(GeneExp), function(x) {
    grep(x, cohort_plot_labels$Kids_First_Biospecimen_ID)
  })
)
colnames(GeneExp) <- cohort_plot_labels[colnames_index, patient_id_column]


saveRDS(GeneExp,file = file.path(
  "GeneExp.rds"
))
