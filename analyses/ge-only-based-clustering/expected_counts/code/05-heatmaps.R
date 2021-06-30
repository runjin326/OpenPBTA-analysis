suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("broom"))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(rlist))
suppressPackageStartupMessages(library(transcripTools))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))

option_list <- list(
  make_option(c("-r", "--cc_rdata"),
               type = "character", default = NULL,
               help = "matrix with features as rows, Sample as column"),
  make_option(c("-i", "--geneexp"),
              type = "character", default = NULL,
              help = "matrix with gene exp data"),
  make_option(c("-f", "--cohort_label_file"),
              type = "character", default = "../../../data/initial_hgat_histology.tsv",
              help = "main manifest for cohort "),
  make_option(c("-v", "--output_geneexp_tsv"),
              type = "character", default = "../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.tsv",
              help = "tsv format of the expressin file with desired cohort"),
  make_option(c("-p", "--patient_id_column"),
              type = "character", default = NULL,
              help = "patient id column id "),
  make_option(c("-n", "--nclust"),
              type = "integer", default = NULL,
              help = "Number of clusters in CC"),
  make_option(c("-c", "--colorfile"),
              type = "character", default = NULL,
              help = "File with colors for annotations"),
  make_option(c("-o", "--outprefix"),
              type = "character", default = NULL,
              help = "Output  prefix "),
  make_option(c("-t", "--print10kout"),
              type = "character", default = NULL,
              help = "YES or NO")
)

opt <- parse_args(OptionParser(option_list = option_list))
cohort_label_file <- opt$cohort_label_file
cc_rdata<-opt$cc_rdata
nclust <- opt$nclust
patient_id_column <- opt$patient_id_column
colorfile <- opt$colorfile
geneexp <- opt$geneexp
colorfile <- opt$colorfile
outprefix <- opt$outprefix
print10kout <- opt$print10kout
out1 <- paste(outprefix, "samplevssample.png", sep="_")
out2 <- paste(outprefix, "samplevsfeature.png", sep="_")
out_cluster_and_annotation <- paste(outprefix, "cluster_and_annotation.tsv", sep="_")
output_geneexp_tsv <- opt$output_geneexp_tsv

# read in histologies file 
cohort_plot_labels <- read_tsv(cohort_label_file) 
cohort_annotation_labels <- read_tsv(cohort_label_file) 

### The input RNAseq dataset
GeneExp <- readRDS(geneexp)

# read CC ckustering groups
load(cc_rdata)
CC_group<-result_CC[[nclust]]$consensusClass %>%
  as.data.frame()
colnames(CC_group)<-"CC"

# read in consensus clustering matrix
CC_consensus_mat <- result_CC[[nclust]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)


# reading color file
anno_color_list<-list.load(colorfile)
anno_colour = lapply(anno_color_list,function(x) unlist(x))


# add annotation to cluster groups
cluster<-as.data.frame(CC_group) %>%
  rownames_to_column() %>%
  # merge the annotation per cohort_participant_id
  dplyr::left_join(cohort_annotation_labels,by=c("rowname"= patient_id_column)) %>%
  dplyr::select("rowname", "CC", "integrated_diagnosis", "CNS_region") %>%
  unique() %>%
  remove_rownames() %>%
  column_to_rownames("rowname")

cluster$CC <- as.character(cluster$CC)

# the column name of the expression data is already Kids_First_Biospecimen_ID - no need to change anything - just write it out
write.table(GeneExp, output_geneexp_tsv, sep="\t")

#VST  transformation
GeneExp <- varianceStabilizingTransformation(round(as.matrix(GeneExp)),blind = TRUE,fitType="mean")

cluster_and_annotation_df <- tibble::rownames_to_column(cluster, "Sample.Names")
write.table(cluster_and_annotation_df, file.path("output",out_cluster_and_annotation), sep="\t", row.names=F)

# Creating heatmap with scaled gene expression and annotations
pheatmap::pheatmap(CC_consensus_mat,
    color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(9),
    annotation_col = cluster,
    annotation_colors = anno_colour,
    cluster_rows = result_CC[[nclust]]$consensusTree,
    cluster_cols=result_CC[[nclust]]$consensusTree,
    show_rownames = F,
    filename = file.path("plots",out1),
    fontsize = 10,cellwidth = 10, cellheight = 10)


# get most variable
#means <- rowMeans(GeneExp)
vars <- apply(GeneExp, 1, var)
varorder <- order(vars, decreasing = T)
top100_GeneExpzscored <- GeneExp[varorder, ] %>% head(100)
top10k_GeneExpzscored <- GeneExp[varorder, ] %>% head(10000)
#top100_GeneExpzscored <- mostVar(GeneExpzscored, 100, i_want_most_var = TRUE)


if (print10kout ==  "YES")
  GeneExp_top10k <- gsub(".rds", "_top10k.tsv", geneexp)
  write.table(top10k_GeneExpzscored, GeneExp_top10k, sep="\t")


# Creating z-score matrix from top100
top100_GeneExpzscored_means <-rowMeans(top100_GeneExpzscored, na.rm = TRUE)
top100_GeneExpzscored_sd <- apply(top100_GeneExpzscored, 1, sd, na.rm = TRUE)
# subtract mean
top100_GeneExpzscored_zscored <- sweep(top100_GeneExpzscored, 1, top100_GeneExpzscored_means, FUN = "-")
# divide by SD remove NAs and Inf values from zscore for genes with 0 in normData
top100_GeneExpzscored <- sweep(top100_GeneExpzscored_zscored, 1,top100_GeneExpzscored_sd, FUN = "/") %>% na_if(Inf) %>% na.omit()


# heatmap using the input data and consensusClass from CC
pheatmap::pheatmap(
    top100_GeneExpzscored, color = colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(13),
    breaks = c(10,8,6,4,2,1,0.1,-0.1,-1,-1.4,-1.8,-2, -3, -4),
    annotation_col = cluster,
    annotation_colors = anno_colour,
    cluster_cols=result_CC[[nclust]]$consensusTree,
    show_rownames = T,
    filename = file.path("plots",out2),
    fontsize = 10,cellwidth = 10, cellheight = 10)
