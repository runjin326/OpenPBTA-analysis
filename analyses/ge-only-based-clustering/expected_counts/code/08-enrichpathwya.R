# 1. This  script takes in a input with  the  following columns -
  # feature(gene  symbol), cluster number and vtscore
# 2. The script converts gene symbol to entrez ids
# 3. Using "enrichpathway", the  script generates overerrepresented features
# 4. Print out only feature  name and q-value in the output

suppressPackageStartupMessages(library(broom))
library(matrixStats)
library(tidyr)
library(readr)
library(tibble)
library(dplyr)
library(optparse)

library(ReactomePA)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

option_list <- list(
  make_option(c("-i", "--vtestfile"),
               type = "character", default = NULL,
               help = "File with features that show vtest score >/< 2 in each  cluster"),
  make_option(c("-c", "--clust"),
              type = "integer", default = NULL,
              help = "Number of clusters in vtest files"),
  make_option(c("-o", "--outprefix"),
              type = "character", default = NULL,
              help = "Prefix for output files")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Reading in input gene exp file with top 10k features
vtest = read.csv(opt$vtestfile, sep = '\t', header = TRUE)

# Making  a list with all cluster numbers
k <- seq(1, opt$clust, 1)

# For every cluster, convert gene symbols to entrezIDs  and
  # use enrichpathway to get overrepresented features
for (n in k) {
  vtest_highscores_cluster = subset(vtest, cluster==n)
  geneset = vtest_highscores_cluster[['feature']]
  geneentrezids <- as.numeric(mapIds(org.Hs.eg.db, geneset, 'ENTREZID', 'SYMBOL'))
  enrichOUT <- enrichPathway(gene=geneentrezids, organism = "human", pvalueCutoff = 0.05, readable = T)
  enrichOUT = as.data.frame(enrichOUT)[c("Description","qvalue")]
  top_features <- subset(enrichOUT, qvalue < 0.05)
  write.table(top_features, paste(opt$outprefix,"cluster",as.character(n),"overrep.tsv",sep = "_"),
              sep="\t", row.names=FALSE)
}
