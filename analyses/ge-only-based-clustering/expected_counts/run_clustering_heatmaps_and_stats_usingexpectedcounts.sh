## Subset expected counts RDS file
Rscript code/01-subset-cohort-samples.R \
-i ../../../data/pbta-gene-counts-rsem-expected_count-collapsed.combined.rds  \
-p ../../../data/pbta-histologies.tsv \
-c ../../../data/initial_hgat_histology.tsv \
-o  ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count.rds \
-t ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count.tsv

## Batch correction for samples of interest
Rscript code/02-batch-correct.R \
--metadata ../../../data/initial_hgat_histology.tsv \
--mat ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count.rds \
--output_prefix initial_hgat-pbta-gene-counts-rsem-expected_count 

## Filter genes in matrix to protein coding only
Rscript code/03-filter-protein-coding.R \
--input ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.rds \
--gencode_version 27 \
--output ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds

## The following scripts will be run after the parameters are selected based on: 
## run_clustering_CC_using_expectedcounts.sh
## Clustering using CC using UQ.pgQ2
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l initial_hgat_histology.tsv -d "spearman" -a "pam" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST  \
--gene_expression_file  initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# Heatmap using cluster number from CC
Rscript code/05-heatmaps.R \
-r output/CC/Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR/initial_hgat/CC.Rdata \
-p "Kids_First_Biospecimen_ID" -c input/color_withPCAcolor_3clusters.yaml \
-o CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR  \
-i ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds  \
-f ../../../data/initial_hgat_histology.tsv -t YES --nclust 5


# VTEST to pick immportant features
python3 code/06-vtest.py \
-i ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.tsv \
-c output/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_cluster_and_annotation.tsv  \
-o stats/vtest/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_VTEST \
-v stats/vtest/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_VTEST.tsv \
-t CC


## CHoose top features per cluster
python3 code/07-vtest-above-vscore2.py \
-i ../../../data/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding_top10k.tsv \
-c  output/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_cluster_and_annotation.tsv \
-t CC -v output/overrep_analysis/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_VTEST.tsv

## vtest for enrichpathway - For  each cluster output only the samples that >2 in one cluster and
#     <2 in all other clusters. Repeat with  bottom features except use <2 as cutoff
Rscript code/08-enrichpathway.R \
-i output/overrep_analysis/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_VTEST.highvtscore.tsv \
-c 5 -o output/overrep_analysis/highvtest_pathwayanalysis

## enrcich pathway
Rscript code/08-enrichpathway.R \
-i output/overrep_analysis/CC_based_heatmap_Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR_VTEST.lowvtscore.tsv \
-c 5 -o output/overrep_analysis/lowvtest_pathwayanalysis
