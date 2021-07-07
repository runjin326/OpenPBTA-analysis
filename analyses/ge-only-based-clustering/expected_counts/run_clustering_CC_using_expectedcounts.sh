# spearman HC
# olfd run - Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID" -l hgat_all_primary.tsv -d "spearman" -a "hc"  -o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "hc" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_HC_UQ-pgQ2 \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# spearman PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "pam" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_PAM_UQ-pgQ2  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# euclidean HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "hc" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_HC_UQ-pgQ2  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# euclidian PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "pam" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_PAM_UQ-pgQ2  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# euclidean KM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "km" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_KM_UQ-pgQ2  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# manhattan PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "pam" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_PAM_UQ-pgQ2  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ

# manhattan HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "hc" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_HC_UQ-pgQ2 \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform UQ



# VST normalization
# spearman HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "hc" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_HC_expct_counts_VST \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# spearman PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "pam" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_VST  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# euclidean HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "hc" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_HC_expct_counts_VST  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# euclidian PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "pam" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_PAM_expct_counts_VST  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# euclidean KM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "km" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# manhattan HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "hc" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_HC_VST \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst

# manhattan PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "pam" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_PAM_expct_counts_VST  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform vst




## Using edgeR  normalization
# spearman HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "hc" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_HC_expct_counts_edgeR \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# spearman PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "spearman" -a "pam" \
-o output/CC/Distance_spearman_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# euclidean HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "hc" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_HC_expct_counts_edgeR  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# euclidian PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "pam" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# euclidean KM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "euclidean" -a "km" \
-o output/CC/Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_edgeR  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# manhattan HC
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "hc" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_HC_edgeR \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

# manhattan PAM
Rscript code/04-clustering-CC.R --cohort_id "initial_hgat" -p "Kids_First_Biospecimen_ID"   \
-l input/initial_hgat_histology.tsv -d "manhattan" -a "pam" \
-o output/CC/Distance_manhattan_finalLinkage_average_clusterAlg_PAM_expct_counts_edgeR  \
--gene_expression_file  input/initial_hgat-pbta-gene-counts-rsem-expected_count_corrected.proteincoding.rds \
--transform edgeR

#
