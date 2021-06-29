
# This script takes in a geneexp file for top 10K genes with cohort_participant_ids
# It also takes cluster number assigedn to each gene from clustering
# The output includes two files -
#       - genes with vtest higher  than 2 and lower than 2 in all other  clusters
#       - genes with vtest lower than -2 and highar than -2 in all other clusters

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as plt
import matplotlib.pyplot as plt
import pyreadr
import argparse
import sys
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--geneexpfile', required = True,
                    help = 'path to the geneexpfile file with cohort participant id as row names')
parser.add_argument('-c', '--clusterfile', required = True,
                    help = 'path to file with cohort_participantid and cluster number')
parser.add_argument('-t', '--clustertype', required = True,
                    help = 'column name with cluster type in cluster file')
parser.add_argument('-v', '--vtest_outname', required = True,
                    help = "Output name with vtest results")
args = parser.parse_args()

# Reading cluster file and renaming column name
clusteringtype = args.clustertype
clusters = pd.read_csv(args.clusterfile, sep="\t")
clusters = clusters.rename(columns={
    "Sample.Names": "cohort_participant_id"})[["cohort_participant_id", clusteringtype]]


mutsig_features = pd.read_csv(args.geneexpfile, sep="\t").T
mutsig_features = mutsig_features.reset_index().rename(columns={
    "index": "cohort_participant_id"}).set_index("cohort_participant_id")


## Merging cluster labels and geneexpfile scores
mutsig_with_cluster = clusters.merge(mutsig_features, on='cohort_participant_id')
mutsig_with_cluster = mutsig_with_cluster.iloc[:,1:] # Removing sample name


mutsig_with_cluster = mutsig_with_cluster.astype(float)
mutsig_with_cluster[clusteringtype] = mutsig_with_cluster[clusteringtype].astype(int)

## Mean and var(sigma sq) for all features from geneexpfile score datafram
mutsig_means = mutsig_with_cluster.iloc[:, 1:].mean(axis = 0)
mutsig_var = mutsig_with_cluster.iloc[:, 1:].var()


## This part of the code calculates the formula from slide 12 in these
###. slides - http://eric.univ-lyon2.fr/~ricco/cours/slides/en/classif_interpretation.pdf
mutsig_groupedbycluster = mutsig_with_cluster.groupby(clusteringtype)
vt_sign_values = pd.DataFrame({})


# iterate over each group
for group_name, df_group in mutsig_groupedbycluster:
    mutsig = df_group.iloc[:, 1:]
    num_samples = mutsig.shape[0]
    numerator = mutsig.mean() - mutsig_means
    denominator = np.sqrt(((len(mutsig_with_cluster) - num_samples)/(len(mutsig_with_cluster) -1)) * (mutsig_var/num_samples))
    vt_sign_values[group_name] = numerator/denominator



# Renaming and choosing improtant features based on v.test scores
vt_sign_values = vt_sign_values.reset_index().rename(columns={"index": "feature"})
vt_sign_values_plotinput_reformatted = vt_sign_values.set_index('feature').stack().reset_index(name='vt_score').rename(columns={'level_1':'cluster'})

vtest_highscore_out = args.vtest_outname.replace('.tsv', '.highvtscore.tsv')
vtest_lowscore_out = args.vtest_outname.replace(".tsv", ".lowvtscore.tsv")

# Using  vtest score, choosing only top and bottom in each cluster
    #  For loop around each gene feature,
    # Sort based on vt_score
    #  For top most features, last feature should have vtscore of >2 and
            # others need to be <2
    # For bottom features, first feature should have vtscore of <-2 and
            # other s need to be > -2
groupedby_features = vt_sign_values_plotinput_reformatted.groupby('feature')
up_regulated_df = pd.DataFrame(columns = ['feature', 'cluster', 'vt_score'])
down_regulated_df = pd.DataFrame(columns = ['feature', 'cluster', 'vt_score'])
for feat, df in tqdm(groupedby_features):
    df = df.sort_values('vt_score')
    print(feat, df)
    if df['vt_score'].iloc[-1] > 2 and df['vt_score'].iloc[-2] < 2:
        up_regulated_df =  up_regulated_df.append(df.iloc[-1:], ignore_index=True)
    elif df['vt_score'].iloc[0] < -2 and df['vt_score'].iloc[1] > -2:
        down_regulated_df = down_regulated_df.append(df.iloc[0], ignore_index=True)


# Printing up and down features into separate outfiles
up_regulated_df.to_csv(vtest_highscore_out, sep="\t", index=None)
down_regulated_df.to_csv(vtest_lowscore_out, sep="\t", index=None)


#
