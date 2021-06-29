### This script
##     1. reads in a tsv file of GeneExp
##     2. Reads in gencode GTF file
##     3. Filters for genes that only showed as "protein_coding" in GeneExp file and writes as output

import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import sys
import pyreadr

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputGenExp", required=True, help="eInput GeneExp file")
parser.add_argument(
    "-g", "--gtf_file", required=True, help="Gencode file in GTF format"
)
parser.add_argument("-o", "--outfilename", required=True, help="Out file name")
args = parser.parse_args()

out_text_name = args.outfilename.replace("rds", "tsv")

#Reading in gene Exp matrix
geneExp = pd.read_csv(args.inputGenExp, sep="\t")

# Reading in GTF  file and changing format to dataframe
gr = pr.read_gtf(args.gtf_file)
gencode_df = gr.df

# filtering only for protein coding
gencode_df_proteincoding = gencode_df[gencode_df["gene_type"] == "protein_coding"]
# Creating a gene dict with just protein coding from GTF dataframe
# Key is the gene name and value is just zzero
# Dict is mush faster that numpy array to search
gencode_gene_list = np.unique(gencode_df_proteincoding["gene_name"])
gencode_df_proteincoding_dict = dict(zip(gencode_gene_list,[0]*len(gencode_gene_list)))


# Filtering gene Exp matrix with genes from GTF dict
geneExp_filtered = geneExp.loc[[i in gencode_df_proteincoding_dict for i in geneExp.index]]
# Writing gene Exp matrix to output file
geneExp_filtered.to_csv(out_text_name, sep="\t", index=True)

# Also writing as RDS file for clustering input 
pyreadr.write_rds(args.outfilename, geneExp_filtered)



#
