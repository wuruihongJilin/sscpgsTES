# sscpgsTES
 single-sample classifier based on the paired genes for predicting tumor ecosystem subtypes of hepatocellular carcinoma
# Data input
Exp, a dataframe with mRNA expression profiles, samples in columns, genes in rows, rownames corresponding to Gene symbols. It requires log2 transformed RPKM/FPKM/TPM level data for RNA-seq dataset, and log2 level data ( after RMA or quantile normalization) for microarray data.
