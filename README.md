# sscpgsTES
 single-sample classifier based on the paired genes for predicting tumor ecosystem subtypes of hepatocellular carcinoma
# Data input
Exp, a dataframe with mRNA expression profiles, samples in columns, genes in rows, rownames corresponding to Gene symbols. It requires log2 transformed RPKM/FPKM/TPM level data for RNA-seq dataset, and log2 level data ( after RMA or quantile normalization) for microarray data.

# an example
## for CHCC HBV dataset


rm(list=ls())
##source the program get_sscTTME_subtypes.R
source("code/get_sscTTME_subtypes.R")

##create a directory for saving the result
savedir <- "E:/sscpgsTES/"
dir.create(savedir)

##load the mRNA expression data frame

log2expdata <- readRDS(paste0(datadir,"CHCC_HBV_proteincoding_UQFPKMtolog2p1_19539_159Tumor_maxMAD.RData"))
cursavename <- "CHCC_HBV"#the dataset name

## load the paired genes of this classifier
pairgenelist <- readRDS(paste0(datadir,"SSC_MAD1500pless0.001degree1_orderbyPvalueTOP130_Fscore_pairgenelist_degreeless1.rds"))

## predict the subtypes, the result file will be stored in the savedir
results <- get_sscTTME_subtypes(log2expdata,pairgenelist,savedir,savename=cursavename)


