# sscpgsTES
 single-sample classifier based on the paired genes for predicting tumor ecosystem subtypes of hepatocellular carcinoma
# Data input
log2expdata, a dataframe with mRNA expression profiles, samples in columns, genes in rows, rownames corresponding to Gene symbols. It requires log2 transformed RPKM/FPKM/TPM level data for RNA-seq dataset, and log2 level data ( after RMA or quantile normalization) for microarray data.

# An example using CHCC-HBV dataset
rm(list=ls())

source("code/get_sscTTME_subtypes.R")#source the program get_sscTTME_subtypes.R

savedir <- "E:/sscpgsTES/"

dir.create(savedir)#create a directory for saving the result

log2expdata <- readRDS("CHCC_HBV_proteincoding_UQFPKMtolog2p1_19539_159Tumor_maxMAD.RData")#load the mRNA expression data frame
cursavename <- "CHCC_HBV"#the dataset name

pairgenelist <- readRDS("SSC_MAD1500pless0.001degree1_orderbyPvalueTOP130_Fscore_pairgenelist_degreeless1.rds")#load the paired genes of this classifier

results <- get_sscTTME_subtypes(log2expdata,pairgenelist,savedir,savename=cursavename)#predict the subtypes, the result file will be stored in the savedir

# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this classifer and datasets to: 
the Corresponding authors: Yunpeng Zhang; Shangwei Ning; Xia Li. College of Bioinformatics Science and Technology, Harbin Medical University, Harbin, Heilongjiang; E-mail: zhangyp@hrbmu.edu.cn (Yunpeng Zhang); ningsw@ems.hrbmu.edu.cn (Shangwei Ning); lixia@hrbmu.edu.cn (Xia Li)
or the first author, Ruihong Wu, by wuruihong@jlu.edu.cn
