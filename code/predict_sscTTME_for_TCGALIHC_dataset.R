
## for TCGA_LIHC dataset

rm(list=ls())

dirfront <- "E:/01wrh/write/sscTTME/"
setwd(dirfront)

source("code/get_sscTTME_subtypes.R")
savedir <- paste0(dirfront,"results/")

##
datadir <- paste0(dirfront,"data/")
if(!dir.exists(savedir)){
  dir.create(savedir)
  
}

log2expdata <- readRDS(paste0(datadir,"TCGA_LIHC_log2p1FPKM_proteincoding_HGNCsymbolslog2exp_atgenelevel_log2p1FPKM_19626_370Tumor_maxMAD.Rdata"))


pairgenelist <- readRDS(paste0(datadir,"SSC_MAD1500pless0.001degree1_orderbyPvalueTOP130_Fscore_pairgenelist_degreeless1.rds"))

results <- get_sscTTME_subtypes(log2expdata,pairgenelist,savedir,savename="TCGA_LIHC")
  
