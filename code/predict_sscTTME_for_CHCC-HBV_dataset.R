
## for CHCC HBV dataset


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



log2expdata <- readRDS(paste0(datadir,"CHCC_HBV_proteincoding_UQFPKMtolog2p1_19539_159Tumor_maxMAD.RData"))

pairgenelist <- readRDS(paste0(datadir,"SSC_MAD1500pless0.001degree1_orderbyPvalueTOP130_Fscore_pairgenelist_degreeless1.rds"))

results <- get_sscTTME_subtypes(log2expdata,pairgenelist,savedir,savename="CHCC_HBV")
  
