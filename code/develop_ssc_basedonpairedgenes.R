#' developing the single sample classifier for hepatocellular carcinoma by making gene pairs between each pair of subtypes (1 tp 7)
#' then wilcox.test(paired) was used to test the difference between the first column and the second column of the paired genes
#' # finally, the subtype obtaining the highest number of votes will be assigned to the tested sample



rm(list=ls())
library(reshape2)

#########
dirfront <- "E:/01wrh/write/sscTTME/"
setwd(dirfront)

savedir <- paste0(dirfront,"results/")
dir.create(savedir)

savefront <- "SSC_MAD1500pless0.001degree1_orderbyPvalueTOP130_Fscore"# the name used to save the results

###read the results of the differentiall analysis
datadir <- paste0(dirfront,"data/")
degs <- read.csv(paste0(datadir,"CHCCHBV_159tumorNMF_7cluster_true_degswilcox.res.csv"))# the difference of each gene between each pair of subtypes
###
head(degs)

length(unique(degs$gene))

###########
degs$gene <- sub("[.]","-",degs$gene)
dim(degs)

#######
# the expression data of the discovery data set, CHCC-HBV cohort
expnew <- readRDS(paste0(datadir,"CHCC_HBV_proteincoding_UQFPKMtolog2p1_19539_159Tumor_maxMAD.RData"))
# the  true class
curclassfile <- data.frame(read.csv(file=paste0(datadir,"NMF_cluster7_groupConsensus.new_finalused.csv"),
                                    row.names = 1))

head(curclassfile)

table(curclassfile$x)
# 1  2  3  4  5  6  7
# 23 25 37 17 14 20 23

class <- curclassfile$x
names(class) <- rownames(curclassfile)

df <- reshape2::melt(data=cbind(gene=rownames(expnew),expnew))
dim(df)
head(df)

classvector <- class[match(df$variable,names(class))]

#############for  all pairs including the group i vs group j, and the group j vs group i
colnames(degs)
head(degs$gene.y)
head(degs$gene.x)
dim(degs)# 31500    30
head(degs[,13:20])

a1 <- cbind(gene=degs$gene,group1=degs$group1,group2=degs$group2,p=degs$p,
            p.adj=degs$p.adj,mean1 = degs$x.x.Mean,mean2=degs$x.y.Mean,
            diff1vs2_mean = degs$x.x.Mean-degs$x.y.Mean,sd1=degs$sd.x,sd2=degs$sd.y,
            median1=degs$x.x.Median,median2=degs$x.y.Median,diff1vs2_median = degs$x.x.Median-degs$x.y.Median)
head(a1)
a1[1:10,]
a2 <- cbind(gene=degs$gene,group1=degs$group2,group2=degs$group1,p=degs$p,
            p.adj=degs$p.adj,mean1 = degs$x.y.Mean,mean2=degs$x.x.Mean,
            diff1vs2 = degs$x.y.Mean-degs$x.x.Mean,sd1=degs$sd.y,sd2=degs$sd.x,
            median1=degs$x.y.Median,median2=degs$x.x.Median,diff1vs2_median = degs$x.y.Median-degs$x.x.Median)
head(a2)
a2[1:10,]
degs2 <- data.frame(rbind(a1,a2))
dim(degs2)

head(degs2)
degs2[,4:ncol(degs2)] <- apply(degs2[,4:ncol(degs2)],2,as.numeric)
degs2$meandividedsd_diff <- degs2$mean1/degs2$sd1-degs2$mean2/degs2$sd2
degs2[1:10,]
degs2$mediandividedsd_diff <- degs2$median1/degs2$sd1-degs2$median2/degs2$sd2


degs2[1:10,]

########### for each cluster get the degs number


unigroup <- unique(degs2$group1)# the unique class, 1 to 7 here
degs2 <- degs2[degs2$p < 0.001,]#limits the genes with p less than 0.001
dim(degs2)
max(degs2$p)

#############

##set the cutoff
degreecutoff  <- 1 # Each gene is allowed to form only one gene pair
topcutoff_degs <- 130#the top cutoff, if the number of degs(P < 0.001) between some subtype-pairs are more than this cutoff, only
#the top 'topcutoff_degs' will be used as candidates to form the gene pairs

initialabsolutecut <- 1 # the difference of expression level between a pair of genes at the log2 level
pcutoff <- 0.001## the p value cutoff used to determine whether the expression of the genes in first column of paired genes are significantly higher or lower than the second column
ratiodiff <- 0.75

#decisonclass <- matrix(NA,ncol(expnew),42)# 7 subtypes, therefore 42 combinations
plist <- vector(mode="list",length=7)
h=1
p1_mtx <- matrix(NA,ncol(expnew),42)
rownames(p1_mtx) <- colnames(expnew)
colnames(p1_mtx) <- seq(1:ncol(p1_mtx))


pairnumber <- matrix(NA,42,3)
colnames(pairnumber) <- c("No.genepair","No.upgeneinGroup1","No.downgeneinGroup1")
rownames(pairnumber) <- 1:nrow(pairnumber)
pairgenelist <- vector(mode="list",length=42)


for (i in 1:length(unigroup)){

  for (j in 1:length(unigroup)){
    print(paste0(i,"_",j))

    if(i!=j){

      curdegsUP <- degs2[degs2$group1==unigroup[i] & degs2$group2==unigroup[j] & degs2$diff1vs2_median > 0,]
      dim(curdegsUP)

      curdegsDown <- degs2[degs2$group1==unigroup[i] & degs2$group2==unigroup[j] & degs2$diff1vs2_median < 0,]
      dim(curdegsDown)


      if(nrow(curdegsUP) > topcutoff_degs){
        curdegsUP <- curdegsUP[order(-log10(curdegsUP$p),decreasing = TRUE),]
        curdegsUP <- curdegsUP[1:topcutoff_degs,]


      }else{

        curdegsUP <- curdegsUP[order(-log10(curdegsUP$p),decreasing = TRUE),]

      }
      if(nrow(curdegsDown) > topcutoff_degs){
        curdegsDown <- curdegsDown[order(-log10(curdegsDown$p),decreasing = TRUE),]
        curdegsDown <- curdegsDown[1:topcutoff_degs,]

      }else{
        curdegsDown <- curdegsDown[order(-log10(curdegsDown$p),decreasing = TRUE),]


      }
      curdegsDown$gene

      dim(curdegsDown)
      dim(curdegsUP)
      curdegsUP$gene
      curdegsDown$gene



      ###

      zup1 <- matrix(curdegsUP$median1,nrow(curdegsUP),1)

      rownames(zup1) <- curdegsUP$gene

      zdown1 <- matrix(curdegsDown$median1,nrow(curdegsDown),1)

      rownames(zdown1) <- curdegsDown$gene

      zupMtx1 <- matrix(NA,nrow(zup1),nrow(zdown1))
      zupMtx1[,1:ncol(zupMtx1)] <- zup1

      zdownMtx1 <- matrix(NA,nrow(zdown1),nrow(zup1))
      zdownMtx1[,] <- zdown1
      zdownMtx1 <- t(zdownMtx1)


      diffmtx12org <- (zupMtx1-zdownMtx1)
      rownames(diffmtx12org) <- curdegsUP$gene
      colnames(diffmtx12org) <- curdegsDown$gene
      df12 <- reshape2::melt(diffmtx12org)
      df12$Var1 <- as.vector(df12$Var1)
      df12$Var2 <- as.vector(df12$Var2)


      diffmtx12 <- (zupMtx1-zdownMtx1) > initialabsolutecut###
      diffmtx12[isTRUE(diffmtx12)] <- 1
      diffmtx12[isFALSE(diffmtx12)] <- 0


      rownames(diffmtx12) <- curdegsUP$gene
      colnames(diffmtx12) <- curdegsDown$gene

      #head(diffmtx12[,1:10])

      colsum <- apply(diffmtx12,2,sum)

      table(colsum)

      diffmtx12

      ###diff21 ###############
      zdown2 <- matrix(curdegsUP$median2,nrow(curdegsUP),1)

      rownames(zdown2) <- curdegsUP$gene

      zup2 <- matrix(curdegsDown$median2,nrow(curdegsDown),1)

      rownames(zup2) <- curdegsDown$gene

      zupMtx2 <- matrix(NA,nrow(zup2),nrow(zdown2))
      zupMtx2[,] <- (as.vector(zup2[,1]))################## this is according to the columns always
      zupMtx2 <- t(zupMtx2)



      zdownMtx2 <- matrix(NA,nrow(zdown2),nrow(zup2))
      zdownMtx2[,] <- zdown2


      diffmtx21org <- (zupMtx2-zdownMtx2)
      rownames(diffmtx21org) <- curdegsUP$gene
      colnames(diffmtx21org) <- curdegsDown$gene
      df21 <- reshape2::melt(diffmtx21org)
      df21$Var1 <- as.vector(df21$Var1)
      df21$Var2 <- as.vector(df21$Var2)

      df <- data.frame(cbind(upgene_inGroupi_ivsj=df12$Var1,
                             downgene_inGroupi_ivsj=df12$Var2,
                             AminusB_inGroupi=data.frame(df12$value),
                             AminusB_inGroupj=data.frame(df21$value)),
                       stringsAsFactors = FALSE)

      dim(df)
      head(df)

      loc <- which(df$df12.value > initialabsolutecut & df$df21.value > initialabsolutecut)
      length(loc)#
      temptes <- df[loc,]

      dim(temptes)


      #for each pair, we get the ratio with A>B in each group

      curexp1 <- expnew[,rownames(curclassfile)[curclassfile$x==i]]
      curexp2 <- expnew[,rownames(curclassfile)[curclassfile$x==j]]
      dim(curexp1)
      dim(curexp2)

      sub12_diff <- matrix(NA,nrow(temptes),ncol(curexp1))
      sub21_diff <- matrix(NA,nrow(temptes),ncol(curexp2))


      for(yyy in 1:ncol(curexp1)){
        sub1 <- curexp1[temptes$upgene_inGroupi_ivsj,yyy]
        sub2 <- curexp1[temptes$downgene_inGroupi_ivsj,yyy]
        sub12_diff[,yyy] <- sub1-sub2
      }
      dim(sub12_diff)
      ratio1 <- apply(sub12_diff>0,1,sum)/ncol(sub12_diff)
      ratio1

      for(yyy in 1:ncol(curexp2)){
        sub1 <- curexp2[temptes$upgene_inGroupi_ivsj,yyy]
        sub2 <- curexp2[temptes$downgene_inGroupi_ivsj,yyy]
        sub21_diff[,yyy] <- sub2-sub1
      }
      dim(sub21_diff)
      ratio2 <- apply(sub21_diff>0,1,sum)/ncol(sub21_diff)
      ratio2

      ratio <- cbind(ratio1,ratio2)
      ratio <- cbind(ratio,Fscore=(ratio[,1]+ratio[,2])/2)
      temptes.2 <- cbind(temptes,ratio)

      temptes.2 <- temptes.2[order(temptes.2$Fscore,decreasing = TRUE),]

      temptes.2 <- temptes.2[temptes.2$ratio1>ratiodiff & temptes.2$ratio2>ratiodiff,]

      dim(temptes.2)
      keptpair <- temptes.2[1,]



      for(sss in 2:nrow(temptes.2)){
        v1 <- sum(keptpair$upgene_inGroupi_ivsj %in% temptes.2[sss,]$upgene_inGroupi_ivsj)
        v2 <- sum(keptpair$downgene_inGroupi_ivsj %in% temptes.2[sss,]$downgene_inGroupi_ivsj)

        if(v1 < degreecutoff & v2 < degreecutoff){

          keptpair <- rbind(keptpair,temptes.2[sss,])

        }



      }



      ####
      keptpair <- unique(keptpair)
      dim(keptpair)#
      keptpair


      length(unique(keptpair$upgene_inGroupi_ivsj))

      length(unique(keptpair$downgene_inGroupi_ivsj))


      #########

      res <- keptpair
      res$upgene_inGroupi_ivsj <- as.vector(res$upgene_inGroupi_ivsj)
      res$downgene_inGroupi_ivsj <- as.vector(res$downgene_inGroupi_ivsj)

      head(res)

      p1 <- matrix(NA,ncol(expnew),1)


      for (ttt in 1:ncol(expnew)){
        v1 <- expnew[res$upgene_inGroupi_ivsj,ttt]
        v2 <- expnew[res$downgene_inGroupi_ivsj,ttt]
        p1[ttt,1] <- wilcox.test(v1,v2,alternative ="greater",paired = TRUE)$p.value


      }


      namethis <- paste0("class_",i,"_vs_",j)

      pairgenelist[[h]] <- list(genepair=res,
                                uniUpGene=sort(unique(res$upgene_inGroupi_ivsj)),
                                uniDownGene=sort(unique(res$downgene_inGroupi_ivsj)))

      names(pairgenelist)[h] <- namethis



      pairnumber[h,1:3] <- c(nrow(res),length(unique(res$upgene_inGroupi_ivsj)),length(unique(res$downgene_inGroupi_ivsj)))
      rownames(pairnumber)[h] <- namethis

      p1_mtx[,h] <- p1
      colnames(p1_mtx)[h] <- namethis


      h=h+1

    }




  }

}

names(pairgenelist)



decisonclass2 <- matrix(NA,nrow(p1_mtx),ncol(p1_mtx))
#colnames(decisonclass2) <- colnames(decisonclass)

decisonclass2[p1_mtx <= pcutoff] <- 1
decisonclass2[p1_mtx > pcutoff] <- 0
classvector <- as.numeric(gsub("_.*","",sub("class_","",colnames(p1_mtx))))
decisonclass2 <- t(t(decisonclass2)*classvector)
head(decisonclass2)

uniclass <- sort(unique(classvector))
class_votnumber <- matrix(NA,nrow(decisonclass2),length(uniclass))
colnames(class_votnumber) <- uniclass

for (i in 1:length(uniclass)){
  class_votnumber[,i] <- apply(decisonclass2==uniclass[i],1,sum)

}


#############
#votenumber <- apply(decisonclass,1,table)

predict.cluster <- matrix(NA,ncol(expnew),1)
multipleclass <- NA

predict.cluster.new <- matrix(NA,ncol(expnew),1)

for(kkk in 1:nrow(class_votnumber)){
  curnames <- colnames(class_votnumber)
  curvote <- class_votnumber[kkk,which.max(class_votnumber[kkk,])]

  curclassname <- curnames[class_votnumber[kkk,]==curvote]

  predict.cluster[kkk] <- paste(curclassname,collapse = ",")

  multipleclass[kkk] <-  length(curclassname)

  ############ if multiple class, get the pairwise p value among the clusters, the get the one with the minimum p value


  if(length(curclassname) > 1){
    name <- c()
    for(ppp in 1:length(curclassname)){

      for(hhh in setdiff(1:length(curclassname),ppp)){

        name <- c(name,paste0("class_",curclassname[ppp],"_vs_",curclassname[hhh]))
      }
    }

    tempvalue <- p1_mtx[kkk,colnames(p1_mtx) %in% name]


    predict.cluster.new[kkk] <- stringr::str_split(names(which.min(tempvalue)),"_",simplify = TRUE)[2]




  }else{

    predict.cluster.new[kkk] <- predict.cluster[kkk]
  }

}


multipleclass
tempres <- cbind(multi=predict.cluster[multipleclass>1],decisonclass2[multipleclass>1,])

table(predict.cluster)
#predict.cluster
# 1   1,2 1,2,3   1,3   1,5     2   2,3   2,7     3   3,4   3,5     4     5   5,7     6   6,7     7
# 13    11     1     1     1    16     1     2    37     2     1    15    13     1    17     7    20
table(predict.cluster.new)
# predict.cluster.new
# 1  2  3  4  5  6  7
# 22 24 37 17 14 22 23
table(curclassfile$x)
# 1  2  3  4  5  6  7
# 23 25 37 17 14 20 23

#table(predict.cluster,curclassfile$x)

table(predict.cluster.new,curclassfile$x)# the crosstable


sum(diag(table(predict.cluster.new,curclassfile$x)))/length(predict.cluster.new)#0.9559748, the accuracy


names(pairgenelist)


pairnumber
####################

saveRDS(pairgenelist,paste0(savedir,savefront,"_pairgenelist_degreeless",degreecutoff,".rds"))
write.csv(pairnumber,paste0(savedir,savefront,"_pairgenelist_statics_degreeless",degreecutoff,".csv"))
write.csv(p1_mtx,paste0(savedir,savefront,"_ppairgenelist_degreeless",degreecutoff,"wixcox_p1_mtx.csv"))
write.csv(cbind(sampleID=rownames(curclassfile),
                Trueclass=curclassfile$x,
                predict.cluster=predict.cluster[,1],
                No.class=multipleclass,
                predict.cluster.modify=predict.cluster.new[,1]),
          paste0(savedir,savefront,"_pairgenelist_clusterpredict.csv"))


## showing the number of gene pairs
barplot(pairnumber[,1])

