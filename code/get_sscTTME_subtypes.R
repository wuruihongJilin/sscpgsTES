


#' predicting the transcriptomics sscTTME subtypes (1 to 7) of hepatocellular carcinoma based on the paired genes
#'
#' @param log2expdata the mRNA expression dataset, and the rowname is HGNC-symbol, the column is sampleID,
#' and the values are log2-transformed expression  at the levels such as UQFPKM/FPKM/TPM/fRMA
#' @param pairgenelist the paired genes between different subtypes, they are used to differ the subtypes
#' @param savedir the directory used to store the results
#' @param savename the filename for current dataset to be saved
#'
#' @return a list with two dataframe, sscTTME and statics
#' @export
#'
#' @examples
get_sscTTME_subtypes <- function(log2expdata,pairgenelist,savedir,savename){


  pcutoff <- 0.001# the p value cutoff used to determine whether the expression of the genes in first column of paired genes are significantly higher or lower than
  #the second column


  ###################
  res_all <- NULL#storing the final subtypes


  p1_mtx <- matrix(NA,ncol(log2expdata),length(pairgenelist))
  rownames(p1_mtx) <- colnames(log2expdata)
  colnames(p1_mtx) <- seq(1:ncol(p1_mtx))

  colnames(p1_mtx) <- names(pairgenelist)


  for (iii in 1:length(pairgenelist)){#for each of subtype pairs
    print(paste0("comparing: ",iii))

    res <- pairgenelist[[iii]]$genepair

    head(res)

    p1 <- matrix(NA,ncol(log2expdata),1)


    for (ttt in 1:ncol(log2expdata)){#for each sample

      v1 <- log2expdata[res$upgene_inGroupi_ivsj,ttt]
      v2 <- log2expdata[res$downgene_inGroupi_ivsj,ttt]
      if(sum(!is.na(v1-v2))>5){

        p1[ttt,1] <- wilcox.test(v1,v2,alternative ="greater",paired = TRUE)$p.value

      }


    }

    p1_mtx[,iii] <- p1

  }

  #head(p1_mtx[1:10,1:6])
  ########
  decisonclass2 <- matrix(NA,nrow(p1_mtx),ncol(p1_mtx))
  colnames(decisonclass2) <- colnames(p1_mtx)

  decisonclass2[p1_mtx <= pcutoff] <- 1# vote 1 if the P <= pcutoff
  decisonclass2[p1_mtx > pcutoff] <- 0#vote 0 if the P > pcutoff
  classvector <- as.numeric(gsub("_.*","",sub("class_","",colnames(p1_mtx))))
  decisonclass2 <- t(t(decisonclass2)*classvector)
  head(decisonclass2)

  uniclass <- sort(unique(classvector))#1 2 3 4 5 6 7
  class_votnumber <- matrix(NA,nrow(decisonclass2),length(uniclass))#storing the vote number for each class (column) for each tumor (row)

  colnames(class_votnumber) <- uniclass


  for (zzz in 1:length(uniclass)){
    class_votnumber[,zzz] <- apply(decisonclass2==uniclass[zzz],1,sum)

  }


  #############

  predict.cluster <- matrix(NA,ncol(log2expdata),1)
  multipleclass <- NA

  predict.cluster.new <- matrix(NA,ncol(log2expdata),1)

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




  p1_mtx <- format(p1_mtx,digits = 3)
  names(predict.cluster.new) <- colnames(log2expdata)

  colnames(class_votnumber) <- paste0("NO.vote.",colnames(class_votnumber))
  colnames(p1_mtx) <- paste0(colnames(p1_mtx),".Pwilcox_greater")
  res <- cbind(datasetID=savename,sampleID=rownames(p1_mtx),predict.cluster=predict.cluster.new[,1],predict.cluster_withmultiple=predict.cluster[,1],class_votnumber, p1_mtx)
  rownames(res) <- NULL


  ###### get the number and percentage of each subtypes in the input dataset

  res_all <- data.frame(res)#

  statics_all <- NULL


  tempnumber <- t(as.matrix(table(res_all$predict.cluster)))
  colnames(tempnumber) <- paste0("No.",colnames(tempnumber))

  tempsum <- apply(tempnumber,1,sum)

  ratio <- tempnumber/tempsum * 100
  colnames(ratio) <-sub("No.","Percent.",colnames(ratio))

  statics_all <- cbind(datasetID=savename,tempnumber,No.Total =tempsum,ratio)

  statics_all <- data.frame(statics_all)

  ### save the results

  write.csv(res,paste0(savedir,savename,"_predicted_sscTTME_subtypes.csv"))#
  write.csv(statics_all,paste0(savedir,savename,"_predicted_sscTTME_subtypes_statics.csv"))
  print(paste0("The results are in ",savedir))
  ## return the results
  classres <- list(sscTTME=res,statics=statics_all)


}
