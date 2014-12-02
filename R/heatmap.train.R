heatmap.train <- function(diffGeneList, trainingData, trainingLabel){
  nPath <- length(trainingLabel) - 1
  
  bgPosB <- NULL; edPosB <- NULL
  for (i in 1:length(trainingLabel[[1]])){
    bgPosB <- c(bgPosB, trainingLabel[[1]][[i]][1])
    edPosB <- c(edPosB, trainingLabel[[1]][[i]][length(trainingLabel[[1]][[i]])])
  }
  bgPosS <- NULL; edPosS <- NULL
  for (i in 2:length(trainingLabel)){
    bgPosS <- c(bgPosS, trainingLabel[[i]][1])
    edPosS <- c(edPosS, trainingLabel[[i]][length(trainingLabel[[i]])])
  }
  
  pdf("signature_heatmap_trainingset.pdf")
  for (i in 1:nPath){
    tmp <- match(diffGeneList[[i]], row.names(trainingData))
    path <- trainingData[tmp, ]
    sig <- rowMeans(path[ ,bgPosS[i]:edPosS[i]]) - rowMeans(path[ ,bgPosB[i]:edPosB[i]])
    ord <- order(sig)
    path_ord <- path[ord, c(bgPosB[i]:edPosB[i], bgPosS[i]:edPosS[i])]
    cc <- rep(c(1,2),c((edPosB[i] - bgPosB[i] + 1), (edPosS[i] - bgPosS[i] + 1)))
    heatmap(as.matrix(path_ord),scale="row",Rowv=NA,Colv=NA,ColSideColors=as.character(cc),col=bluered(128),margins = c(10,10), main=paste(names(trainingLabel)[i+1],"signature",sep=" "))
  }
  invisible(dev.off())
}