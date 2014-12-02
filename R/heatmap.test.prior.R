heatmap.test.prior <- function(diffGeneList, testData, trainingLabel,testLabel=NULL, coef_test, geneList=NULL){
  if (is.null(geneList)){
    nPath <- length(trainingLabel) - 1
    pathName <- names(trainingLabel)[-1]
  } else {
    nPath <- length(geneList)
    pathName <- names(geneList)
  }
  
  pdf("signature_heatmap_testset_prior.pdf")
  if (!is.null(testLabel)){
    cc <- as.numeric(as.factor(testLabel))
  }
  for (i in 1:nPath){
    tmp <- match(diffGeneList[[i]], row.names(testData))
    path <- testData[tmp, ]
    if (!is.null(testLabel)){
      heatmap(as.matrix(path[,order(coef_test[,i])]),Colv=NA, scale="row",ColSideColors=as.character(cc[order(coef_test[,i])]), col=bluered(128),margins = c(10,10), main=paste(pathName[i],"signature",sep=" "))
    } else {
      heatmap(as.matrix(path[,order(coef_test[,i])]),Colv=NA, scale="row", col=bluered(128),margins = c(10,10), main=paste(pathName[i],"signature",sep=" ")) 
    }
  }
  invisible(dev.off())
}
