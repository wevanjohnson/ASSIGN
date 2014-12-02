scatter.plot.test <- function(coef_test, trainingLabel, testLabel=NULL, geneList=NULL){
  if (is.null(geneList)){
    nPath <- length(trainingLabel) - 1
    pathName <- names(trainingLabel)[-1]
  } else {
    nPath <- length(geneList)
    pathName <- names(geneList)
  }
  
  if (!is.null(testLabel)){
    cc <- as.factor(testLabel)
  } else {
    cc <- rep(1, NROW(coef_test))
  }
  
  pdf("pathway_activity_scatterplot_testset.pdf")
  for(i in 1:nPath){
    ord <- order(coef_test[,i])
    cc_ord <- cc[ord]
    plot(coef_test[ord,i], col=cc_ord, pch=19, cex=0.7, xlab="Test sample", ylab=paste(pathName[i], "pathway activity", sep=" "), main=paste("Evaluate in vitro", pathName[i],  "signature in cancer samples",sep=" "))
    if (!is.null(testLabel)){
      legend("topleft", legend=unique(testLabel), col=unique(cc), pch=19, cex=0.7)
    }
  }
  invisible(dev.off())
}