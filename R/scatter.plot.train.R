scatter.plot.train <- function(coef_train, trainingData, trainingLabel){
  nPath <- length(trainingLabel) - 1
  trainL <- character(length = length(unlist(trainingLabel)))
  for (i in 1:(nPath+1)){
    if (i == 1){
      x <- unique(trainingLabel[[i]])
      names(x) <- paste("control", 1:length(x), sep="")
      for (j in 1:length(x)){
        trainL[x[[j]]] <- rep(names(x)[j], length(x[[j]]))
      }
    } else {
      trainL[trainingLabel[[i]]] <- rep(names(trainingLabel)[i], length(trainingLabel[[i]]))
    }
  }
  trainL <- trainL[trainL != ""]
  
  pdf("pathway_activity_scatterplot_trainingset.pdf")
  for (i in 1:nPath){
    HMEC_samples <- 1:ncol(trainingData)
    Pathway_strength_HMEC <- coef_train[,i]
    plot(HMEC_samples, Pathway_strength_HMEC, col=as.factor(trainL),xlab="HMEC sample", ylab=paste(names(trainingLabel)[i+1], "pathway activity", sep=" "), main=paste("Cross-validation in HMEC", names(trainingLabel)[i+1], "pathway", sep=" "), pch=19, cex=0.7)
    legend("topleft", legend=unique(trainL), pch=19, cex=0.7, col=as.numeric(as.factor(unique(trainL))))
  }  
  invisible(dev.off())
}