assign.cv.output <- function(processed.data, mcmc.pos.mean.trainingData, trainingData, trainingLabel, adaptive_B=FALSE, adaptive_S=FALSE, mixture_beta=TRUE, outputDir){
  cat("Outputing results...\n")
  
  if (mixture_beta){
    coef_train = mcmc.pos.mean.trainingData$kappa_pos
  } else {
    coef_train = mcmc.pos.mean.trainingData$beta_pos
  }
  
  cwd <- getwd()
  setwd(outputDir)
  rownames(coef_train) <- names(processed.data$trainingData_sub)
  colnames(coef_train) <- names(trainingLabel)[-1]
  write.csv(coef_train, file="pathway_activity_trainingset.csv")
  
  #heatmaps of each pathway
  heatmap.train(diffGeneList=processed.data$diffGeneList, trainingData, trainingLabel)
  
  #provide test labels for model validation
  scatter.plot.train(coef_train, trainingData, trainingLabel)
  
  setwd(cwd)
}