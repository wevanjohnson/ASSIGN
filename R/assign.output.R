assign.output <- function(processed.data, mcmc.pos.mean.testData, trainingData, testData, trainingLabel, testLabel, geneList, adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE, outputDir){
  cat("Outputing results...\n")
  
  if (mixture_beta){
    coef_test = mcmc.pos.mean.testData$kappa_pos
  } else {
    coef_test = mcmc.pos.mean.testData$beta_pos
  }
  
  cwd <- getwd()
  setwd(outputDir)
  
  if (is.null(geneList)){
    pathName <- names(trainingLabel)[-1]
  } else {
    pathName <- names(geneList)
  }
  
  rownames(coef_test) <- colnames(processed.data$testData_sub)
  colnames(coef_test) <- pathName
  write.csv(coef_test, file="pathway_activity_testset.csv")
  
  #heatmaps of each pathway
  if (!is.null(trainingData) & !is.null(trainingLabel)){
    heatmap.train(diffGeneList=processed.data$diffGeneList, trainingData, trainingLabel)
  }
  heatmap.test.prior(diffGeneList=processed.data$diffGeneList, testData, trainingLabel,testLabel, coef_test, geneList)
  if (adaptive_S){
    heatmap.test.pos(testData=processed.data$testData_sub, Delta_pos = mcmc.pos.mean.testData$Delta_pos, trainingLabel, testLabel, Delta_cutoff = 0.95,coef_test, geneList)
  }  
  
  #provide test labels for model validation
  scatter.plot.test(coef_test, trainingLabel, testLabel, geneList)
  if (!is.null(testLabel)){
    box.plot.test(coef_test, trainingLabel, testLabel, geneList)
  }
  
  setwd(cwd)
}



