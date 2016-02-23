#' Prediction/validation output for test data
#' 
#' The assign.output function outputs the summary results and plots for
#' prediction/validation for the test dataset.
#' 
#' The assign.output function is suggested to run after the assign.preprocess,
#' assign.mcmc and assign.summary functions. For the prediction/validation in
#' the test dataset, The Y argument in the assign.mcmc function is the output
#' value "testData_sub" from the assign.preprocess function.
#' 
#' @param processed.data The list object returned from the assign.preprocess
#' function.
#' @param mcmc.pos.mean.testData The list object returned from the assign.mcmc
#' function. Notice that for prediction/validation in the test dataset, the Y
#' argument in the assign.mcmc function should be set as the test dataset.
#' @param trainingData The genomic measure matrix of training samples (i.g.,
#' gene expression matrix). The dimension of this matrix is probe number x
#' sample number.
#' @param testData The genomic measure matrix of test samples (i.g., gene
#' expression matrix). The dimension of this matrix is probe number x sample
#' number.
#' @param trainingLabel The list linking the index of each training sample to a
#' specific group it belongs to.
#' @param testLabel The vector of the phenotypes/labels of the test samples.
#' @param geneList The list that collects the signature genes of one/multiple
#' pathways. Every component of this list contains the signature genes
#' associated with one pathway.
#' @param adaptive_B Logicals. If TRUE, the model adapts the
#' baseline/background (B) of genomic measures for the test samples. The
#' default is TRUE.
#' @param adaptive_S Logicals. If TRUE, the model adapts the signatures (S) of
#' genomic measures for the test samples. The default is FALSE.
#' @param mixture_beta Logicals. If TRUE, elements of the pathway activation
#' matrix are modeled by a spike-and-slab mixuture distribution. The default is
#' TRUE.
#' @param outputDir The path to the directory to save the output files. The
#' path needs to be quoted in double quotation marks.
#' @return The assign.output returns one .csv file containing one/multiple
#' pathway activity for each individual test samples, scatter plots of pathway
#' activity for each individual pathway in all the test samples, and heatmap
#' plots for the gene expression of the prior signature and posterior signtures
#' (if adaptive_S equals TRUE) of each individual pathway in the test samples.
#' @author Ying Shen
#' @examples
#' 
#' \dontshow{
#' setwd(tempdir())
#' tempdir <- tempdir()
#' 
#' data(trainingData1)
#' data(testData1)
#' data(geneList1)
#' 
#' trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, ras=1:10, 
#' src=1:10), bcat = 11:19, e2f3 = 20:28, myc= 29:38, ras = 39:48, src = 49:55)
#' testLabel1 <- rep(c("subtypeA","subtypeB"),c(53,58))
#' 
#' processed.data <- assign.preprocess(trainingData=trainingData1, 
#' testData=testData1, trainingLabel=trainingLabel1, geneList=geneList1)
#' mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, 
#' Bg = processed.data$B_vector, X=processed.data$S_matrix, 
#' Delta_prior_p = processed.data$Pi_matrix, iter = 20, 
#' adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE)
#' mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=10, iter=20, 
#' adaptive_B=FALSE, adaptive_S=FALSE,mixture_beta=TRUE)
#' }
#' assign.output(processed.data=processed.data, 
#' mcmc.pos.mean.testData=mcmc.pos.mean, trainingData=trainingData1, 
#' testData=testData1, trainingLabel=trainingLabel1, testLabel=testLabel1, 
#' geneList=NULL, adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE, outputDir=tempdir)
#' 
#' @export assign.output
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



