assign.preprocess <- function(trainingData=NULL, testData, anchorGenes=NULL, excludeGenes=NULL,
                              trainingLabel, geneList=NULL, n_sigGene=NA, theta0=0.05, theta1=0.9){
  cat("Runing ASSIGN development version: truncated_S\n")
  cat("Performing QC on the input data...\n")
  dat <- qc(trainingData, testData, geneList)
  if (!is.null(trainingLabel)){
    if (identical(names(trainingLabel$control), names(trainingLabel)[-1]) == FALSE){warning("Control Labels DO NOT match the experimental Labels!\nPlease make sure that you specify the correct indice for control and experimental samples in the trainingLabel!")}
  }
  
  cat("Generating starting/prior values for model parameters...\n")
  if (is.null(trainingData) & is.null(geneList)){
    stop("trainingData and geneList are both set NULL. Need one of them for the analysis!")
  }
  if (!is.null(trainingData) & !is.null(trainingLabel)){
    x <- data_prep_s1(n_sigGene, trainingData=dat$trainingData, testData=dat$testData,
                      trainingLabel, anchorGenes, excludeGenes, geneList=dat$geneList,
                      theta0, theta1)
    return(list(trainingData_sub=x$trainingData_sub, testData_sub=x$testData_sub,
                B_vector=x$trainingBaseline_sub, S_matrix=x$S_matrix,
                Delta_matrix=x$Delta_matrix, Pi_matrix=x$Pi_matrix, diffGeneList=x$diffGeneList))
  }
  else if (is.null(trainingData) & is.null(trainingLabel) & !is.null(geneList)) {
    x <- data_prep_s2(geneList=dat$geneList, anchorGenes, testData=dat$testData, theta0, theta1)
    return(list(testData_sub=x$testData_sub, B_vector=x$B_vector, S_matrix=x$S_matrix,
                Delta_matrix=x$Delta_matrix, Pi_matrix=x$Pi_matrix, diffGeneList=dat$geneList))
  }
}
