# training data and label are available, but geneList is not available.

data_prep_s1 <- function(n_sigGene, trainingData, testData, trainingLabel, geneList, theta0, theta1)
{
  if (is.null(geneList)){  
    geneSelection <- bayes.gene.selection(n_sigGene, dat=trainingData, trainingLabel,iter=500, burn_in=100, sigmaZero = 0.1, sigmaNonZero = 1, alpha_tau = 1, beta_tau = 0.01, p = 0.01)
    diffGeneList <- geneSelection$diffGeneList
    signature <- geneMatch_sub2(dat=trainingData, diffGeneList, trainingLabel)
  } else {
    diffGeneList <- geneList
    signature <- geneMatch_sub2(dat=trainingData, diffGeneList, trainingLabel)
  }
  
  tmp1 <- match(unique(unlist(diffGeneList)), row.names(trainingData))
  
  trainingBaseline_sub <- rowMeans(trainingData[tmp1, unique(unlist(trainingLabel$control))])
  trainingData_sub <- trainingData[tmp1, ]
  testData_sub <- testData[tmp1, ]
  
  if (is.null(geneList)){
    Pi_matrix <- geneSelection$prior_p[tmp1, ]
    Pi_matrix <- Pi_matrix * signature$Delta_matrix #??Do we force the non-signatue genes to be 0?? 
    Pi_matrix <- theta0 + Pi_matrix * theta1
  } else {
    Pi_matrix <- signature$Delta_matrix
    Pi_matrix <- theta0 + Pi_matrix * theta1
  }
  return(list(trainingBaseline_sub = trainingBaseline_sub, S_matrix = signature$S_matrix, Delta_matrix = signature$Delta_matrix, Pi_matrix = Pi_matrix, trainingData_sub = trainingData_sub, testData_sub = testData_sub, diffGeneList = diffGeneList))
}
