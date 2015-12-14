# training data and label are available, but geneList is not available.

data_prep_s1 <- function(n_sigGene, trainingData, testData, trainingLabel, geneList, anchorGenes, excludeGenes, theta0, theta1)
{
  if (is.null(geneList)){  
    geneSelection <- bayes.gene.selection(n_sigGene, dat=trainingData, trainingLabel,iter=5000, burn_in=500, sigmaZero = 0.1, sigmaNonZero = 1, alpha_tau = 1, beta_tau = 0.01, p = 0.01)
    diffGeneList <- geneSelection$diffGeneList
  } else {
    diffGeneList <- geneList
  }
  
  #Add anchor genes to list if they aren't in there
  if(!is.null(anchorGenes)){
    for (j in 1:length(names(anchorGenes))){
      #if an anchor gene is not in the diffGeneList
      if(length(intersect(diffGeneList[[j]], anchorGenes[[j]])) != length(anchorGenes[[j]])){
        #for each anchor gene, if it's not there, replace last gene
        replaced <- 0
        for (k in 1:length(anchorGenes[[j]])){
          if(!(anchorGenes[[j]][k] %in% diffGeneList[[j]])){
            diffGeneList[[j]][length(diffGeneList[[1]])-replaced] <- anchorGenes[[j]][k]
            replaced <- replaced + 1
          }
        }
      }
    }
  }

  #Remove exclude genes if they are in there
  if(!is.null(excludeGenes)){
    for (j in 1:length(names(excludeGenes))){
      #if an exclude gene is in geneList
      if(length(intersect(diffGeneList[[j]], excludeGenes[[j]])) != 0){
        #remove any exclude genes
        diffGeneList[[j]] <- diffGeneList[[j]][!(diffGeneList[[j]] %in% excludeGenes[[j]])]
      }
    }
  }
  
  signature <- geneMatch_sub2(dat=trainingData, diffGeneList, trainingLabel)
  
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
  
  #change the Pi_matrix values for the anchor genes to 1
  if(!is.null(anchorGenes)){
    for (j in 1:length(names(anchorGenes))){
      for (k in 1:length(anchorGenes[[j]])){
        Pi_matrix[anchorGenes[[j]][k],names(anchorGenes)[j]] <- 1
      }
    }
  }
  return(list(trainingBaseline_sub = trainingBaseline_sub, S_matrix = signature$S_matrix, Delta_matrix = signature$Delta_matrix, Pi_matrix = Pi_matrix, trainingData_sub = trainingData_sub, testData_sub = testData_sub, diffGeneList = diffGeneList))
}
