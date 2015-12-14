# only the geneList and test data is available
data_prep_s2 <- function(geneList, anchorGenes, testData, theta0, theta1){
  S <- matrix(0,nrow=nrow(testData),ncol=length(geneList))
  
  for (i in 1:length(geneList)){
    tmp <- match(geneList[[i]], row.names(testData))
    tmp <- tmp[!is.na(tmp)]
    S[tmp,i] <- 1
  }
  
  #Check for anchor Genes
  if(!is.null(anchorGenes)){
    for (j in 1:length(names(anchorGenes))){
      #if an anchor gene is not in the diffGeneList
      if(length(intersect(diffGeneList[[j]], anchorGenes[[j]])) != length(anchorGenes[[j]])){
        #FAIL with error message
        stop("All anchor genes must be listed in the geneList")
      }
    }
  }
  
  rownames(S) <- row.names(testData)
  colnames(S) <- names(geneList)
  
  S1 <- S[apply(S,1,sum)!=0,]
  S2 <- runif(length(S1),0,1); dim(S2) <- dim(S1)
  rownames(S2) <- rownames(S1); colnames(S2) <- colnames(S1)
  test <- testData[apply(S,1,sum)!=0,]
  
  B <- runif(NROW(S1),0,1)
  Pi_matrix <- as.matrix(ifelse(S1==0, theta0, theta1))
  
  #change the Pi_matrix values for the anchor genes to 1
  if(!is.null(anchorGenes)){
    for (j in 1:length(names(anchorGenes))){
      for (k in 1:length(anchorGenes[[j]])){
        Pi_matrix[anchorGenes[[j]][k],j] <- 1
      }
    }
  }
  
  return(list(testData_sub = test, Delta_matrix = S1, Pi_matrix = Pi_matrix, S_matrix = S2, B_vector = B))
}