# only the geneList and test data is available
data_prep_s2 <- function(geneList, testData, theta0, theta1){
  S <- matrix(0,nrow=nrow(testData),ncol=length(geneList))
  
  for (i in 1:length(geneList)){
    tmp <- match(geneList[[i]], row.names(testData))
    tmp <- tmp[!is.na(tmp)]
    S[tmp,i] <- 1
  }
  
  rownames(S) <- row.names(testData)
  colnames(S) <- names(geneList)
  
  S1 <- S[apply(S,1,sum)!=0,]
  S2 <- runif(length(S1),0,1); dim(S2) <- dim(S1)
  rownames(S2) <- rownames(S1); colnames(S2) <- colnames(S1)
  test <- testData[apply(S,1,sum)!=0,]
  
  B <- runif(NROW(S1),0,1)
  Pi_matrix <- ifelse(S1==0, theta0, theta1)
  
  return(list(testData_sub = test, Delta_matrix = S1, Pi_matrix = Pi_matrix, S_matrix = S2, B_vector = B))
}