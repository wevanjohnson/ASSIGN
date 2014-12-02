geneMatch_sub2 <- function(dat, diffGeneList, trainingLabel)
{  
  nPath <- length(trainingLabel) - 1
  bgPosB <- NULL; edPosB <- NULL
  for (i in 1:length(trainingLabel[[1]])){
    bgPosB <- c(bgPosB, trainingLabel[[1]][[i]][1])
    edPosB <- c(edPosB, trainingLabel[[1]][[i]][length(trainingLabel[[1]][[i]])])
  }
  bgPosS <- NULL; edPosS <- NULL
  for (i in 2:length(trainingLabel)){
    bgPosS <- c(bgPosS, trainingLabel[[i]][1])
    edPosS <- c(edPosS, trainingLabel[[i]][length(trainingLabel[[i]])])
  }	
  
  geneName2 <- diffGeneList
  geneName1 <- unique(unlist(geneName2))
  tmp <- match(geneName1, row.names(dat))
  dat1 <- dat[tmp, ]
  S1 <- matrix(0, nrow=nrow(dat1), ncol=length(bgPosB))
  S2 <- matrix(0, nrow=nrow(dat1), ncol=length(bgPosB))
  for (i in 1:length(bgPosB)){
    pathBase <- rowMeans(dat1[,bgPosB[i]:edPosB[i]])
    pathSig <- rowMeans(dat1[,bgPosS[i]:edPosS[i]]) - pathBase
    tmp2 <- match(geneName2[[i]], row.names(dat1))
    S1[, i] <- pathSig
    S2[tmp2,i] <- 1
  }
  rownames(S1) <- row.names(dat1)
  colnames(S1) <- names(trainingLabel)[-1]
  rownames(S2) <- row.names(dat1)
  colnames(S2) <- names(trainingLabel)[-1]
  return(list(S_matrix=S1, Delta_matrix=S2))
}