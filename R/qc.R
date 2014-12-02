qc <- function(trainingData, testData, geneList){
  if (sum(is.na(testData)>0)){stop("Found missing values in the test data")}
  if (!is.null(trainingData)){
    if (sum(is.na(trainingData)>0)){stop("Found missing values in the training data")}
    probe <- intersect(row.names(trainingData), row.names(testData))
    trainingData <- trainingData[match(probe, row.names(trainingData)),]
    testData <- testData[match(probe, row.names(testData)),]
  } else {
    probe <- row.names(testData)
  }
  for (i in 1:length(geneList)){
    geneList[[i]] <- geneList[[i]][geneList[[i]] %in% probe]
  }
  rtlist <- list(testData=testData, geneList=geneList)
  if (!is.null(trainingData)) {rtlist <- c(rtlist, list(trainingData=trainingData))}
  return(rtlist)
}