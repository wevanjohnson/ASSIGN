# Tests the ASSIGN preprocess function on the lung dataset
library(ASSIGN)
library(testthat)
data(testData1)
data(trainingData1)

context("Tests the ASSIGN preprocess function on the lung data")

test_that("test assign.preprocess on the lung data with multiple parameters", {
	set.seed(0)
	trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, ras=1:10, src=1:10),bcat = 11:19, e2f3 = 20:28, myc= 29:38, ras = 39:48, src = 49:55)

	testLabel1 <- rep(c("Adeno","Squamous"),c(53,58))

	processed.data <- assign.preprocess(trainingData=trainingData1,testData=testData1,trainingLabel=trainingLabel1,geneList=NULL, n_sigGene=rep(200,5))


	line750 = c(0.1438143, 0.0190799, 0.0703004, 0.3838708, 0.8276032)
	expect_equal(round(as.numeric(processed.data$"S_matrix"[750,]),7), line750)

})



