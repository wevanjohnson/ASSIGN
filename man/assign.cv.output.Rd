\name{assign.cv.output}
\alias{assign.cv.output}
\title{
Cross validation output
}
\description{
The assign.cv.output function outputs the summary results and plots for the cross validation done on the training dataset. 
}
\usage{
assign.cv.output(processed.data, mcmc.pos.mean.trainingData, trainingData, 
trainingLabel, adaptive_B=FALSE, adaptive_S=FALSE, mixture_beta=TRUE, outputDir)
}

\arguments{
  \item{processed.data}{
The list object returned from the assign.preprocess function.
}
  \item{mcmc.pos.mean.trainingData}{
The list object returned from the assign.mcmc function. Notice that for cross validation, the Y argument in the assign.mcmc function should be set as the training dataset.
}
  \item{trainingData}{
The genomic measure matrix of training samples (i.g., gene expression matrix). The dimension of this matrix is probe number x sample number. The default is NULL.
}
  \item{trainingLabel}{
The list linking the index of each training sample to a specific group it belongs to. 
}
  \item{adaptive_B}{
Logicals. If TRUE, the model adapts the baseline/background (B) of genomic measures for the test samples. The default is FALSE. 
}
  \item{adaptive_S}{
Logicals. If TRUE, the model adapts the signatures (S) of genomic measures for the test samples. The default is FALSE. 
}
  \item{mixture_beta}{
Logicals. If TRUE, elements of the pathway activation matrix are modeled by a spike-and-slab mixuture distribution. The default is TRUE. 
}
  \item{outputDir}{
The path to the directory to save the output files. The path needs to be quoted in double quotation marks. 
}
}
\details{
The assign.cv.output function is suggested to run after the assign.preprocess, assign.mcmc and assign.summary function. For the cross validation, The Y argument in the assign.mcmc function is the output value "trainingData_sub" from the assign.preprocess function. 
}
\value{
The assign.cv.output returns one .csv file containing one/multiple pathway activity for each individual training samples, scatter plots of pathway activity for each individual pathway in all the training samples, and heatmap plots for the gene expression signatures for each individual pathways. 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Ying Shen 
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontshow{
setwd(tempdir())
tempdir <- tempdir()
data(trainingData1)
data(testData1)
data(geneList1)

trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, 
                                      ras=1:10, src=1:10), 
                       bcat = 11:19, e2f3 = 20:28, myc= 29:38, 
                       ras = 39:48, src = 49:55)
                       
processed.data <- assign.preprocess(trainingData=trainingData1, 
testData=testData1, trainingLabel=trainingLabel1, geneList=geneList1)
mcmc.chain <- assign.mcmc(Y=processed.data$trainingData_sub, 
Bg = processed.data$B_vector, X=processed.data$S_matrix, 
Delta_prior_p = processed.data$Pi_matrix, iter = 20, adaptive_B=TRUE, 
adaptive_S=FALSE, mixture_beta=TRUE)
mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=10, iter=20, 
adaptive_B=FALSE, adaptive_S=FALSE, mixture_beta=TRUE)
}
assign.cv.output(processed.data=processed.data, 
mcmc.pos.mean.trainingData=mcmc.pos.mean, trainingData=trainingData1,
trainingLabel=trainingLabel1, 
adaptive_B=FALSE, adaptive_S=FALSE, mixture_beta=TRUE, outputDir=tempdir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
