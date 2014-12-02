\name{assign.summary}
\alias{assign.summary}
\title{
Summary of the model parameters estimated by the Gibbs sampling algorithm
}
\description{
The assign.summary function computes the posterior mean of the model parameters estimated in every iteration during the Gibbs sampling. 
}
\usage{
assign.summary(test, burn_in=1000, iter=2000, adaptive_B = TRUE, adaptive_S = FALSE, 
mixture_beta = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{
The list object returned from the assign.mcmc function. The list components are the MCMC chains of the B, S, Delta, beta, gamma, and sigma.
}
  \item{burn_in}{
The number of burn-in iterations. These iterations are discarded when computing the posterior means of the model parameters. The default is 1000.
}
  \item{iter}{
The number of total iterations. The default is 2000.
}
  \item{adaptive_B}{
Logicals. If TRUE, the model adapts the baseline/background (B) of genomic measures for the test samples. The default is TRUE.}
  \item{adaptive_S}{
Logicals. If TRUE, the model adapts the signatures (S) of genomic measures for the test samples. The default is FALSE. 
}
 \item{mixture_beta}{
Logicals. If TRUE, elements of the pathway activation matrix are modeled by a spike-and-slab mixuture distribution. The default is TRUE.
}
}
\details{
The assign.summary function is suggested to run after the assign.convergence function, which is used to check the convergency of the MCMC chain. If the MCMC chain does not converge to a stationary phase, more iterations are required in the assign.mcmc function. The number of burn-in iterations is usually set to be half of the number of total iterations, meaning that the first half of the MCMC chain is discarded when computing the posterior means.
}
\value{
\item{beta_pos}{The N x K matrix of the posterior mean of the pathway activation level in test samples (transposed matrix A). Columns:K pathways; rows: N test samples}
\item{sigma_pos}{The G x 1 vector of the posterior mean of the variance of gene.}
\item{kappa_pos}{The N x K matrix of posterior mean of pathway activation level in test samples (transposed matrix A) (adjusted beta_pos scaling between 0 and 1). Columns:K pathways; rows: N test samples}
\item{gamma_pos}{The N x K matrix of the posterior probability of pathways being activated in test samples.}
\item{S_pos}{The G x K matrix of the posterior mean of pathway signature genes.}
\item{Delta_pos}{The G x K matrix of the posterior probability of genes being significant in the associated pathways.}
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

\examples{
\dontshow{
data(trainingData1)
data(testData1)
data(geneList1)
trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, 
                                      ras=1:10, src=1:10), 
                       bcat = 11:19, e2f3 = 20:28, myc= 29:38, 
                       ras = 39:48, src = 49:55)
                       
processed.data <- assign.preprocess(trainingData=trainingData1, 
testData=testData1, trainingLabel=trainingLabel1, geneList=geneList1)
mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, Bg = processed.data$B_vector, 
X=processed.data$S_matrix, Delta_prior_p = processed.data$Pi_matrix, iter = 20, 
adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE)
mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=10, iter=20, 
adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
