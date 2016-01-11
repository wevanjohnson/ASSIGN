#' Summary of the model parameters estimated by the Gibbs sampling algorithm
#' 
#' The assign.summary function computes the posterior mean of the model
#' parameters estimated in every iteration during the Gibbs sampling.
#' 
#' The assign.summary function is suggested to run after the assign.convergence
#' function, which is used to check the convergency of the MCMC chain. If the
#' MCMC chain does not converge to a stationary phase, more iterations are
#' required in the assign.mcmc function. The number of burn-in iterations is
#' usually set to be half of the number of total iterations, meaning that the
#' first half of the MCMC chain is discarded when computing the posterior
#' means.
#' 
#' @param test The list object returned from the assign.mcmc function. The list
#' components are the MCMC chains of the B, S, Delta, beta, gamma, and sigma.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters. The
#' default is 1000.
#' @param iter The number of total iterations. The default is 2000.
#' @param adaptive_B Logicals. If TRUE, the model adapts the
#' baseline/background (B) of genomic measures for the test samples. The
#' default is TRUE.
#' @param adaptive_S Logicals. If TRUE, the model adapts the signatures (S) of
#' genomic measures for the test samples. The default is FALSE.
#' @param mixture_beta Logicals. If TRUE, elements of the pathway activation
#' matrix are modeled by a spike-and-slab mixuture distribution. The default is
#' TRUE.
#' @return \item{beta_pos}{The N x K matrix of the posterior mean of the
#' pathway activation level in test samples (transposed matrix A). Columns:K
#' pathways; rows: N test samples} \item{sigma_pos}{The G x 1 vector of the
#' posterior mean of the variance of gene.} \item{kappa_pos}{The N x K matrix
#' of posterior mean of pathway activation level in test samples (transposed
#' matrix A) (adjusted beta_pos scaling between 0 and 1). Columns:K pathways;
#' rows: N test samples} \item{gamma_pos}{The N x K matrix of the posterior
#' probability of pathways being activated in test samples.} \item{S_pos}{The G
#' x K matrix of the posterior mean of pathway signature genes.}
#' \item{Delta_pos}{The G x K matrix of the posterior probability of genes
#' being significant in the associated pathways.}
#' @author Ying Shen
#' @examples
#' 
#' \dontshow{
#' data(trainingData1)
#' data(testData1)
#' data(geneList1)
#' trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, 
#'                                       ras=1:10, src=1:10), 
#'                        bcat = 11:19, e2f3 = 20:28, myc= 29:38, 
#'                        ras = 39:48, src = 49:55)
#'                        
#' processed.data <- assign.preprocess(trainingData=trainingData1, 
#' testData=testData1, trainingLabel=trainingLabel1, geneList=geneList1)
#' mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, Bg = processed.data$B_vector, 
#' X=processed.data$S_matrix, Delta_prior_p = processed.data$Pi_matrix, iter = 20, 
#' adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE)
#' mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=10, iter=20, 
#' adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta = TRUE)
#' }
#' 
#' @export assign.summary
assign.summary <- function(test, burn_in=1000, iter=2000, adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE){
  beta_pos <- test$beta_mcmc[burn_in:iter, , ,drop=FALSE]
  beta_pos_mean <- apply(beta_pos, 3:2, mean)
  
  sigma_pos <- sqrt(1 / test$tau2_mcmc[burn_in:iter, ])
  sigma_pos_mean <- apply(sigma_pos, 2, mean)
  
  rtlist <- list(beta_pos = beta_pos_mean, sigma_pos=sigma_pos_mean)
  
  if (adaptive_B == TRUE) {
    B_pos <- test$B_mcmc[burn_in:iter, ]
    B_pos_mean <- apply(B_pos, 2, mean)
    
    rtlist <- c(rtlist, list(B_pos=B_pos_mean))
  }
  
  
  if (mixture_beta == TRUE) {
    kappa_pos <- test$kappa_mcmc[burn_in:iter, , ,drop=FALSE]
    kappa_pos_mean <- apply(kappa_pos, 3:2, mean)
    
    gamma_pos <- test$gamma_mcmc[burn_in:iter, , ,drop=FALSE]
    gamma_pos_mean <- apply(gamma_pos, 3:2, mean)
    
    rtlist <- c(rtlist, list(gamma_pos=gamma_pos_mean, kappa_pos=kappa_pos_mean))
  }
  
  if (adaptive_S == TRUE) {
    S_pos <- test$S_mcmc[burn_in:iter, , ,drop=FALSE]
    S_pos_mean <- apply(S_pos, 2:3, mean)
    
    Delta_pos <- test$Delta_mcmc[burn_in:iter, , ,drop=FALSE]
    Delta_pos_mean <- apply(Delta_pos, 2:3, mean)
    
    rtlist <- c(rtlist, list(S_pos=S_pos_mean, Delta_pos=Delta_pos_mean))
  }
  
  return(rtlist)	
}
