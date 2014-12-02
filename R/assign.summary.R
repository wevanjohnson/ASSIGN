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