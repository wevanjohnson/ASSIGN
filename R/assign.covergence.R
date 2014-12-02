assign.convergence <- function(test, burn_in=0, iter=2000, parameter=c("B","S","Delta", "beta", "kappa", "gamma", "sigma"), whichGene, whichSample, whichPath){
  if (parameter == "B"){
    B_pos <- test$B_mcmc[burn_in:iter, whichGene]
    plot(B_pos, type="l", xlab="Iteration", ylab="Posterior B", main="Convergency of posterior B")
    return(B_pos)
  }
  
  if (parameter == "S"){
    S_pos <- test$S_mcmc[burn_in:iter, whichGene, whichPath]
    plot(S_pos, type="l", xlab="Iteration", ylab=paste("Posterior S of pathway ", whichPath, " gene ", whichGene, sep=""), main="Convergency of posterior S")
    return(S_pos)
  }
  
  if (parameter == "Delta"){
    Delta_pos <- test$Delta_mcmc[burn_in:iter, whichGene, whichPath]
    plot(Delta_pos, type="l", xlab="Iteration", ylab=paste("Posterior Delta of pathway ", whichPath, " gene ", whichGene, sep=""), main="Convergency of posterior Delta")
    return(Delta_pos)
  }
  
  
  if (parameter == "beta"){
    beta_pos <- test$beta_mcmc[burn_in:iter,whichPath,whichSample]
    plot(beta_pos, type="l", xlab="Iteration", ylab=paste("Posterior beta of sample ", whichSample, " pathway ", whichPath, sep=""), main="Convergency of posterior beta")
    return(beta_pos)
  }
  
  if (parameter == "kappa"){
    kappa_pos <- test$kappa_mcmc[burn_in:iter,whichPath,whichSample]
    plot(kappa_pos, type="l", xlab="Iteration", ylab=paste("Posterior kappa of sample ", whichSample, " pathway ", whichPath, sep=""), main="Convergency of posterior kappa")
    return(kappa_pos)
  }
  
  if (parameter == "gamma"){
    gamma_pos <- test$gamma_mcmc[burn_in:iter,whichPath,whichSample]
    plot(gamma_pos, type="l", xlab="Iteration", ylab=paste("Posterior gamma of sample ", whichSample, " pathway ", whichPath, sep=""), main="Convergency of posterior gamma")
    return(gamma_pos)
  }
  
  if (parameter == "sigma"){
    sigma_pos <- 1/test$tau2_mcmc[burn_in:iter, whichGene]
    plot(sigma_pos, type="l",xlab="Iteration", ylab=paste("Posterior sigma^2 of gene ", whichGene, sep=""), main="Convergency of posterior sigma^2")
    return(sigma_pos)
  }
}