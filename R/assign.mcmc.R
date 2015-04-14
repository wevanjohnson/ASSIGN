assign.mcmc <- function(Y, Bg, X, Delta_prior_p, iter=2000, adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE, sigma_sZero = 0.01, sigma_sNonZero = 1, p_beta = 0.01, sigma_bZero = 0.01, sigma_bNonZero = 1, alpha_tau = 1, beta_tau = 0.01, Bg_zeroPrior=TRUE, S_zeroPrior=TRUE, ECM = FALSE)
{  
  cat("Start Gibbs sampling...\n")
  
  Y <- as.matrix(Y)
  n <- NROW(Y)
  k <- NCOL(Y)
  S <- as.matrix(X)
  m <- NCOL(S)
  
  #prior of B
  if (adaptive_B == FALSE){
    mu_B_0 <- Bg
  } else {
    if (Bg_zeroPrior == TRUE){
      mu_B_0 <- rep(0,n)
    } else {
      mu_B_0 <- Bg
    }
  }
  
  s_B_0 <- rep(10^2, n)
  s_B_0_inv <- 1 / (s_B_0)
  sB0_inv_muB0 <-  s_B_0_inv * mu_B_0
  
  # prior of S 
  if (adaptive_S == FALSE){
    S_0 <- S
  } else {
    if (S_zeroPrior == TRUE){
      S_0 <- matrix(0,n,m)
    } else {
      S_0 <- S
    }
  }
  
  sigma_s1 <- sigma_sZero
  sigma_s2 <- sigma_sNonZero
  
  # prior of beta
  p_beta <- p_beta
  sigma_b1 <- sigma_bZero
  sigma_b2 <- sigma_bNonZero
  odds_beta <- (1-p_beta)/p_beta
  onesMK <- matrix(rep(1,m*k),m,k)
  
  # prior of delta
  p_delta <- Delta_prior_p
  odds <- (1-p_delta)/p_delta
  onesNM <- matrix(rep(1,n*m),n,m)
  
  
  # prior of tau
  u <- alpha_tau
  v <- beta_tau
  
  P_B <- matrix(nrow=iter, ncol=n)
  P_S <- array(dim=c(iter,n,m))
  P_delta <- array(dim=c(iter,n,m))
  P_delta_pr <- array(dim=c(iter,n,m))
  P_beta <- array(dim=c(iter, m,k))
  P_gamma <- array(dim=c(iter, m,k))
  P_kappa <- array(dim=c(iter, m,k))
  P_gamma_pr <- 	array(dim=c(iter, m,k))
  P_tau2 <- matrix(nrow=iter, ncol=n)
  
  #initial values
  if (adaptive_B == TRUE){
    B_temp <- rtnorm(n,0,1,lower=0)
  } else {
    B_temp <- mu_B_0
  }
  if (adaptive_S == TRUE){
    S_temp <- S 
  } else {
    S_temp <- S_0
  }
  delta_temp <- matrix(rbern(onesNM, p_delta), n, m)
  delta_pr_temp <- p_delta
  beta_temp <- matrix(0,m,k)
  if (mixture_beta == TRUE){
    gamma_temp <- matrix(rbern(onesMK,p_beta), m, k)
    gamma_pr_temp <- matrix(p_beta, m, k)
  } else {
    gamma_temp <- matrix(rep(1, m*k), m, k)
    gamma_pr_temp <- matrix(rep(1, m*k), m, k)
  }
  kappa_temp <- beta_temp * gamma_temp
  tau_temp <- rep(u/v, n)
  
  #update first row
  P_B[1, ] <- B_temp
  if (adaptive_S == TRUE){
    P_S[1, , ] <- S_temp
    P_delta[1, , ] <- delta_temp
    P_delta_pr[1, , ] <- delta_pr_temp
  }
  P_beta[1, , ] <- beta_temp
  if (mixture_beta == TRUE){
    P_gamma[1, , ] <- gamma_temp
    P_gamma_pr[1, , ] <- gamma_pr_temp
  }
  P_kappa[1, , ] <- kappa_temp
  P_tau2[1, ] <- tau_temp	
  
  for (i in 2:iter) {
    if(i %% 100 == 0){cat("iteration", i, "\n", sep=" ")}
    
    #update B
    if (adaptive_B == TRUE){
      tmp0 <- S_temp %*% beta_temp
      
      s_B_1 <- 1 / (k * tau_temp + s_B_0_inv)
      J <- apply(Y - tmp0, 1, sum)
      mu_B_1 <- s_B_1 * (tau_temp * J + sB0_inv_muB0)
      if(ECM == TRUE) {
        B_temp <- mu_B_1
      } else {
        B_temp <- rnorm(n, mu_B_1, sqrt(s_B_1))
      }
    }
    
    B_rep = matrix(rep(B_temp,k), n, k)
    Y_minus_B_rep <- Y - B_rep
    
    
    for (s in 1:m){
      if (m == 1){
        tmp1 <- S_temp[ ,s]^2*tau_temp
        tmp2 <- S_temp[ ,s]*tau_temp
        s_beta_0_inv <- 1 / ifelse(gamma_temp[s, ]==1, sigma_b2^2, sigma_b1^2)
        s_beta_1 <- 1/(sum(tmp1) + s_beta_0_inv)
        mu_beta_1 <- s_beta_1 * (colSums(tmp2*(Y-B_rep)))
      } else {
        tmp1 <- S_temp[ ,s]^2 * tau_temp	
        tmp2 <- S_temp[ ,s] * tau_temp	 
        tmp3 <- S_temp[,-s,drop=FALSE] %*% beta_temp[-s,,drop=FALSE]
        s_beta_0_inv <- 1 / ifelse(gamma_temp[s, ]==1, sigma_b2^2, sigma_b1^2)
        s_beta_1 <- 1/ (s_beta_0_inv + sum(tmp1))
        mu_beta_1 <- s_beta_1 * (t(tmp2) %*% (Y_minus_B_rep - tmp3))
      }
      if (ECM == TRUE) {
        beta_temp[s, ] <- mu_beta_1 
      } else {
        if (mixture_beta == TRUE){
          lower <- ifelse(gamma_temp[s, ]==1, 0, -Inf)
          upper <- ifelse(gamma_temp[s,]==1, 1, Inf)
          beta_temp[s, ] <- rtnorm(k, mu_beta_1, sqrt(s_beta_1), lower, upper)
          kappa_temp[s, ] <- beta_temp[s, ] * gamma_temp[s, ]
          
        } else {
          beta_temp[s, ] <- rtnorm(k, mu_beta_1, sqrt(s_beta_1),lower=0,upper=1) 
        }
      }
    }
    
    
    if (mixture_beta == TRUE){
      for (s in 1:m)
      {
        gamma_b_div_a <- (pnorm(1) - pnorm(0)) * odds_beta * (sigma_b2 / sigma_b1) * exp(-1/2 *(beta_temp[s,]^2/sigma_b1^2 - beta_temp[s,]^2/sigma_b2^2))
        gamma_pr_temp[s, ] <- ifelse(beta_temp[s, ] < 0,  0, 1/(1+gamma_b_div_a))
        gamma_temp[s, ] <- rbern(k, gamma_pr_temp[s, ])
      }
    }
    
    #update S
    if (adaptive_S == TRUE){
      #mu_S_0 <- ifelse(delta_temp==1, S_0, 0)
      #s_S_inv_0 <- 1 / ifelse(delta_temp==1, sigma_s2^2, sigma_s1^2)
      #tmp4 <- s_S_inv_0 * mu_S_0
      #s_S_inv_1 <- 1 / (as.matrix(tau_temp) %*% t(as.matrix(apply(beta_temp^2, 1, sum))) + s_S_inv_0)
      
      #for (j in 1:m){
      #  E1 <-  Y_minus_B_rep - S_temp[,-j,drop=FALSE] %*% beta_temp[-j,,drop=FALSE]
      #  mu_S_1 <- s_S_inv_1[, j] *(tau_temp*(E1 %*% beta_temp[j, ]) + tmp4[, j]) 
      #  S_temp[,j] <- rnorm(n, mu_S_1, sqrt(s_S_inv_1[, j]))
        
      #}
      S_temp <- S
      
      # update delta
      
      delta_b_div_a <- (sigma_s2 / sigma_s1) * odds * exp(-1/2 *(S_temp^2/sigma_s1^2 - (S_temp - S)^2/sigma_s2^2))
      delta_pr_temp <- 1/(1+delta_b_div_a)
      delta_temp <- matrix(rbern(onesNM, 1/(1+delta_b_div_a)), n, m)
    }
    
    # update tau
    tmp5 <- S_temp %*% beta_temp
    
    res <- Y_minus_B_rep - tmp5
    un <- u + k/2
    vn <- apply((res)^2, 1, sum)/2 + v
    tau_temp <- rgamma(n,un,vn)
    
    
    # update
    if (adaptive_B == TRUE){
      P_B[i, ] <- B_temp
    }
    if (adaptive_S == TRUE){
      P_S[i, , ] <- S_temp
      P_delta[i, , ] <- delta_temp
      P_delta_pr[i, , ] <- delta_pr_temp
    }
    P_beta[i, , ] <- beta_temp
    if (mixture_beta == TRUE){
      P_gamma[i, , ] <- gamma_temp
      P_gamma_pr[i, , ] <- gamma_pr_temp
    }
    P_kappa[i, , ] <- kappa_temp
    P_tau2[i, ] <- tau_temp
  }
  
  rtlist <- list(beta_mcmc = P_beta,  tau2_mcmc = P_tau2)
  if (adaptive_B == TRUE){
    rtlist <- c(rtlist, list(B_mcmc = P_B))
  }
  
  if (mixture_beta == TRUE){
    rtlist <- c(rtlist, list(gamma_mcmc = P_gamma, kappa_mcmc = P_kappa))
  }
  
  if (adaptive_S == TRUE){
    rtlist <- c(rtlist, list(S_mcmc = P_S, Delta_mcmc = P_delta))
  }
  
  return(rtlist)
  
}
