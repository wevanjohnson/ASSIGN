#' ASSIGN All-in-one function
#' 
#' The assign.wrapper function integrates the assign.preprocess, assign.mcmc,
#' assign.summary, assign.output, assign.cv.output functions into one wrapper
#' function.
#' 
#' The assign.wrapper function is an all-in-one function which outputs the
#' necessary results for basic users. For users who need more
#' intermediate results for model diagnosis, it is better to run the
#' assign.preprocess, assign.mcmc, assign.convergence, assign.summary functions
#' separately and extract the output values from the returned list objects of
#' those functions.
#' 
#' @param trainingData The genomic measure matrix of training samples (i.g.,
#' gene expression matrix). The dimension of this matrix is probe number x
#' sample number. The default is NULL.
#' @param testData The genomic measure matrix of test samples (i.g., gene
#' expression matrix). The dimension of this matrix is probe number x sample
#' number.
#' @param trainingLabel The list linking the index of each training sample to a
#' specific group it belongs to. See examples for more information.
#' @param testLabel The vector of the phenotypes/labels of the test samples.
#' The default is NULL.
#' @param geneList The list that collects the signature genes of one/multiple
#' pathways. Every component of this list contains the signature genes
#' associated with one pathway. The default is NULL.
#' @param anchorGenes A list of genes that will be included in the signature
#' even if they are not chosen during gene selection.
#' @param excludeGenes A list of genes that will be excluded from the signature
#' even if they are chosen during gene selection.
#' @param n_sigGene The vector of the signature genes to be identified for one
#' pathway. n_sigGene needs to be specified when geneList is set NULL. The
#' default is NA. See examples for more information.
#' @param adaptive_B Logicals. If TRUE, the model adapts the
#' baseline/background (B) of genomic measures for the test samples. The
#' default is TRUE.
#' @param adaptive_S Logicals. If TRUE, the model adapts the signatures (S) of
#' genomic measures for the test samples. The default is FALSE.
#' @param mixture_beta Logicals. If TRUE, elements of the pathway activation
#' matrix are modeled by a spike-and-slab mixuture distribution. The default is
#' TRUE.
#' @param outputDir The path to the directory to save the output files. The
#' path needs to be quoted in double quotation marks.
#' @param p_beta p_beta is the prior probability of a pathway being activated
#' in individual test samples. The default is 0.01.
#' @param theta0 The prior probability for a gene to be significant, given that
#' the gene is NOT defined as "significant" in the signature gene lists
#' provided by the user. The default is 0.05.
#' @param theta1 The prior probability for a gene to be significant, given that
#' the gene is defined as "significant" in the signature gene lists provided by
#' the user. The default is 0.9.
#' @param iter The number of iterations in the MCMC. The default is 2000.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters. The
#' default is 1000.
#' @param sigma_sZero Each element of the signature matrix (S) is modeled by a
#' spike-and-slab mixuture distribution. Sigma_sZero is the variance of the
#' spike normal distribution. The default is 0.01.
#' @param sigma_sNonZero Each element of the signature matrix (S) is modeled by
#' a spike-and-slab mixuture distribution. Sigma_sNonZero is the variance of
#' the slab normal distribution. The default is 1.
#' @param S_zeroPrior Logicals. If TRUE, the prior distritribution of signature
#' follows a normal distribution with mean zero. The default is TRUE.
#' @param pctUp By default, ASSIGN bayesian gene selection chooses the
#' signature genes with an equal fraction of genes that increase with pathway
#' activity and genes that decrease with pathway activity. Use the pctUp
#' parameter to modify this fraction. Set pctUP to NULL to select the most
#' significant genes, regardless of direction. The default is 0.5
#' When running ASSIGN with a balanced signature, specifies the
#' percent of the genes in the signature that increase with pathway activity.
#' The default is 0.5.
#' @param geneselect_iter The number of iterations for bayesian gene selection. The
#' default is 500.
#' @param geneselect_burn_in The number of burn-in iterations for bayesian gene selection.
#' The default is 100
#' @return The assign.wrapper returns one/multiple pathway activity for each
#' individual training sample and test sample, scatter plots of pathway
#' activity for each individual pathway in the training and test data,
#' heatmap plots for gene expression signatures for each individual pathway,
#' heatmap plots for the gene expression of the prior and posterior
#' signtures (if adaptive_S equals TRUE) of each individual pathway in the test
#' data
#' @author Ying Shen and W. Evan Johnson
#' @examples
#' 
#' \dontshow{
#' setwd(tempdir())
#' tempdir <- tempdir()
#' }
#' data(trainingData1)
#' data(testData1)
#' data(geneList1)
#' 
#' trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10, ras=1:10,
#' src=1:10), bcat = 11:19, e2f3 = 20:28, myc= 29:38, ras = 39:48, src = 49:55)
#' testLabel1 <- rep(c("subtypeA","subtypeB"),c(53,58))
#' 
#' assign.wrapper(trainingData=trainingData1, testData=testData1,
#' trainingLabel=trainingLabel1, testLabel=testLabel1, geneList=geneList1,
#' adaptive_B=TRUE, adaptive_S=FALSE, mixture_beta=TRUE,
#' outputDir=tempdir, p_beta=0.01, theta0=0.05, theta1=0.9,
#' iter=20, burn_in=10)
#' 
#' @export assign.wrapper
assign.wrapper<-function (trainingData = NULL, testData, trainingLabel, testLabel = NULL,
          geneList = NULL, anchorGenes = NULL, excludeGenes = NULL, n_sigGene = NA,
          adaptive_B = TRUE, adaptive_S = FALSE, mixture_beta = TRUE, outputDir, p_beta = 0.01,
          theta0 = 0.05, theta1 = 0.9, iter = 2000, burn_in = 1000, sigma_sZero = 0.01,
          sigma_sNonZero = 1, S_zeroPrior=FALSE, pctUp=0.5, geneselect_iter=500, geneselect_burn_in=100)
{
  if (is.null(geneList)) {
    pathName <- names(trainingLabel)[-1]
  }
  else {
    pathName <- names(geneList)
  }
  processed.data <- assign.preprocess(trainingData, testData, anchorGenes, excludeGenes,
                                      trainingLabel, geneList, n_sigGene, theta0, theta1, pctUp=pctUp,
                                      geneselect_iter=geneselect_iter, geneselect_burn_in=geneselect_burn_in)
  if (!is.null(trainingData)) {
    cat("Estimating model parameters in the training dataset...\n")
    mcmc.chain.trainingData <- assign.mcmc(Y = processed.data$trainingData_sub,
                                           Bg = processed.data$B_vector, X = processed.data$S_matrix,
                                           Delta_prior_p = processed.data$Pi_matrix, iter = iter,
                                           sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, S_zeroPrior=S_zeroPrior,
                                           adaptive_B = FALSE, adaptive_S = FALSE, mixture_beta = TRUE)
    mcmc.pos.mean.trainingData <- assign.summary(test = mcmc.chain.trainingData, 
                                                 burn_in = burn_in, iter = iter, adaptive_B = FALSE, 
                                                 adaptive_S = FALSE, mixture_beta = TRUE)

  }
  cat("Estimating model parameters in the test dataset...\n")
  mcmc.chain.testData <- assign.mcmc(Y = processed.data$testData_sub, 
                                     Bg = processed.data$B_vector, X = processed.data$S_matrix,
                                     Delta_prior_p = processed.data$Pi_matrix, iter = iter,
                                     sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, S_zeroPrior=S_zeroPrior,
                                     adaptive_B = adaptive_B, adaptive_S = adaptive_S, mixture_beta = mixture_beta,
                                     p_beta = p_beta)
  mcmc.pos.mean.testData <- assign.summary(test = mcmc.chain.testData, 
                                           burn_in = burn_in, iter = iter, adaptive_B = adaptive_B, 
                                           adaptive_S = adaptive_S, mixture_beta = mixture_beta)

  cat("Outputing results...\n")
  if (mixture_beta) {
    if (!is.null(trainingData)) {
      coef_train = mcmc.pos.mean.trainingData$kappa_pos
    }
    coef_test = mcmc.pos.mean.testData$kappa_pos
  }
  else {
    if (!is.null(trainingData)) {
      coef_train = mcmc.pos.mean.trainingData$beta_pos
    }
    coef_test = mcmc.pos.mean.testData$beta_pos
  }
  cwd <- getwd()
  dir.create(outputDir,showWarnings = F)##moom added this to create the output folder if doesn't exist already.
  setwd(outputDir)
  param=as.matrix(paste(pathName,"analysis was run using the following parameters :",
        "n_sigGene=",n_sigGene, "adaptive_B=",adaptive_B,"adaptive_S=", adaptive_S,
        "mixture_beta=",mixture_beta,"p_beta=",p_beta,"theta0=", theta0, "theta1=",theta1, 
        "iter=",iter, "burn_in=",burn_in,"The output files are located at:",outputDir,sep=' '))###moom added this 
  write.table(param,"parameters.txt",col.names=F,sep='\t')###moom added this 
  if (!is.null(trainingData)) {
    rownames(coef_train) <- colnames(processed.data$trainingData_sub)
    colnames(coef_train) <- pathName
    write.csv(processed.data$S_matrix, file = "signature_gene_list_prior.csv")###moom added this to include the gene list and prior coefficient
    write.csv(coef_train, file = "pathway_activity_trainingset.csv")
  }
  rownames(coef_test) <- colnames(processed.data$testData_sub)
  colnames(coef_test) <- pathName
  write.csv(coef_test, file = "pathway_activity_testset.csv")
  if (!is.null(trainingData)) {
    heatmap.train(diffGeneList = processed.data$diffGeneList, 
                  trainingData, trainingLabel)
  }
  heatmap.test.prior(diffGeneList = processed.data$diffGeneList, 
                     testData, trainingLabel, testLabel, coef_test, geneList)
  if (adaptive_S) {
    heatmap.test.pos(testData = processed.data$testData_sub, 
                     Delta_pos = mcmc.pos.mean.testData$Delta_pos, trainingLabel, 
                     testLabel, Delta_cutoff = 0.95, coef_test, geneList)
    ####Added by moom####
    ##Evan please double check if this is informative and needed
    pdf("Signature_convergence.pdf") 
    plot(mcmc.chain.testData$S_mcmc)
    abline(h=0,col="red")
    dev.off()
    ##
    dimnames(mcmc.pos.mean.testData$Delta_pos)=dimnames(processed.data$S_matrix)    
    deltas<-cbind(processed.data$S_matrix,processed.data$Delta_matrix,mcmc.pos.mean.testData$S_pos,mcmc.pos.mean.testData$Delta_pos)
    colnames(deltas)=c(paste("Prior change in expression",pathName,sep=":"),paste("Prior probability of inclusion",pathName,sep=":"),paste("Posterior change in expression",pathName,sep=":"),paste("Posterior probability of inclusion",pathName,sep=":"))
    delta_in=NULL
    for(i in 1:ncol(deltas)){delta_in[i]=(strsplit(colnames(deltas),":")[[i]][2])}
    write.csv(round(deltas[,order(delta_in)],digits = 4),"posterior_delta.csv",quote=F)
    #####End: added by moom###
      }
  if (!is.null(trainingData)) {
    scatter.plot.train(coef_train, trainingData, trainingLabel)
  }
  scatter.plot.test(coef_test, trainingLabel, testLabel, geneList)
  if (!is.null(testLabel)) {
    box.plot.test(coef_test, trainingLabel, testLabel, geneList)
  }
  if (!is.null(trainingData)) {
    output.data <- list(processed.data = processed.data, 
                        mcmc.pos.mean.trainingData = mcmc.pos.mean.trainingData, 
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  else {
    output.data <- list(processed.data = processed.data, 
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  save(output.data, file = "output.rda")
  setwd(cwd)
}
