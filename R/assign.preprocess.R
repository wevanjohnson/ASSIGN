#' Input data preprocessing
#' 
#' The assign.preprocess function is used to perform quality control on the
#' user-provided input data and generate starting values and/or prior values
#' for the model parameters. The assign.preprocess function is optional. For
#' users who already have the correct format for the input of the assign
#' function, they can skip this step and go directly to the assign.mcmc
#' function.
#' 
#' The assign.preprocess is applied to perform quality control on the
#' user-provided genomic data and meta data, re-format the data in a way that
#' can be used in the following analysis, and generate starting/prior values
#' for the pathway signature matrix. The output values of the assign.preprocess
#' function will be used as input values for the assign.mcmc function.
#' 
#' For training data with 1 control group and 3 experimental groups (10
#' samples/group; all 3 experimental groups share 1 control group), the
#' trainingLabel can be specified as: trainingLabel <- list(control =
#' list(expr1=1:10, expr2=1:10, expr3=1:10), expr1 = 11:20, expr2 = 21:30,
#' expr3 = 31:40)
#' 
#' For training data with 3 control groups and 3 experimental groups (10
#' samples/group; Each experimental group has its corresponding control group),
#' the trainingLabel can be specified as: trainingLabel <- list(control =
#' list(expr1=1:10, expr2=21:30, expr3=41:50), expr1 = 11:20, expr2 = 31:40,
#' expr3 = 51:60)
#' 
#' It is highly recommended that the user use the same expriment name when
#' specifying control indice and exprimental indice.
#' 
#' @param trainingData The genomic measure matrix of training samples (i.g.,
#' gene expression matrix). The dimension of this matrix is probe number x
#' sample number. The default is NULL.
#' @param testData The genomic measure matrix of test samples (i.g., gene
#' expression matrix). The dimension of this matrix is probe number x sample
#' number.
#' @param trainingLabel The list linking the index of each training sample to a
#' specific group it belongs to. See details and examples for more information.
#' @param geneList The list that collects the signature genes of one/multiple
#' pathways. Every component of this list contains the signature genes
#' associated with one pathway. The default is NULL.
#' @param n_sigGene The vector of the signature genes to be identified for one
#' pathway. n_sigGene needs to be specified when geneList is set NULL. The
#' default is NA. See examples for more information.
#' @param theta0 The prior probability for a gene to be significant, given that
#' the gene is NOT defined as "significant" in the signature gene lists
#' provided by the user. The default is 0.05.
#' @param theta1 The prior probability for a gene to be significant, given that
#' the gene is defined as "significant" in the signature gene lists provided by
#' the user. The default is 0.9.
#' @return \item{trainingData_sub}{The G x N matrix of G genomic measures
#' (i.g., gene expession) of N training samples. Genes/probes present in at
#' least one pathway signature are retained. Only returned when the training
#' dataset is available.} \item{testData_sub}{The G x N matrix of G genomic
#' measures (i.g., gene expession) of N test samples. Genes/probes present in
#' at least one pathway signature are retained.} \item{B_vector}{The G x 1
#' vector of genomic measures of the baseline/background. Each element of the
#' B_vector is calculated as the mean of the genomic measures of the control
#' samples in training data.} \item{S_matrix}{The G x K matrix of genomic
#' measures of the signature. Each column of the S_matrix represents a pathway.
#' Each element of the S_matrix is calculated as the mean of genomic measures
#' of the experimental samples minus the mean of the control samples in the
#' training data.} \item{Delta_matrix}{The G x K matrix of binary indicators.
#' Each column of the Delta_matrix represents a pathway. The elements in
#' Delta_matrix are binary (0, insignificant gene; 1, significant gene).}
#' \item{Pi_matrix}{The G x K matrix of probability p of a Bernoulli
#' distribution. Each column of the Pi_matrix represents a pathway. Each
#' element in the Pi_matrix is the probability of a gene to be significant in
#' its associated pathway. } \item{diffGeneList}{The list that collects the
#' signature genes of one/multiple pathways generated from the training samples
#' or from the user provided gene list. Every component of this list contains
#' the signature genes associated with one pathway. }
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
#' }                       
#' processed.data <- assign.preprocess(trainingData=trainingData1, 
#' testData=testData1, trainingLabel=trainingLabel1, geneList=geneList1)
#' 
#' @export assign.preprocess
assign.preprocess <- function(trainingData=NULL, testData, anchorGenes=NULL, excludeGenes=NULL,
                              trainingLabel, geneList=NULL, n_sigGene=NA, theta0=0.05, theta1=0.9, balanced=FALSE, pctUp=0.5){
  cat("Runing ASSIGN development version: truncated_S\n")
  cat("Performing QC on the input data...\n")
  dat <- qc(trainingData, testData, geneList)
  if (!is.null(trainingLabel)){
    if (identical(names(trainingLabel$control), names(trainingLabel)[-1]) == FALSE){warning("Control Labels DO NOT match the experimental Labels!\nPlease make sure that you specify the correct indice for control and experimental samples in the trainingLabel!")}
  }
  
  cat("Generating starting/prior values for model parameters...\n")
  if (is.null(trainingData) & is.null(geneList)){
    stop("trainingData and geneList are both set NULL. Need one of them for the analysis!")
  }
  if (!is.null(trainingData) & !is.null(trainingLabel)){
    x <- data_prep_s1(n_sigGene, trainingData=dat$trainingData, testData=dat$testData,
                      trainingLabel, anchorGenes, excludeGenes, geneList=dat$geneList,
                      theta0, theta1, balanced=balanced, pctUp=pctUp)
    return(list(trainingData_sub=x$trainingData_sub, testData_sub=x$testData_sub,
                B_vector=x$trainingBaseline_sub, S_matrix=x$S_matrix,
                Delta_matrix=x$Delta_matrix, Pi_matrix=x$Pi_matrix, diffGeneList=x$diffGeneList))
  }
  else if (is.null(trainingData) & is.null(trainingLabel) & !is.null(geneList)) {
    x <- data_prep_s2(geneList=dat$geneList, anchorGenes, testData=dat$testData, theta0, theta1)
    return(list(testData_sub=x$testData_sub, B_vector=x$B_vector, S_matrix=x$S_matrix,
                Delta_matrix=x$Delta_matrix, Pi_matrix=x$Pi_matrix, diffGeneList=dat$geneList))
  }
}
