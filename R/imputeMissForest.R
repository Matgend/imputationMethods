# missForest
############
#' @title non-parametric missing value imputation for mixed-type data by missForest(randomForest)
#' @description This function imputes missing data for continuous and discrete traits applying the missForest
#' approach  
#' @usage imputeMissForest(missingData, variance_fraction = 0.8, maxiter = 10, ntree = 100, 
#' mtry = sqrt(ncol(missingData)), tree = NULL)
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors
#' @param maxiter maximum number of iterations to be performed given the stopping criterion is not met beforehand.
#' @param ntree number of trees to grow in each forest.
#' @param mtry number of variables randomly sampled at each split. By default it's the square root of the number of 
#' variables
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond 
#' to the 
#' @param tree phylo object
#' phylogenetic inertia
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values. + parameters used for the imputation
#' @export
imputeMissForest <- function(missingData, variance_fraction = 0, maxiter = 10, ntree = 100, 
                             mtry = sqrt(ncol(missingData)), tree = NULL, hint = NULL){
  
  Nvariables <- ncol(missingData)
  
  #include imputed data in case 
  if(!is.null(hint)){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))
    
    missingData <- cbind(missingData, hint)
  }
  
  # want to include phylogeny information
  if(variance_fraction != 0 & variance_fraction != 2){
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], 
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }
  
  #run missForest
  missForest_imputation <- missForest::missForest(xmis = missingData, 
                                                  maxiter = maxiter, ntree = ntree, mtry = mtry, verbose = TRUE)
  
  #cut eigenvectors columns
  missForest_imputation$ximp <- missForest_imputation$ximp[,1:Nvariables, drop = FALSE]
  
  parameters <- list(maxiter = maxiter, ntree = ntree, mtry = mtry)
  
  return(list(imputedData = missForest_imputation$ximp, parameters = parameters))
}
