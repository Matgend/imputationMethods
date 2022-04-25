#Categorical in dummy
#####################
#' @title Conversion factors in dummy
#' @description This function returns a data.frame of dummy corresponding to a data.frame of factors      
#' representing a categorical variable.
#' @usage generateDummyVariables(NaNData)
#' @param NaNData data.frame of one ore several factors columns 
#' @return a data.frame in which each variable are represented as dummies.
generateDummyVariables <- function(NaNData){
  
  Nstates <- as.numeric(tail(levels(NaNData[,1]), n = 1)) + 1
  
  #in case NaNData is a vector
  if(is.null(dim(NaNData))){
    stop("The column(s) shoud be of class data.frame")
  }
  
  #check if all the columns are factors
  if(sum(sapply(NaNData, is.factor)) != ncol(NaNData)){
    stop("The columns should be of class factor")
  }
  
  dummy <- fastDummies::dummy_cols(NaNData)
  row.names(dummy) <- row.names(NaNData)
  columns_to_remove <- c(1:ncol(NaNData), grep("NA", colnames(dummy)))
  dummy <- dummy[, -columns_to_remove]
  
  if(is.integer(dummy) & Nstates != 1){
    stateInData <- as.numeric(levels(NaNData[,1]))+1
    intermediateDummy <- dummy
    intermediateDummy[which(intermediateDummy == 1)] <- 0
    intermediateMatrix <- matrix(rep(intermediateDummy, Nstates), ncol = Nstates)
    intermediateMatrix[, stateInData] <- dummy
    dummy <- intermediateMatrix
    row.names(dummy) <- row.names(NaNData)
  }
  
  return(dummy)
  
}

# Generating eigenvectors
#########################
#' @title Eigenvectors calculations
#' @description This function calculates the eigenvectors of a phylo object 
#' @usage get_eigenvec(tree, variance_fraction = 0.8, numEigen = NULL)
#' @param tree phylogenetic tree of class "phylo"
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond
#' to the phylogenetic inertia
#' @return a data.frame in which each colum represent an eigenvector
get_eigenvec <- function(tree, variance_fraction){
  
  decomp <- PVR::PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
  eigvec <- as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors
  
  egval <- decomp@Eigen$values #extract eigenvalues
  eigPerc <- egval/(sum(egval)) #calculate % of variance
  eigPercCum <- t(cumsum(eigPerc)) #cumulated variance
  
  #eigenvectors representing more than X% variance
  numEigen <- sum(eigPercCum < variance_fraction)
  
  eigOK <- eigvec[,1:numEigen, drop = T] 
  
  if(numEigen == 0){
    print(paste("The variance_fraction should at leat be equal to ", 
                eigPercCum[1]))
    eigOK <- eigvec[1]
  }
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  #Eigenvectors generated in object 'eigenTobind'
  #rename eigenTobind species column so it matches trait dataset species column
  eigOK <- as.data.frame(eigOK)
  names(eigOK)[1] <- "c1"
  row.names(eigOK) <- decomp@phylo$tip.label
  
  return(eigOK)
}


#' @title Check + conversion columns class. 
#' @description This function checks if character and integer columns are factors and columns composed of float are
#' numeric. If not the case the columns are converted. 
#' @usage checkConvert(missingData)
#' @param missingData array or dataframe of data with missing values
#' @return return the missing dataset with the columns containing integers or characters as factors and columns 
#' containing float as numeric.
#' @export
checkConvert <- function(missingData){
  
  for(c in 1:ncol(missingData)){
    
    if(class(missingData[ ,c]) == "character" | class(missingData[,c]) == "integer"){
      missingData[, c] <- as.factor(missingData[, c])
    }
    
    else if(class(missingData[,c]) == "numeric"){
      
      #check if composed only of integers 
      nbrIntegers <- (sum(missingData[, c] - floor(missingData[,c]) == 0, na.rm = T))
      
      if(nbrIntegers == length(missingData[,c])){
        missingData[, c] <- as.factor(missingData[, c])
      }
    }
  }
  return(missingData)
}
