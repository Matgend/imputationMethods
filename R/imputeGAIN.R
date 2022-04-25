#GAIN
#####
#' @title GAIN
#' @description This function imputes missing data for continuous and discrete traits applying the Generative 
#' Adversarial #' Imputation Networks (GAIN)
#' @usage gainR(missingData, variance_fraction, Data, batch_size, hint_rate = 0.9, alpha = 100, epochs = 10000)
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond 
#'  to the phylogenetic inertia
#' @param tree phylo object
#' @param batch_size integer
#' @param hint_rate numerical
#' @param alpha numerical, hyperparameter
#' @param epochs integer, iterations
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a list of list containing in the "tab" imputedData the imputed Data and in the "tab" parametersGain, the discriminant loss values, the generative loss values, the MSE loss values and the iterations correspond to these values (important to generate a plot).
gainR <- function(missingData, variance_fraction, tree = NULL, batch_size = round(ncol(missingData)*0.2),
                  hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL){

  NbrCol <- ncol(missingData)

  #get the factors columns
  factorsColumns <- names(Filter(is.factor, missingData))

  rNames <- row.names(missingData)
  colNames <- colnames(missingData)

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

  #convert all the values as numeric
  missingData <- apply(missingData, 2, as.numeric)

  #source python function
  
  reticulate::source_python(system.file("python",
                            "mainR.py",
                            package = "GAIN"))
  
  #run Gain
  obj <- Gain()
  gainOuput <- obj$runGain(missingData, batch_size, hint_rate, alpha, epochs)

  #merge list of list
  for (l in 2:length(gainOuput)){
    gainOuput[[l]] <- unlist(gainOuput[[l]], recursive = FALSE)
  }

  #convert imputedData in dataframe
  imputedData <- as.data.frame(gainOuput[[1]])
  imputedData <- imputedData[, 1:NbrCol, drop = FALSE]
  names(imputedData) <- colNames
  row.names(imputedData) <- rNames

  #convert the column that should be factors
  imputedData[ ,factorsColumns] <- lapply(imputedData[ ,factorsColumns], factor)

  #parameters
  lossValues <- gainOuput[c(2:5)]
  names(lossValues) <- c("d_loss", "g_loss", "mse_loss", "epochs")

  return(list(imputedData = imputedData, parameters = lossValues))
}