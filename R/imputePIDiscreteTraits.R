#' @title Imputation of missing data for one discrete trait
#' 
#' @description This function imputes missing data for discrete traits using the R package phytools. The first step
#' is to run the function make.simmap fits a continuous-time reversible Markov model for the evolution of one trait
#' and then simulates stochastic character histories using that model and the tip states on the tree. 
#' for three models(ER, SYM and ARD) and select the one having the smallest AIC.
#' @usage imputeOneDiscreteTrait(trait, Data)
#' @param missingData data.frame of 1 factor column containing NAs
#' @param tree phylo
#' @return a data.frame of 1 factor column with the NAs replaced by values
#' @export
imputeOneDiscreteTrait <- function(missingData, tree){
  
  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(tree$tip.label, row.names(missingData))){
    
    #change order of the rows, match the order of the phylogeny
    missingData <- missingData[match(tree$tip.label, row.names(missingData)), drop = FALSE]
  }
  
  colName <- names(missingData)
  #print(missingData)
  #if only one state represented in the trait
  if(sum(!is.na(unique(missingData[, 1]))) == 1){
    #print("one state")
    state <- missingData[which(!is.na(missingData)), ][1]
    missingData[which(is.na(missingData)), ] <- state
    return(list(imputedData = missingData, parameters = list("noModel")))
  }
  
  # #if trait is ordinal (use random selection imputation?) (now use continuous method)
  # row <- as.numeric(str_extract(colnames(missingData), "(?<=\\/)\\d+"))
  # if(Data$dataframe[row,"class"] == "ordinal"){
  #   print("ordinal trait")
  #   missingData[, 1] <- as.numeric(missingData[, 1]) - 1
  #   #print(missingData)
  #   missingData <- imputeContinousTraits(missingData)
  #   #print(missingData)
  #   missingData[, 1] <- round(missingData[, 1], 0)
  #   return(list(imputedData = missingData, parameters = list("noModel")))
  # }
  
  #add the tip names in the dataframe
  missingData <- cbind(species = row.names(missingData), missingData)
  
  #convert missingData as character
  missingData[,2] <- as.character(missingData[,2])
  missingData[,2][which(is.na(missingData[,2]))] <- "?" #replace NA by "?"because corHMM don't like it
  #Define the rate model
  model <- "ER"
  FitCorHMM <- corHMM::corHMM(phy = tree, data = missingData, model = model, rate.mat = , 
                              rate.cat = 1, get.tip.states = TRUE)
  
  
  #Calculate AIC
  AIC <- FitCorHMM$AIC
  models <- c("SYM", "ARD")
  model <- "ER"
  for (i in 1:length(models)){
    
    FitCorHMMDiffModel <- corHMM::corHMM(phy = tree, data = missingData, model = models[i], 
                                         rate.cat = 1, get.tip.states = TRUE)
    
    #Calculate AIC
    AICDiffModel <- FitCorHMMDiffModel$AIC
    if(AIC > AICDiffModel){
      AIC <- AICDiffModel
      FitCorHMM <- FitCorHMMDiffModel
      model <- models[i]
    }
  }
  
  #Imputation
  MostLikelyState <- apply(FitCorHMM$tip.states, 1, which.max)
  MostLikelyState <- as.data.frame(MostLikelyState)
  colnames(MostLikelyState) <- colName
  
  #Parameters
  parameters <- list(model = model, rate.cat = 1)
  
  #print(MostLikelyState)
  return(list(imputedData = MostLikelyState, parameters = parameters, probabilities = FitCorHMM$tip.states))
}


#' @title Imputation of missing data for discrete traits, mandatory phylogenetic information  
#' @description This function apply the function imputeOneDiscretTrait on several columns at one time.
#' @usage PIDiscreteTrait(trait)
#' @param missingData data.frame of 1 or more factor columns containing NAs
#' @param tree phylo object
#' @return a data.frame of 1 or more factor columns with the NAs replaced by values, parameters used for the 
#' imputation, matrix of probabilities of each states
#' @export
imputePIDiscreteTraits <- function(missingData, Data){
  
  #select the columns with missing values
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))
  #NaNTraits <- missingData[, NaNColIndex, drop = FALSE]
  
  proba <- matrix(NA, nrow = nrow(missingData), ncol = ncol(missingData))
  parameters <- list()
  for(i in 1:length(names(NaNColIndex))){
    imputation <- imputeOneDiscreteTrait(missingData[, names(NaNColIndex)[i], drop = FALSE], Data)
    #print("ok")
    missingData[ ,names(NaNColIndex)[i]] <- imputation$imputedData
    
    if(i == 1){
      proba <- imputation$probabilities
    }
    
    else{
      proba <- cbind(proba, imputation$probabilities)
    }
    parameters <- c(parameters, imputation$parameters)
    #print(missingData[,i])
  }
  proba <- as.data.frame(proba, row.names = rownames(missingData))
  
  return(list(imputedData = missingData, parameters = parameters, probabilities = proba))
}