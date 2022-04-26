# imputationMethods

This package allows the user to impute missing values in dataset containing countinuous and/or discrete data.

To install the package from Github using devtools.

```
install.packages("devtools")
install.packages("visdat") 
library(devtools)

devtools::install_github("Matgend/imputationMethods")
```
```
#Load missingData
NaData <- get(load(system.file("extdata", "missingData.RData", package = "imputationMethods")))
missingData <- NaData$DataNaN$MCAR$`MCAR/AllTraits/13/0.33`

#Load phylogenetic tree
trueData <- get(load(system.file("extdata", "trueDataAndTree.RData", package = "imputationMethods")))
tree <- trueData$TreeList$`0`
plot(tree)


#visualize the missing value distribution
vis_miss(missingData)
```
