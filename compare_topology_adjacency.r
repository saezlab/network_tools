#'\code{createAdjencyMatrix}
#'
#'Creation of an adjency matrix based on the inicial network feed to CARNIVAL
#'@param scafoldNET Data frame of values source, interaction, target
#'@param network Data frame of a subset of interactions from scafoldNET.
#'Same values as scafoldNET, but it can aditionaly contain a column with
#'weights
#'@param weighted Boolean (False by default). Whether the network contains 
#'weights for the interactions or not
#'
#'@return Adjency matrix for the network. 
#'The sorce nodes go in rows, and the target in columns.
#'The values are positive or negative to indicate the sign. 
#'If weights availables, these values are included. 
#'If not, only presence absence (+/-)1/0 is returned. 
#'@import doParallel
#'@import CARNIVAL
#'@import readr
#'@import readxl
#'@import lpSolve
#'
#'@export
#'
#'Rosa HErnansaiz Ballesteros, 2020

createAdjencyMatrix <- function(scafoldNET=NULL, network=NULL, weighted=F){
  colnames(scafoldNET) = c('source', 'interaction', 'target')
  scafoldNET$source = gsub("-", "_", scafoldNET$source, ignore.case = FALSE, perl = FALSE,
       fixed = FALSE, useBytes = FALSE)
  scafoldNET$target = gsub("-", "_", scafoldNET$target, ignore.case = FALSE, perl = FALSE,
                           fixed = FALSE, useBytes = FALSE)
  
  colnames(network)[1:3] = c('source', 'interaction', 'target')
  if (weighted){ colnames(network)[4] = 'weight' }
  
  #Create a matrix from interaction file
  scafoldMTX = matrix(data = 0, 
                      nrow = length(unique(scafoldNET$source)), 
                      ncol = length(unique(scafoldNET$target)),
                      byrow = FALSE,
                      dimnames = list(sort(unique(scafoldNET$source)),
                                      sort(unique(scafoldNET$target))))
  
  # Find duplicated interactions with different sign
  #duplicated(network[,c("source","target")])
  #fill matrix
  for (n in 1:nrow(network)){
    scafoldMTX[network$source[n], network$target[n]] = as.numeric(network$interaction[n])*as.numeric(network$weight[n])
  }
  
  return(scafoldMTX)
  
}



