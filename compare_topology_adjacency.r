#'\code{createAdjacencyMatrix}
#'
#'Creation of an adjacency matrix based on the inicial network feeded to CARNIVAL
#'@param network Data frame of a subset of interactions from scafoldNET.
#'Same values as scafoldNET, but it can aditionaly contain a column with
#'weights
#'@param scafoldNET Data frame of values source, interaction, target
#'@param weighted Boolean (False by default). Whether the network contains 
#'weights for the interactions or not
#'
#'@return Adjency matrix for the network. 
#'The sorce nodes go in rows, and the target in columns.
#'The values are positive or negative to indicate the sign. 
#'If weights availables, these values are included. 
#'If not, only presence absence (+/-)1/0 is returned. 
#'
#'@export
#'
#'Rosa Hernansaiz Ballesteros, 2020

createAdjacencyMatrix <- function(network=NULL, scafoldNET=NULL, weighted=F){
  colnames(scafoldNET) = c('source', 'interaction', 'target')
  scafoldNET$source = gsub("-", "_", scafoldNET$source, ignore.case = FALSE, perl = FALSE,
       fixed = FALSE, useBytes = FALSE)
  scafoldNET$target = gsub("-", "_", scafoldNET$target, ignore.case = FALSE, perl = FALSE,
                           fixed = FALSE, useBytes = FALSE)
  
  colnames(network)[1:3] = c('source', 'interaction', 'target')
  if (weighted){ colnames(network)[4] = 'weight' }
  
  #Create a matrix from interaction file
  adjacency = matrix(data = 0, 
                      nrow = length(unique(scafoldNET$source)), 
                      ncol = length(unique(scafoldNET$target)),
                      byrow = FALSE,
                      dimnames = list(sort(unique(scafoldNET$source)),
                                      sort(unique(scafoldNET$target))))
  
  # Find duplicated interactions with different sign
  #duplicated(network[,c("source","target")])
  #fill matrix
  for (n in 1:nrow(network)){
    adjacency[network$source[n], network$target[n]] = as.numeric(network$interaction[n])*as.numeric(network$weight[n])
  }
  
  return(adjacency)
  
}


#'\code{compareAdjacencies}
#'
#'Comparison of adjacency matrices
#'@param adjMAT1 adjacency matrix with or without weights
#'@param adjMAT2 adjacency matrix with or without weights
#'@param weighted Boolean (False by default). Whether taking into account
#'the match of weights or only the existence of the interaction.
#'
#'@return list of three elements:
#' - adjacency matrix with unique interactions for matrix 1
#' - adjacency matrix with unique interactions for matrix 2
#' - adjacency matrix with shared interactions
#'
#'@export
#'
#'Rosa Hernansaiz Ballesteros, 2020

compareAdjacencies <- function(adjMAT1=NULL, adjMAT2=NULL, weighted=F){

  # check that both matrices have same dimensions, and stop otherwhise
  stopifnot(ncol(adjMAT1)==ncol(adjMAT2), 
            all(colnames(adjMAT1)%in%colnames(adjMAT2)))
  
  stopifnot(nrow(adjMAT1)==nrow(adjMAT2), 
            all(row.names(adjMAT1)%in%row.names(adjMAT2)))
  
  # remove columns and rows that contain 0 for both matrices
  intrm = list()
  for (i in 1:2){
    aux1 = apply(adjMAT1, i, function(x) all(x==0))
    aux1 = aux1[aux1]
    aux2 = apply(adjMAT2, i, function(x) all(x==0))
    aux2 = aux2[aux2]
    
    intrm[[i]] <- intersect(names(aux1), 
                      names(aux2))
  }
  adjMAT1 = adjMAT1[!(row.names(adjMAT1)%in%intrm[[1]]), !(colnames(adjMAT1)%in%intrm[[2]])]
  adjMAT2 = adjMAT2[!(row.names(adjMAT2)%in%intrm[[1]]), !(colnames(adjMAT2)%in%intrm[[2]])]
  
  # Compare matrices
  sources = row.names(adjMAT1)
  targets = colnames(adjMAT1)
  
  adjMtxList = list()
  for (nm in c('sharedMTX', 'uMtx1', 'uMtx2')){
    adjMtxList[[nm]] = matrix(data = 0, 
                              nrow = length(sources), 
                              ncol = length(targets),
                              byrow = FALSE,
                              dimnames = list(sources,
                                              targets))
  }
  
  if(!weighted){
    adjMAT1[adjMAT1 < 0] <- -1
    adjMAT1[adjMAT1 > 0] <- 1
    adjMAT2[adjMAT2 < 0] <- -1
    adjMAT2[adjMAT2 > 0] <- 1
  }

  for(source in sources){
    for(target in targets){
      if( adjMAT1[source, target] == adjMAT2[source, target] ){

        adjMtxList[['sharedMTX']][source, target] = adjMAT1[source, target]
        
      }else{
        
        adjMtxList[['uMtx1']][source, target] = adjMAT1[source, target]
        adjMtxList[['uMtx2']][source, target] = adjMAT2[source, target]
        
      }
    }
    
  }
  
  # Remove columns and rows that contain 0 for both matrices
  
  adjMtxList = lapply(adjMtxList, function(mtx){
    rmr = apply(mtx, 1, function(x) all(x==0))
    rmc = apply(mtx, 2, function(x) all(x==0))
    mtx = mtx[!(row.names(mtx)%in%rmr),
              !(colnames(mtx)%in%rmc)]
    return(mtx)
  })

  return(adjMtxList)
  
}



