#### DistatisNetwork Function ###

#' DistatisNetwork
#' 
#' @description distatis applied on networks with quantitative and qualitative attributes 
#' @param Adj list of adjacency matrices (same node names in the same order, can not analyse isolated nodes)
#' @param attrQual matrix of qualitative attributes 
#' @param attrQuant matrix of quantitative attributes 
#' @param weights default NO for binary network or "freq" frequency weights
#' @return
#' @export
#' @import DistatisR
#' @importFrom PERMANOVA DistBinary
#'
#' @examples DistatisNetwork (Adj, attrQual, attrQuant, weights)

DistatisNetwork <- function(Adj, attrQual, attrQuant, weights="NO"){
  
  # i nomi delle righe devono essere allo stesso posto per tutti e 3 
  # devono avere la stessa dimensione del numero di nodi, eliminare i nodi isolati e quindi il loro record nella matrici attrQual e attrQuant
  # inserire controlla nomi nodi (anche ordine)
  # provide with more options for weight
  # controllare l'output del compromesso
  
  if(length(rownames(Adj[[1]]))>0)
  {
    ind <- combn(seq_len(length(Adj)),2)
    for (i in 1:dim(ind)[2])
    {
      if(any(rownames(Adj[[ind[1,i]]]) != rownames(Adj[[ind[2,i]]])))
      {
        stop("Rownames and colnames should be in the same order in all the 
                    adjacency matrices.")
      }
    }
  }
  
  
  
  # inserire un avvertimento per i nodi isolati oppure eliminarli già all'interno della funzione
  isolated <- NULL
  for (i in 1:length(Adj))
  {
    isol <- which(rowSums(as.array(Adj[[i]]))==0)
    isolated <- c(isolated, isol)
  }
  
  if(length(isolated)>0)
  {
    warning("isolated nodes will be removed")
    for (i in 1:length(Adj))
    {
      Adj[[i]] <- Adj[[i]][-c(isolated),-c(isolated)]
    }
    attrQual <- attrQual[-c(isolated),]
    attrQuant <- attrQuant[-c(isolated),]
  }
  
  
  
  
  name <- NULL
  netD <- list()
  for (i in 1:length(Adj)){
    ### Network - distance matrix ####
    # Weighted Networks
    if(weights=="freq")
    { 
      #Se frequenze: 
      netD[[i]] <- analogue::distance(Adj[[i]], method = "chi.square")
      
    }else if (weights=="NO"){
      # Unweighted Networks  
      netD[[i]] <- PERMANOVA::DistBinary(Adj[[i]], coefficient ="Jaccard", transformation = "1-S")$D
    }
    nameI <- paste0("Net",i)
    name <- c(name,nameI)
  }
  
  ### Qualitative attributes - distance matrix ####
  if(length(attrQual)>0)
  {
    attr.qualD <- as.matrix(analogue::distance(attrQual, method = "chi.square"))
  }
  
  ### Quantitative attributes - distance matrix ####
  if(length(attrQuant)>0)
  { 
    attr.quantD <- as.matrix(dist(attrQuant))
  }
  
  ### Build array
  if(length(attrQual)>0 & length(attrQuant)>0)
  {
    ### Array 3D  ####
    array.dist <- cbind(do.call(cbind, netD),attr.qualD,attr.quantD)
    name <- c(name,"Qual","Quant")
  } else if(length(attrQual)>0 ){
    ### Array 3D  ####
    array.dist <- cbind(do.call(cbind, netD),attr.qualD)
    name <- c(name,"Qual")
  } else if(length(attrQuant)>0 ){
    ### Array 3D  ####
    array.dist <- cbind(do.call(cbind, netD),attr.quantD)
    name <- c(name,"Quant")
  }
  
  if (!dim(array.dist)[2]>=3) {
    stop("At least 3 matrices") 
  }
  
  Array <- array(array.dist,
                 dim=c(dim(attrQual)[1], dim(attrQual)[1], length(name)))
  
  if (length(row.names(attrQual)) > 0){
    dimnames(Array)[[1]] <- dimnames(Array)[[2]] <- row.names(attrQual)
  } else {
    dimnames(Array)[[1]] <- dimnames(Array)[[2]] <- seq(1:dim(attrQual)[1])
  }
  
  dimnames(Array)[[3]] <- name
  
  output <- NULL
  ## Apply DISTATIS 
  Distatis <- DistatisR::distatis(Array)
  
  ### Importance of each variable 
  output$alpha <- Distatis$res4Cmat$alpha
  output$ExpVar <- Distatis$res4Cmat$tau
  output$eigValues <- Distatis$res4Cmat$eigValues
  
  ### Coordiate delle variabili
  output$G <- Distatis$res4Cmat$G
  
  
  ### Compromise matrix (it is a cross-product)
  output$compromise <- Distatis$res4Splus$Splus
  
  ### Coordinate dei nodi nelle prime 3 dimensioni
  output$F <- Distatis$res4Splus$F
  
  ### PartialF serve per plottare Net qual e quant nello spazio dei nodi per ogni nodo
  output$PartialF <- Distatis$res4Splus$PartialF
  
  
  
  
  return(output)
  
}


