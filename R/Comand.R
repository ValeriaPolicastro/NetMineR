#### Comand Function ###

#' comand
#' 
#' @description comand COMmunities in Attributed Network with Distatis. Network 
#' with both quantitative and qualitative attributes. (Isolated nodes will be removed)
#' @param adj The adjacency matrix, the nodes should be in the same order for 
#' all the matrices.
#' @param attrQual matrix of qualitative attributes (as factors objects), the nodes should be in the
#'  same order for all the matrices. 
#' @param attrQuant matrix of quantitative attributes, the nodes should be in the
#'  same order for all the matrices. 
#' @param methodQuant distance method for quantitative variables (default euclidean 
#' see dist function stats) 
#' @param methodQual distance method for qualitative variables (default chi.square 
#' see distance function analogue)
#' @return list of outputs
#' @export
#' @importFrom PERMANOVA DistBinary
#' @importFrom stats dist
#' @importFrom utils combn
#' @importFrom analogue distance
#'
#' @examples 
#' pm <- cbind(c(.1, .001), c(.001, .1))
#' g <- igraph::sample_sbm(100, pref.matrix = pm, block.sizes = c(50, 50))
#' Adj <- igraph::as_adjacency_matrix(g)
#' Qual <- data.frame(V1=factor(c(rep("M",50),rep("F",50))),
#' V2=factor(c(rep("I",30),rep("II",30),rep("III",40))))
#' Quant <- data.frame(X1=c(rnorm(50, 2, sd = 0.5),rnorm(50, 10, sd = 0.5)), 
#' X2= c(runif(50, 1, 4), runif(50, 3, 6)))
#' comand(adj=Adj, attrQual=Qual, attrQuant=Quant)

comand <- function(adj, attrQual, attrQuant, methodQuant="euclidean", methodQual="chi.square"){
  
  # i nomi delle righe devono essere allo stesso posto per tutti e 3 
  # devono avere la stessa dimensione del numero di nodi, eliminare i nodi
  # isolati e quindi il loro record nella matrici attrQual e attrQuant
  # inserire controlla nomi nodi (anche ordine)
  # provide with more options for weight
  # controllare l'output del compromesso
  
  
  if(length(rownames(adj))>0)
  {
    ind <- utils::combn(seq_len(length(adj)),2)
    for (i in 1:dim(ind)[2])
    {
      if(any(rownames(adj[[ind[1,i]]]) != rownames(adj[[ind[2,i]]])))
      {
        stop("Rownames and colnames should be in the same order in all the 
                    adjacency matrices.")
      }
    }
  }
  
  
  
  # inserire un avvertimento per i nodi isolati oppure eliminarli già all'interno della funzione
  
  #for (i in 1:length(adj)) per estensione multilayer possiamo fare una lista di adj matrix
  isolated <- which(rowSums(as.matrix(adj))==0)

  
  if(length(isolated)>0)
  {
    warning("isolated nodes will be removed") 
    
    # from adj and attributes qualitative and quantitative
      adj <- adj[-c(isolated),-c(isolated)]
    
    
    if(length(attrQual)>0)
    {
    attrQual <- attrQual[-c(isolated),]}
    
    
    if(length(attrQuant)>0)
    {
    attrQuant <- attrQuant[-c(isolated),]}
  }
  
  
  
  ### Network - distance matrix ####
    # Weighted Networks
  if(length(adj[adj!=0 & adj!=1])>0)
    { 
      #Se frequenze: 
      netD <- analogue::distance(adj, method = "euclidean")
      
    }else{
      # Unweighted Networks  
      netD <- PERMANOVA::DistBinary(adj, coefficient ="Jaccard", transformation = "1-S")$D
    }
  
  
  ### Qualitative attributes - distance matrix ####
  if(length(attrQual)>0)
  {
    attr.qualD <- as.matrix(analogue::distance(attrQual, method = methodQual))
  }
  
  ### Quantitative attributes - distance matrix ####
  if(length(attrQuant)>0)
  { 
    attr.quantD <- as.matrix(stats::dist(attrQuant, method = methodQuant))
  }
  
 
  
  array.dist <- cbind(netD,attr.qualD,attr.quantD)
  name <- c("Net","Qual","Quant")
  
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
  compromise <- Distatis$res4Splus$Splus
  output$compromise <- compromise
  
  ### Coordinate dei nodi nelle prime 3 dimensioni
  output$F <- Distatis$res4Splus$F
  
  ### PartialF serve per plottare Net qual e quant nello spazio dei nodi per ogni nodo
  output$PartialF <- Distatis$res4Splus$PartialF
  
  
  
  # CLUSTERING:
  ### Recover the compromise distance matrix from the compromise matrix (cross-product)
  diag_compr <- diag(compromise)
  diag_compr <- as.matrix(diag_compr)
  
  one <- matrix(1,nrow = dim(compromise)[1], ncol = 1)
  
  # Recover the distance matrix from a cross-product matrix (vedi equazione 4 paper Distatis 2005)
  # D(I × I) = diag(S) %*% t(1) + 1 %*% t(diag(S)) - 2S
  compr_dist <- diag_compr%*%t(one) + one%*%t(diag_compr) - 2*compromise
  
  res <- 
    NbClust::NbClust(diss=as.dist(compr_dist), distance = NULL,min.nc=2, max.nc=20, 
                 method = "ward.D2", index = c("silhouette"))
  
  output$community <- res$Best.partition
  
  return(output)
  
}


