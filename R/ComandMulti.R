#### Comand Multiplex Function ###

#' comandMulti
#' 
#' @description comand COMmunities in Attributed Network with Distatis. Network 
#' with both quantitative and qualitative attributes. (Isolated nodes will be removed)
#' @param adjList The adjacency matrix, the nodes should be in the same order for 
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
#' AdjL <- list(igraph::as_adjacency_matrix(g),igraph::as_adjacency_matrix(g),igraph::as_adjacency_matrix(g))
#' Qual <- data.frame(V1=factor(c(rep("M",50),rep("F",50))),
#' V2=factor(c(rep("I",30),rep("II",30),rep("III",40))))
#' Quant <- data.frame(X1=c(rnorm(50, 2, sd = 0.5),rnorm(50, 10, sd = 0.5)), 
#' X2= c(runif(50, 1, 4), runif(50, 3, 6)))
#' comandMulti(adj=AdjL, attrQual=Qual, attrQuant=Quant)

comandMulti <- function(adjList, attrQual, attrQuant, methodQuant="euclidean", methodQual="chi.square"){
  
  # i nomi delle righe devono essere allo stesso posto per tutti e 3 
  # devono avere la stessa dimensione del numero di nodi, eliminare i nodi
  # isolati e quindi il loro record nella matrici attrQual e attrQuant
  # inserire controlla nomi nodi (anche ordine)
  # provide with more options for weight
  # controllare l'output del compromesso
  
  warning("Rownames (nodes) must be identical and in the same order within each matrices.") 
 
#   reference_names <- NULL
#     
#   for (l in seq_along(adjList)) {
#       
#     if(length(rownames(adjList[[l]]))>0)
#       {
#       # Extract rownames and colnames
#       rn <- rownames(adjList[[l]])
#       cn <- colnames(adjList[[l]])
#       
#       # Check that rownames and colnames are the same within each matrix
#       if (!identical(rn, cn)) {
#         stop("Rownames and colnames must be identical and in the same order within each adjacency matrix.")
#       }
#       
#       # Store the node order of the first matrix
#       if (l == 1) {
#         reference_names <- rn
#       }
#       
#       # Check that all matrices have the same node order
#       if (!identical(rn, reference_names)) {
#         stop("Rownames and colnames should be in the same order in all adjacency matrices.")
#       }
#     }
#   }
#   
#   
#   
# if(length(rownames(Quant))>0)
#   {
#   # Extract rownames 
#   rn <- rownames(Quant)
#   if (!identical(rn, reference_names)) {
#     stop("Rownames and colnames should be in the same order in all adjacency matrices.")
#   }
#   }
# 
#   
#   if(length(rownames(Qual))>0)
#   {
#     rn <- rownames(Qual)
#     if (!identical(rn, reference_names)) {
#       stop("Rownames and colnames should be in the same order in all adjacency matrices.")
#     }
#   }
#   

# Nodi Isolati:
  isolated <- NULL
  for (l in seq_along(adjList)) {
  isolati <- which(rowSums(as.matrix(adjList[[l]]))==0)
  isolated <- c(isolated, isolati)
  }
  
  
  Uisolated <- unique(isolated)
  if(length(Uisolated)>0)
  {
    warning("Removal of isolated nodes:") 
    print(Uisolated)
    for (l in seq_along(adjList)) {
      adjList[[l]] <- adjList[[l]][-c(Uisolated),-c(Uisolated)]
    }
    
    
    if(length(attrQual)>0)
    {
      attrQual <- attrQual[-c(Uisolated),]}
    
    
    if(length(attrQuant)>0)
    {
      attrQuant <- attrQuant[-c(Uisolated),]}
  }
  
  
  
  ### Network - distance matrix ####
  netD <- list()
  for (l in seq_along(adjList)) {
    
    # Weighted Networks
  if(length(adjList[[l]][adjList[[l]]!=0 & adjList[[l]]!=1])>0)
  { 
  
    #Se frequenze: 
    netD[[l]] <- analogue::distance(adjList[[l]], method = "euclidean")
    }
  }else{
  
    # Unweighted Networks  
    netD[[l]] <- PERMANOVA::DistBinary(adjList[[l]], coefficient ="Jaccard", transformation = "1-S")$D
  
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
  
  ###### CONTINUA DA QUI !!!!!!###########
  # Per multilayer network (ci vogliono minimo 3 matrici) 
  # ### Build array
  # if(length(attrQual)>0 & length(attrQuant)>0)
  # {
  #   ### Array 3D  ####
  #   array.dist <- cbind(do.call(cbind, netD),attr.qualD,attr.quantD)
  #   name <- c(name,"Qual","Quant")
  # } else if(length(attrQual)>0 ){
  #   ### Array 3D  ####
  #   array.dist <- cbind(do.call(cbind, netD),attr.qualD)
  #   name <- c(name,"Qual")
  # } else if(length(attrQuant)>0 ){
  #   ### Array 3D  ####
  #   array.dist <- cbind(do.call(cbind, netD),attr.quantD)
  #   name <- c(name,"Quant")
  # }
  # 
  
  #if (!dim(array.dist)[2]>=3) {
  #  stop("At least 3 matrices") 
  #}

  
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


