## Function for rewiring

#' rewireNet
#' @description
#' function to generate simulation
#' 
#' @param Network network
#' @param Prewire percentuale di perturbazioni
#' @param block logical il TRUE it rewires a block of rows across all the 
#' variables, if FALSE it makes a rewire one variable at the time (default TRUE) 
#' @import igraph
#' @return list
#' @export
#'
#' @examples  pm <- cbind(c(.1, .001), c(.001, .1))
#' g <- igraph::sample_sbm(100, pref.matrix = pm, block.sizes = c(50, 50))
#' Adj <- igraph::as_adjacency_matrix(g)
#' Qual <- data.frame(V1=factor(c(rep("M",50),rep("F",50))),
#' V2=factor(c(rep("I",30),rep("II",30),rep("III",40))))
#' Quant <- data.frame(X1=c(rnorm(50, 2, sd = 0.5),rnorm(50, 10, sd = 0.5)), 
#' X2= c(runif(50, 1, 4), runif(50, 3, 6)))
#' rewireNet(adj=g, attrQual=Qual, attrQuant=Quant)
#' 
#' 
rewireNet <- function(attrQual, attrQuant, Network, Prewire, block=TRUE)
{
  
  attrQual.rew <- attrQual  
  attrQuant.rew <- attrQuant
  BMNet.rew <- Network
  
  if (block== FALSE)
  {
    set.seed(100)
    for (i in 1:length(attrQual)){
      for (j in 1:dim(attrQual[[i]])[2]){
        
        P <- ( Prewire * length(attrQual[[i]][,j]) )/2
        
        # Qualitative
        N1 <- sample(1:length(attrQual[[i]][,j]), ceiling(P))
        N2 <- sample(1:length(attrQual[[i]][,j]), ceiling(P))
        
        #print(N1)
        #print(N2)
        
        N1.1 <- attrQual[[i]][N1,j] 
        N2.1 <- attrQual[[i]][N2,j] 
        
        
        attrQual.rew[[i]][N1,j] <- N2.1
        attrQual.rew[[i]][N2,j] <- N1.1
        
        # Quantitative
        N1.quant <- sample(1:length(attrQuant[[i]][,j]), ceiling(P))
        N2.quant <- sample(1:length(attrQuant[[i]][,j]), ceiling(P))
        #print(N1.quant)
        #print(N2.quant)
        
        N1.1.quant <- attrQuant[[i]][N1.quant,j] 
        N2.1.quant <- attrQuant[[i]][N2.quant,j]  
        
        attrQuant.rew[[i]][N1.quant,j] <- N2.1.quant
        attrQuant.rew[[i]][N2.quant,j] <- N1.1.quant
        
      }
      
      PNet <- ( Prewire * igraph::gsize(BMNet.rew[[i]]))
      # Network
      BMNet.rew[[i]] <- igraph::rewire(BMNet.rew[[i]], with = keeping_degseq(ceiling(PNet),
                                                                     loop=FALSE))
      
    }
  }else{
    print("block method")
    attrQual.rew <- attrQual  
    attrQuant.rew <- attrQuant
    BMNet.rew <- Network
    
    set.seed(100)
    for (i in 1:length(attrQual)){
      
      P <- (Prewire * (dim(attrQual[[i]])[1]))/2
      
      # Qualitative
      N1 <- sample(1:(dim(attrQual[[i]])[1]), ceiling(P))
      N2 <- sample(1:(dim(attrQual[[i]])[1]), ceiling(P))
      
      #print(N1)
      #print(N2)
      
      N1.1 <- attrQual[[i]][N1,] 
      N2.1 <- attrQual[[i]][N2,] 
      
      
      attrQual.rew[[i]][N1,] <- N2.1
      attrQual.rew[[i]][N2,] <- N1.1
      
      # Quantitative
      N1.quant <- sample(1:(dim(attrQual[[i]])[1]), ceiling(P))
      N2.quant <- sample(1:(dim(attrQual[[i]])[1]), ceiling(P))
      #print(N1.quant)
      #print(N2.quant)
      
      N1.1.quant <- attrQuant[[i]][N1.quant,] 
      N2.1.quant <- attrQuant[[i]][N2.quant,]  
      
      attrQuant.rew[[i]][N1.quant,] <- N2.1.quant
      attrQuant.rew[[i]][N2.quant,] <- N1.1.quant
      
      
      
      
      
      PNet <- ( Prewire * igraph::gsize(BMNet.rew[[i]]))
      # Network
      BMNet.rew[[i]] <- igraph::rewire(BMNet.rew[[i]], with = keeping_degseq(ceiling(PNet),
                                                                     loop=FALSE))
    }
    
  }
  
  
  return(list("attrQual" = attrQual.rew,
              "attrQuant" = attrQuant.rew,
              "BMNet" = BMNet.rew))
}
