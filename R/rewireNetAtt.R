## Function for rewiring ##

#' rewireNetAtt
#' @description function to rewire network with attributes quantitative and 
#' qualitative
#' 
#' @param Network network
#' @param attrQual attributes qualitative 
#' @param attrQuant attributes quantitative
#' @param Prewire percentage of perturbation
#' @param seed to set a seed (default seed 1)
#'
#' @return list of network and attributes perturbated
#' @export
#'
#' @examples
#' 
#' 
rewireNetAtt <- function(Network, attrQual, attrQuant, Prewire, seed=1)
{
  
  attrQual.rew <- attrQual  
  attrQuant.rew <- attrQuant
  BMNet.rew <- Network
  
  set.seed(seed)
  for (i in 1:length(attrQual)){

      
      P <- ( Prewire * length(attrQual[[i]][,1]) )/2
      
      # Qualitative
      N1 <- sample(1:length(attrQual[[i]][,1]), ceiling(P))
      N2 <- sample(1:length(attrQual[[i]][,1]), ceiling(P))
      # print(N1)
     # print(N2)
      
      N1.1 <- attrQual[[i]][N1,] 
      N2.1 <- attrQual[[i]][N2,] 
      
      attrQual.rew[[i]][N1,] <- N2.1
      attrQual.rew[[i]][N2,] <- N1.1
      
      # Quantitative
      N1.quant <- sample(1:length(attrQuant[[i]][,1]), ceiling(P))
      N2.quant <- sample(1:length(attrQuant[[i]][,1]), ceiling(P))
      print(N1.quant)
      print(N2.quant)
      
      N1.1.quant <- attrQuant[[i]][N1.quant,] 
      N2.1.quant <- attrQuant[[i]][N2.quant,]  
      
      attrQuant.rew[[i]][N1.quant,] <- N2.1.quant
      attrQuant.rew[[i]][N2.quant,] <- N1.1.quant
      

    
    PNet <- ( Prewire * gsize(BMNet.rew[[i]]))/2
    # Network
    BMNet.rew[[i]] <- rewire(BMNet.rew[[i]], with = keeping_degseq(ceiling(PNet),
                                                                   loop=FALSE))
    
  }
  output <- list("attrQual" = attrQual.rew,
              "attrQuant" = attrQuant.rew,
              "BMNet" = BMNet.rew)
return(output)
}
