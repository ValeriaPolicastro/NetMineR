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
#' @importFrom igraph rewire
#'
#' @return list of network and attributes perturbated
#' @export
#'
#' @examples
#' pm <- cbind(c(.1, .001), c(.001, .1))
#' g <- sample_sbm(100, pref.matrix = pm, block.sizes = c(50, 50))
#' Qual <- c(rep("M",50),rep("F",50))
#' Quant <- c(rnorm(50, 2, sd = 0.5), rnorm(50, 10, sd = 0.5))
#' 
#' rewireNetAtt(Network=g, attrQual=Qual, attrQuant=Quant, Prewire=0.5 )
#' 
rewireNetAtt <- function(Network, attrQual, attrQuant, Prewire, seed=1)
{
  
  set.seed(seed)

      P <- ( Prewire * length(attrQual) )/2
      
      # Qualitative
      N1 <- sample(1:length(attrQual), ceiling(P))
      N2 <- sample(1:length(attrQual), ceiling(P))
   print(N1)
   print(N2)
      
      N1.1 <- attrQual[N1] 
      N2.1 <- attrQual[N2] 
      
      attrQual[N1] <- N2.1
      attrQual[N2] <- N1.1
      
      # Quantitative
      N1.quant <- sample(1:length(attrQuant), ceiling(P))
      N2.quant <- sample(1:length(attrQuant), ceiling(P))
     print(N1.quant)
     print(N2.quant)
      
      N1.1.quant <- attrQuant[N1.quant] 
      N2.1.quant <- attrQuant[N2.quant]  
      
      attrQuant[N1.quant] <- N2.1.quant
      attrQuant[N2.quant] <- N1.1.quant
      

    
    PNet <- ( Prewire * gsize(Network))/2
    # Network
    Network <- igraph::rewire(Network, with = keeping_degseq(ceiling(PNet),
                                                                   loop=FALSE))
    
  output <- list("attrQual" = attrQual,
              "attrQuant" = attrQuant,
              "Network" = Network)
return(output)
}
