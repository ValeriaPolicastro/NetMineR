#' Example Data
#'
#' 3 data types: Adjacency matrix, Quantitative and Qualitative attributes.
#' @format ## `exampleL_data`
#' A list of 3 objects:
#' \describe{
#' \item{Adj}{Adjacency matrix.}
#'   \item{Quant}{Quantitative attributes.}
#'   \item{Qual}{Qualitative attributes}
#' }

"exampleL_data"


#pm <- cbind(c(.1, .001), c(.001, .1))
# g <- igraph::sample_sbm(100, pref.matrix = pm, block.sizes = c(50, 50))
#  Adj <- igraph::as_adjacency_matrix(g)
#  Qual <- data.frame(V1=factor(c(rep("M",50),rep("F",50))),
#                        V2=factor(c(rep("I",30),rep("II",30),rep("III",40))))
#  Quant <- data.frame(X1=c(rnorm(50, 2, sd = 0.5),rnorm(50, 10, sd = 0.5)), 
#                        X2= c(runif(50, 1, 4), runif(50, 3, 6)))

# exampleL_data <- list(Adj, Qual, Quant)