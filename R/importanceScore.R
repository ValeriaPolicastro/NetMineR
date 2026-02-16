
#' importanceScore
#'
#' @param graph graph
#' @param memb membership object
#'
#' @returns list
#' @export
#'
#' @examples
#' karate <- make_graph("Zachary")
#' membK <- membership(cluster_walktrap(karate))
#' importanceScore(graph=karate,memb=membK)

importanceScore <- function(graph, memb){
  
   # risoluzione problema nomi in igraph
    if (length(V(graph)$name)<1)
    {
      V(graph)$name <- as.character(V(graph))  
    }
      
      
     
  comm_vertices <- split(V(graph), memb)
    output <- list()
    
    for (i in 1: length(comm_vertices)){
      
      print(i)
      if (length(comm_vertices [[i]])>2){
        sg <- induced_subgraph(graph, vids =  comm_vertices [[i]])
        
  
        
        # find max-degree nodes
        degs <- degree(sg)
        max_deg <- max(degs)
        max_nodes_degree <- V(sg)[degs == max_deg]
       
        
        
        closeness <- closeness(sg)
        max_c <- max(closeness)
        max_nodes_closeness <- V(sg)[closeness == max_c]
        
        
        eigen <- eigen_centrality(sg)$vector
        max_e <- max(eigen)
        max_nodes_eigen <- V(sg)[eigen == max_e]
      } else{
        print("less then 2 nodes")
        max_nodes_degree <- comm_vertices [[i]]
        max_nodes_closeness <- comm_vertices [[i]] 
        max_nodes_eigen <- comm_vertices [[i]]
      }
      
      
      
      output[[i]]<- list(max_degree=max_nodes_degree, 
                         max_closeness=max_nodes_closeness, 
                         max_eigen=max_nodes_eigen)
      
    }
    return(output)  
    
    }
    
    

  

