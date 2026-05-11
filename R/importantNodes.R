#' importantNodes
#'
#' @param graph graph
#' @param memb membership object
#'
#' @returns list and plot
#' @export
#'
#' @examples
#' karate <- make_graph("Zachary")
#' membK <- membership(cluster_walktrap(karate))
#' importantNodes(graph=karate,memb=membK)

importantNodes <- function(graph, memb){
  
  # risoluzione problema nomi in igraph
  if (length(V(graph)$name)<1)
  {
    V(graph)$name <- as.character(V(graph))  
  }
  
  
  
  comm_vertices <- split(V(graph), memb)
  outputComm <- NULL
  df <- NULL
  for (i in 1:length(comm_vertices)){
    
    print(i)
    if (length(comm_vertices[[i]])>2){
      sg <- induced_subgraph(graph, vids =  comm_vertices [[i]])
      
      
      
      # find max-degree nodes
      degs <- degree(sg)
      max_deg <- max(degs)
      max_nodes_degree <- names(V(sg)[degs == max_deg])
      d <- cbind(max_deg, max_nodes_degree, i)
      rownames(d) <- rep("degree",  dim(d)[1])
      
      closeness <- closeness(sg)
      max_c <- max(closeness)
      max_nodes_closeness <- names(V(sg)[closeness == max_c])
      c <- cbind(max_c, max_nodes_closeness, i)
      rownames(c) <- rep("closeness",  dim(c)[1])
      
      #Seed random si modifica ad ogni esecuzione:
      set.seed(100)
       eigen <- eigen_centrality(sg)$vector
       max_e <- max(eigen)
      max_nodes_eigen <- names(V(sg)[eigen == max_e])
      e <- cbind(max_e, max_nodes_eigen, i)
       rownames(e) <-rep( "eigen_centrality",  dim(e)[1])
      
      
      be <- betweenness(sg)
      max_be <- max(be)
      max_nodes_be <- names(V(sg)[be == max_be])
      b <- cbind(max_be, max_nodes_be, i)
      rownames(b) <-rep( "betweenness",  dim(b)[1])
      
      
    } else{
      print("less then 2 nodes")
      
    }
    

    outputC <- rbind(d,c,e,b)
    

    dat <- data.frame(table(outputC[,2]), rep(paste0("Community",i),length(table(outputC[,2]))))
    colnames(dat) <- c("Nodes","Value","Community")
    df <- rbind(df, dat)
    
    outputComm <- rbind(outputComm,outputC) 
    
  }
  
  
   plot <- ggplot2::ggplot(df, aes(x = reorder(Nodes, Value), y = Value, fill = Value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ Community, scales = "free_y") +
    theme_minimal() +
    labs(title = "Communities",
         x = "Nodes",
         y = "Value") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme(legend.position = "none")
 
   print(plot)
 
 
 
  colnames(outputComm) <- c("measure", "node", "community")
  
  
  return(output=list(data_summary=df, 
                     nodes_measure=outputComm))  
  
}
