coveringEdges2 <- function(combs){

  ed <- list()

  #combs <- lapply(combs, function(x){
  #  sort(c(x,0))
  #})

  for (node in 1:length(combs)){
    print(node)

    containedidx <- lapply(combs, function(x){
      length(x) == sum(x %in% combs[[node]]) && length(x) != 0 && !identical(x, combs[[node]])
    })

    contained <- which(unlist(containedidx) == T)
    #contained <- which(containedidx == T)

    contained <- contained[order(sapply(combs[contained],length),decreasing=T)]

    if (length(contained) > 0){
      un <- c()
      pais <- c()
      for (i in 1:length(contained)){
        un <- sort(unique(c(un, combs[[contained[i]]])))
        #print(un)
        pais <- c(pais, contained[i])
        if (identical(un, combs[[node]])){
          break;
        }
      }

      for (p in pais){
        # find the id in combs
        ed[[length(ed) + 1]] <- c(p, node)
      }
    }
  }

  return(ed)
}



cb <- list()

cb[[1]] <- c(0)
cb[[2]] <- c(0,1)
cb[[3]] <- c(0,2)
cb[[4]] <- c(0,3)
cb[[5]] <- c(0,1,2)
cb[[6]] <- c(0,1,2,3)

combs = cb

edges <- coveringEdges2(cb)

no=unlist(lapply(combs, paste, collapse = "-"))
no[which(no == "")] <- "0";
dd <- do.call(rbind.data.frame, edges)
dd <- data.frame(origem=no[dd[,1]], destino=no[dd[,2]])
colnames(dd) <- c("origem", "destino")

library("igraph")
grafo1 <- graph_from_data_frame(d = dd, vertices = no, directed = F)
#layout <- layout_with_kk(grafo1)
layout <- layout_as_tree(grafo1, root = c(6))
#layout[3,2] <- 0

#layout <- layout_with_sugiyama(grafo1,maxiter = 1000, hgap = 100)

#plot(grafo1, layout=layout_as_tree(grafo1, root = c(1)), vertex.size=30)
plot(grafo1, layout=layout, vertex.size=15)




