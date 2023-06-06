require("igraph")

rgumbel = function(n){
  -log(-log(runif(n)))
}


findMST= function(logS){
  G<- graph_from_adjacency_matrix(-logS, mode = "undirected",weighted = TRUE,diag=FALSE)
  G_mst<- mst(graph = G)
  A_T_<- as.matrix(get.adjacency(G_mst))
  A_T_
}


drawT_fast= function(logS){
  n<- nrow(logS)
  gumbelMat<- matrix(0,n,n)
  gumbelMat[lower.tri(gumbelMat,diag = FALSE)]<- rgumbel((n)/2*(n-1))
  gumbelMat<- gumbelMat+ t(gumbelMat)
  A_T_ <- findMST(logS+gumbelMat)
  A_T_
}