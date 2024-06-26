Graph = function (Con) {
  adj = st_is_within_distance(Con$geometry, dist = 100, sparse = FALSE)
  diag(adj) = FALSE
  if (is.list(apply(adj, 1, which))) {
    from = rep(c(1:nrow(Con)), apply(adj, 1, sum))
    to = do.call("c", apply(adj, 1, which))
    edges = data.frame(from = from, to = to)
    
    Centroids = st_centroid(Con[,"geometry"])
    
    net = sfnetwork(nodes = Centroids, edges = edges, directed = FALSE,  edges_as_lines = TRUE) %>% activate("edges") %>% filter(!edge_is_multiple()) %>% as.undirected()
    
    return(net)
  } else {
    net== list()
  }
}