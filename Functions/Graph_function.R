Connectance = function (Con) {
  adj = st_is_within_distance(Con$geometry, dist = 100, sparse = FALSE)
  #creates an adjacency matrix, where a cell is TRUE if the distance between
  #each other is less than 100 m
  #diag(adj) = FALSE #a geometry shouldn't consider adjacency to itself
  N <- nrow(adj)
  upper_half <- adj[upper.tri(adj, diag = FALSE)]
  L <- sum(upper_half == TRUE)

  Gamma = L/((N*(N-1))/2)
  return(Gamma)
  
}
