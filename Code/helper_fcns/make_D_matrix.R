MakeDMatrix <- function(pop.mesh, loc.name, m) {
  # Derive matrix with population density estimates
  #
  # Args:
  #   pop.mesh: Dataframe with mesh id, mesh location (i.e., area it is in),
  #                and mean
  #   loc.name: Column name (string) of dataframe that is the mesh location
  #   m: Number of mesh points
  #
  # Returns:
  #   A list that has entries 
  #       D -- the relative density matrix (number unique locations x m)
  #       mesh.weights -- approx number of people that live at each mesh point 
  D <- matrix(0, nrow=length(unique(pop.mesh[, loc.name])), ncol=m)
  for(cell in 1:nrow(D)) {
    location.indices <- pop.mesh[, loc.name]==cell
    D[cell, pop.mesh$id[location.indices]] <- pop.mesh$pop[location.indices]
  }
  return(list(D=D/rowSums(D),
              mesh.weights=colSums(D)))
}