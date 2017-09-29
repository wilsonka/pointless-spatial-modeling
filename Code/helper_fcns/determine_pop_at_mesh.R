CalcEucDist <- function(pt1, comparisonpts) {
  sqrt((as.numeric(pt1[2]) - comparisonpts[,2])^2 + 
    (as.numeric(pt1[1]) - comparisonpts[,1])^2)
}

FindPopAtMesh <- function(mesh.pts, pop, n.digits.round=2) {
  # Determine population at mesh points
  #
  # Args:
  #   mesh.pts: Dataframe with 2 columns (x: long, y: lat)
  #   pop: Dataframe with 3 columns (x: long, y: lat, pop)
  #   n.digits.round: number of decimal places to round to
  #
  # Return:
  #   Dataframe with long and lat from mesh and mean population at that point
  
  # Obtain rounded mesh locations
  points.rounded <- data.frame(mx=mesh.pts[, 1], 
                               my=mesh.pts[, 2], 
                               rx=round(mesh.pts[, 1], n.digits.round), 
                               ry=round(mesh.pts[, 2], n.digits.round))
  
  # Obtain rounded population locations
  pop$rpx <- round(pop[, 1], n.digits.round)
  pop$rpy <- round(pop[, 2], n.digits.round)
  
  x.unique <- unique(pop[,1])
  y.unique <- unique(pop[,2])
  
  deltax <- max(x.unique[-1] - x.unique[-length(x.unique)])
  deltay <- max(y.unique[-1] - y.unique[-length(y.unique)])
  max.dist <- sqrt((deltax/2)^2 + (deltay/2)^2)
  
  pop.mesh <- points.rounded
  pop.mesh$pop <- 0
  for(i in 1:nrow(mesh.pts)) {
    rx <- pop.mesh[i,]$rx
    ry <- pop.mesh[i,]$ry
    tmp <- pop[pop$rpx %in% c(rx-1*10^(-n.digits.round), rx,
                              rx+1*10^(-n.digits.round)) & 
                 pop$rpy %in% c(ry-1*10^(-n.digits.round), ry,
                                ry+1*10^(-n.digits.round)), ]
    if(nrow(tmp) == 0) {
      pop.mesh$pop[i] <- 0
    } else {
      distances <- CalcEucDist(pop.mesh[i, 1:2], 
                               tmp[, 1:2])
      if(all(distances > max.dist)) {
        pop.mesh$pop[i] <- 0
      } else {
        pop.mesh$pop[i] <- tmp$pop[which.min(distances)]
      }
      
    }
  }
  return(pop.mesh)
}