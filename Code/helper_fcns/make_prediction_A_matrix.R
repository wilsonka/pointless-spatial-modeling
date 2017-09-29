MakePredictionAMatrix <- function(mesh.true, adm0, dims=c(120, 120)) {
  # Make a smaller version of the A matrix so that predicting the risk surface
  #   is less computationally involved
  #
  # Args:
  #   mesh.true: The mesh
  #   adm0: Spatial polygons dataframe that is for the boundary of Kenya
  #   
  # Return: list with prediction matrix 
  #   and points where predictions will be made
  smallproj <- inla.mesh.projector(mesh.true, dims=dims)
  x <- smallproj$x
  y <- smallproj$y
  all.temp <- (expand.grid(x, y , KEEP.OUT.ATTRS=F))
  names(all.temp) <- c("long", "lat")
  coordinates(all.temp) <- ~ long + lat
  all.temp@proj4string <- adm0@proj4string 
  rel.points <- all.temp[adm0, ]
  rel.points.clip <- as.data.frame(rel.points)
  A.pred <- inla.spde.make.A(mesh.true, loc=as.matrix(rel.points.clip))
  return(list(A.pred=A.pred,
              rel.points.clip=rel.points.clip))
}