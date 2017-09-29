DeterminePointsInPolys <- function(df.pts, df.areas) {
  # Determines which points are in which polygon (i.e., area)
  # 
  # Args:
  #     df.pts: Dataframe with long and lat of the points 
  #             we would like to "categorize"
  #     df.areas: Datarame that has been fortified from spatial object,
  #               which will tell us the areas
  #
  # Return
  #     vector with entries corresponding to df.pts that has id of the area 
  #       the pt belongs to
  area.loc <- rep(0, nrow(df.pts))
  for (id in unique(df.areas$id)) {
    poly.area <- df.areas[df.areas$id==id, ]
    for (piece in unique(poly.area$piece)){
      poly.subarea <- poly.area[poly.area$piece==piece, ]
      if (any(point.in.polygon(df.pts[, 1], df.pts[, 2],
                               poly.subarea[, 1], poly.subarea[, 2]) > 1)) {
        stop("Add some noise to df.pts-- some points are on vertices or boundaries")
      }
      area.loc <- ifelse(
        point.in.polygon(df.pts[, 1], df.pts[, 2], 
                         poly.subarea[, 1], poly.subarea[, 2]) == 1,
                         id, area.loc)
    }
  }
  return(area.loc)
}
