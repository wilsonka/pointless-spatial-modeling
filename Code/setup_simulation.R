###################################
##### Read in simulation data #####
###################################

# Load adm0, adm0.5, adm1, 
#   kenya.data (location and sample size of each cluster), and 
#   pop.dat.kenya (population data for Kenya)
load("Data/Kenya/kenya_data.RData")

kenya.df1 <- fortify(adm0)
kenya.df1$id <- as.numeric(kenya.df1$id)

kenya.df47 <- fortify(adm1)
kenya.df47$id <- as.numeric(kenya.df47$id)

kenya.df8 <- fortify(adm0.5)
kenya.df8$id <- as.numeric(kenya.df8$id) + 1

###################################
##### Set-up variables for Sim #### 
###################################
# - create mesh for first 5 scenarios: "mesh.true"
# - create D matrix used for scenario 2 and 3: "D"
#     - a 47 x m matrix
# - create D matrix used for scenario 4 and 5: "D8"
#     - see note above
# - add columns to kenya.data indicating which polygons the EAs reside in

# Create internal mesh points
set.seed(10)
points.kenya <- spsample(adm0, 2500, "regular")
boundary.kenya <- inla.nonconvex.hull(points.kenya@coords, 
                                      convex=0.1, 
                                      resolution=c(120,120))

# Determine population near these mesh points
pop.at.mesh <- FindPopAtMesh(data.frame(points.kenya@coords), pop.dat.kenya)

# If not ordered, the mesh is slightly different which affects simulation
pop.at.mesh <- pop.at.mesh[order(pop.at.mesh$mx, pop.at.mesh$my), ]

# Make the mesh
mesh.true <- inla.mesh.2d(pop.at.mesh[, 1:2], offset=c(1, 2), 
                          max.edge = 1*c(1, 5))

### adm1
# Determine which polys the mesh points live in
meshlocs <- DeterminePointsInPolys(mesh.true$loc, kenya.df47)

tapply(rep(1, length(meshlocs)), meshlocs, sum)

mesh.df <- data.frame(mlong=mesh.true$loc[, 1],
                      mlat=mesh.true$loc[, 2],
                      loc=meshlocs)
mesh.df$id <- 1:nrow(mesh.df)

names(pop.at.mesh)[1:2] <- c("mlong","mlat")

### adm0.5
# Determine which polys the mesh points live in

# There is slight misalignment of the boundaries, which should completely 
#  overlap. Thus, we manually enter the provinces based on which of the 47
#  disctricts the mesh point is in.
meshlocs8 <- DeterminePointsInPolys(mesh.true$loc, kenya.df8)
mesh.df$loc8X <- meshlocs8

mapping47to8 <- as.vector(apply(table(mesh.df[, c("loc", "loc8X")]),
                                1, which.max) - 1)
mesh.df$loc8 <- mapping47to8[mesh.df$loc + 1]

### Population at mesh including polygon location
pop.at.mesh <- plyr::join(mesh.df, pop.at.mesh, type="right")

### Create D matrix
# adm1
D.tmp <- MakeDMatrix(pop.at.mesh, "loc", mesh.true$n)
mesh.df$weight.unscaled <- D.tmp$mesh.weights
D <- D.tmp$D

mesh.df$weight.scaled <- colSums(D)
nmesh.area <- 1 / tapply(rep(1, nrow(mesh.df)), mesh.df$loc, sum)
mesh.df$weight.comp <- nmesh.area[mesh.df$loc + 1] # comparison weight


# adm0.5
D.tmp8 <- MakeDMatrix(pop.at.mesh, "loc8", mesh.true$n)
mesh.df$weight.unscaled8 <- D.tmp8$mesh.weights
D8 <- D.tmp8$D
mesh.df$weight.scaled8 <- colSums(D8)
nmesh.area8 <- 1 / tapply(rep(1, nrow(mesh.df)), mesh.df$loc8, sum)
mesh.df$weight.comp8 <- nmesh.area8[mesh.df$loc8 + 1] # comparison weight
