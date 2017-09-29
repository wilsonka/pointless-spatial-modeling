###################################
###### Read in scotland data ######
###################################

load("Data/Scotland/scot_data.RData")
scot.coords <- coordinates(scot.BNG)/1000  # change from m --> km
scot.coords <- scot.coords[order(as.numeric(rownames(scot.coords))), ]

scot.df <- fortify(scot.BNG,region="ID")  
# long is Eastings
# lat is Northings
scot.df$id <- as.numeric(scot.df$id)
scot.df$long <- scot.df$long/1000  # change to km
scot.df$lat <- scot.df$lat/1000  # change to km
scot.df$Obs <- scot.dat$Observed[scot.df$id]
scot.df$Exp <- scot.dat$Expected[scot.df$id]
scot.df$SIR <- scot.df$Obs/scot.df$Exp

rownames(scot.coords) <- as.numeric(rownames(scot.coords))

###################################
## Set-up variables for analysis ##
###################################
# - create mesh: "mesh.true"
# - create D matrix: "D"
#     - note: this is a 56 x m matrix
# - TMB objects

# Create Mesh points
set.seed(1)
points.scotland <- spsample(scot.BNG, 2000, "regular")

mesh.pts <- data.frame(points.scotland@coords)/1000
pop.at.mesh <- FindPopAtMesh(mesh.pts, pop.dat.scot, n.digits.round=0)

# Create mesh
mesh.true <- inla.mesh.2d(pop.at.mesh[, 1:2], offset=c(1, 2)*100, 
                          max.edge = 100*c(1, 5))
#mesh.true$n  # 2417 mesh points

# Mesh points in which areas
meshlocs <- DeterminePointsInPolys(mesh.true$loc, scot.df)
mesh.df <- data.frame(mx=mesh.true$loc[, 1], my=mesh.true$loc[, 2],
                      loc=meshlocs)
mesh.df$id <- 1:nrow(mesh.df)

# Create df that has population at mesh including polygon location
pop.at.mesh.new <- plyr::join(mesh.df, pop.at.mesh, type="right") 
Dtmp <- MakeDMatrix(pop.at.mesh.new, "loc", mesh.true$n)
mesh.df$weight.unscaled <- Dtmp$mesh.weights
D <- Dtmp$D

mesh.df$weight.scaled <- colSums(D)
nmesh.area <- 1 / tapply(rep(1, nrow(mesh.df)), mesh.df$loc, sum)
mesh.df$weight.comp <- nmesh.area[mesh.df$loc + 1]  # comparison weight

# Prep for TMB & HMC
spde <- inla.spde2.matern(mesh.true)

Data.tmb <- list("y_i"=scot.dat$Observed,
                 "e_i"=scot.dat$Expected,
                 "M0"=spde$param.inla$M0,
                 "M1"=spde$param.inla$M1,
                 "M2"=spde$param.inla$M2,
                 "A"=inla.as.sparse(D) )
Params.tmb <- list("alpha"=0, 
                   "theta"=spde$param.inla$theta.initial, 
                   "S_j"=rep(0, spde$n.spde))

