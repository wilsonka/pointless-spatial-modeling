######################
### General Set-up ###
######################

# Parameters
alpha <- 0
sigma2 <- 0.25
lambda2 <- 1-sigma2
erange <- sqrt(8)/exp(0.5)
theta2 <- log(sqrt(8)/erange)
kappa <- exp(theta2)
tau2 <- 1/(4*pi*kappa^2*lambda2)
theta1 <- log(sqrt(tau2))

# Set-up Field
spde <- inla.spde2.matern(mesh.true)
Q.true <- inla.spde.precision(spde, c(theta1, theta2))
field.true <- as.numeric(inla.qsample(1, Q=Q.true, seed=50L))

# Connecting field to data
proj.survey <- inla.mesh.projector(mesh.true, 
                                   loc=as.matrix(kenya.data[, c("LONGNUM", 
                                                                "LATNUM")]))
field.survey <- inla.mesh.project(proj.survey, field=field.true)
kenya.data$field <- field.survey
kenya.data$mu <- alpha + kenya.data$field

# Make matrix that will be used to obtain predictions of the surface
pred.info <- MakePredictionAMatrix(mesh.true, adm0, dims=c(110, 110))
A.pred <- pred.info$A.pred
stkgrid <- inla.stack(
  data=list(y=NA,
            N=NA),
  A=list(A.pred, 1),
  effects=list(i=1:spde$n.spde, alpha=rep(1, nrow(A.pred))),
  tag='prd.gr')
rel.points.clip <- pred.info$rel.points.clip


######################
#### Simulate Data ###
######################


# Census
# data simulated on 1kmx1km grid
pop.data.kenya.points <- data.frame(rasterToPoints(pop.dat.kenya, spatial = F))
names(pop.data.kenya.points) <- c("long", "lat", "pop")
coordinates(pop.data.kenya.points) <- ~ long + lat
pop.data.kenya.points@proj4string <- adm0@proj4string
pop.pt.insideX <- pop.data.kenya.points[adm0, ]
pop.pt.inside <- as.data.frame(pop.pt.insideX)
coordinates(pop.pt.inside) <- ~ long + lat
pop.pt.inside@proj4string <- adm0@proj4string


### 47 Areas
poplocs47 <- over(pop.pt.inside, adm1)[,1]


# Surveys
set.seed(430)
kenya.data$ybar <- rnorm(400, kenya.data$mu, sd=sqrt(sigma2/kenya.data$N))


# for truth defined at mesh points
proj.to.pop <- inla.mesh.projector(mesh.true, 
                                   loc=as.matrix(pop.pt.inside@coords))

field.at.pop <- inla.mesh.project(proj.to.pop, 
                                  field=field.true)
pop.inside <- pop.pt.inside
pop.inside$mu <- alpha + field.at.pop



pop.inside$locs47 <- poplocs47
popi <- tapply(pop.inside$pop, pop.inside$locs47, sum)
householdsizei <- popi/3.9
mui <- tapply(pop.inside$mu * pop.inside$pop, pop.inside$locs47, sum) / popi

set.seed(999)
census.df47 <- data.frame(id=1:47, ybar=rnorm(47, mui, sd=sqrt(sigma2/householdsizei)),
                          N=householdsizei)

### 8 Areas
poplocs8 <- mapping47to8[poplocs47]
pop.inside$locs8 <- poplocs8
popi8 <- tapply(pop.inside$pop, pop.inside$locs8, sum)
householdsizei8 <- popi8/3.9
mui8 <- tapply(pop.inside$mu * pop.inside$pop, pop.inside$locs8, sum) / popi8

set.seed(9999)
census.df8 <- data.frame(id=1:8, ybar=rnorm(8, mui8, sd=sqrt(sigma2/householdsizei8)),
                         N=householdsizei8)

######################
####### Models #######
######################
f.spde <- y ~ -1 + alpha + f(i,model=spde)

# Option 1: all coordinates
A <- inla.spde.make.A(mesh.true, loc=as.matrix(kenya.data[, c("LONGNUM", 
                                                              "LATNUM")]))
source("Code/simulation_models/run_simulation_points.R")

# Option 2: 47 areas
source("Code/simulation_models/run_simulation_47areas.R")

# Option 3: Survey + 47
D.both <- rbind(A, D)
source("Code/simulation_models/run_simulation_both.R")

# Option 4: 8 Provinces
source("Code/simulation_models/run_simulation_8areas.R")

# Option 5: Survey + 8
D.both8 <- rbind(A, D8)
source("Code/simulation_models/run_simulation_8areas_both.R")



######################
#### Descriptives ####
######################

# Kenya info
mesh.true$n  # m

# Surveys
summary(kenya.data$N)  # N_{ij}
sum(kenya.data$N) # Total number of households surveyed

# Census
summary(census.df47$N)
sum(census.df47$N)

summary(census.df8$N)
sum(census.df8$N)

# Field Info
1 / (4 * pi * exp(theta1)^2 * exp(theta2)^2) # marginal variance
sqrt(8) / exp(theta2)  # practical range
mean(field.true)  # 0.182
var(field.true)  # 0.672


######################
## Results: Numbers ##
######################

# Table 01
# points
ObtainResults(res.spde.points, field.true, !is.na(mesh.df$loc))
# 47 areas
ObtainResults(res.spde.area, field.true, !is.na(mesh.df$loc))
# Both
ObtainResults(res.spde.both, field.true, !is.na(mesh.df$loc))
# 8 areas
ObtainResults(res.spde.area8, field.true, !is.na(mesh.df$loc))
# Both
ObtainResults(res.spde.both8, field.true, !is.na(mesh.df$loc))
