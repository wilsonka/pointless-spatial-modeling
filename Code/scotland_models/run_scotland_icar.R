## ICAR model (no non-spatial term)

data.inla <- data.frame(obs=scot.dat$Observed, exp=scot.dat$Expected)
scot.nbhd <- read.csv("Data/Scotland/scot_neighborhood.csv", header=F)
scot.nbhd <- scot.nbhd[,-(1:2)]  # first entry is the region 
                                 # second entry is # of neighbors

# make neighborhood object -- matrix of 0s and 1s
scot.nbhd.mat <- matrix(0, nrow=56, ncol=56)
for(i in 1:56) {
  scot.nbhd.mat[i, scot.nbhd[i, ][!is.na(scot.nbhd[i, ])]] <- 1
}

data.inla$region.struct <- 1:nrow(data.inla)
res.icar <- inla(obs ~ f(region.struct, model="besag", graph=scot.nbhd.mat,
                         adjust.for.con.comp = T, param=c(0.5,0.0005)), 
                 E=exp, family="poisson", data=data.inla,
                 control.predictor=list(compute=TRUE),
                 control.compute=list(dic=TRUE))
summary(res.icar)

obs.vals.icar <- data.frame(SIR.icar=res.icar$summary.fitted.values$`0.5quant`,
                       SIR.icar0.025=res.icar$summary.fitted.values$
                         `0.025quant`,
                       SIR.icar0.975=res.icar$summary.fitted.values$
                         `0.975quant`)
obs.vals.icar$id <-  scot.dat$District
obs.vals.icar <- obs.vals.icar[order(obs.vals.icar$id), ]

scot.df.icar <- join(obs.vals.icar, scot.df, by = "id")

