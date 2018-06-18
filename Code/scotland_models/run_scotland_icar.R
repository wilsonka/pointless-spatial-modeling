## ICAR model (with non-spatial term)

data.inla <- data.frame(obs=scot.dat$Observed, exp=scot.dat$Expected)

data.inla$region.struct <- 1:nrow(data.inla)
data.inla$region.unstruct <- 1:nrow(data.inla)

res.icar <- inla(obs ~ f(region.struct, model="besag", graph="Data/Scotland/scotland.graph",
                         adjust.for.con.comp = T, param=c(1, 0.2/0.59),
                         scale.model=T) +
                   f(region.unstruct, model="iid", param=c(1, 0.14)), 
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

