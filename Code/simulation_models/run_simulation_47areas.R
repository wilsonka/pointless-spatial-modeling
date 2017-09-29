## 47 Areas

stk.dat.coords.area <- inla.stack(data=list(y=census.df47$ybar,
                                            N=census.df47$N),
                                  A=list(D, 1),
                                  effects=list(i=1:mesh.true$n, 
                                               alpha=rep(1, nrow(D))),
                                  tag='dat')
stk.all.area <- inla.stack(stk.dat.coords.area, stkgrid)

# 0 < sigma2 < 1 --> assign it a beta(2,5) prior:
prior.function = function(log_precision) {
  log(30) - log_precision*(2+4) + 4*log(exp(log_precision)-1)
}


lprec = seq(0.001, 6, length.out = 1000)
prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
                    sep = "", collapse = " ")

# res.spde.area.pr <- inla(f.spde,
#                          control.comput=list(dic=TRUE),
#                          data=inla.stack.data(stk.all.area),
#                          control.family=list(
#                            hyper=list(
#                              prec=list(
#                                prior = prior.table
#                              ))),
#                          control.predictor=list(A=inla.stack.A(stk.all.area),
#                                                 compute=TRUE),
#                          scale=N,  # scale is on precision scale
#                          control.fixed=list(mean=list(alpha=0),
#                            prec=list(alpha=1/100)))
# 
res.spde.area <- inla(f.spde,
                      control.comput=list(dic=TRUE),
                      data=inla.stack.data(stk.dat.coords.area),
                      control.family=list(
                        hyper=list(
                          prec=list(
                            prior = prior.table
                            ))),
                      control.predictor=list(
                        A=inla.stack.A(stk.dat.coords.area), compute=TRUE),
                      scale=N,  # scale is on precision scale
                      control.fixed=list(mean=list(alpha=0),
                                         prec=list(alpha=1/100)))
