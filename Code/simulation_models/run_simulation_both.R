## Both

stk.dat.coords.both <- inla.stack(data=list(y=c(kenya.data$ybar,
                                                census.df47$ybar),
                                            N=c(kenya.data$N,
                                                census.df47$N)),
                                  A=list(D.both, 1), tag='dat',
                                  effects=list(i=1:mesh.true$n,
                                               alpha=rep(1, nrow(D.both))))

stk.all.both <- inla.stack(stk.dat.coords.both, stkgrid)

# 0 < sigma2 < 1 --> assign it a beta(2,5) prior:
prior.function = function(log_precision) {
  log(30) - log_precision*(2+4) + 4*log(exp(log_precision)-1)
}


lprec = seq(0.001, 6, length.out = 1000)
prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
                    sep = "", collapse = " ")

# res.spde.both.pr <- inla(f.spde,
#                          control.comput=list(dic=TRUE),
#                          data=inla.stack.data(stk.all.both),
#                          control.family=list(
#                            hyper=list(
#                              prec=list(
#                                prior = prior.table
#                              ))),
#                          control.predictor=list(A=inla.stack.A(stk.all.both),
#                                                 compute=TRUE),
#                          scale=N,  # scale is on precision scale
#                          control.fixed=list(mean=list(alpha=0),
#                                             prec=list(alpha=1/100)))  # alpha ~ N(0,100)
# 
res.spde.both <- inla(f.spde,
                      control.comput=list(dic=TRUE),
                      data=inla.stack.data(stk.dat.coords.both),
                      control.family=list(
                        hyper=list(
                          prec=list(
                            prior = prior.table
                          ))),
                      control.predictor=list(
                        A=inla.stack.A(stk.dat.coords.both), compute=TRUE),
                      scale=N,  # scale is on precision scale
                      control.fixed=list(mean=list(alpha=0),
                                         prec=list(alpha=1/100)))  # alpha ~ N(0,100)
