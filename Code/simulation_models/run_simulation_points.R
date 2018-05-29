## Points

stk.dat.coords <- inla.stack(data=list(y=kenya.data$ybar,
                                       N=kenya.data$N),
                             A=list(A, 1), 
                             effects=list(i=1:mesh.true$n, 
                                          alpha=rep(1, nrow(A))),
                             tag='dat')

stk.all <- inla.stack(stk.dat.coords, stkgrid)

# 
# 0 < sigma2 < 1 --> assign it a beta(2,5) prior:
prior.function = function(log_precision) {
  log(30) - log_precision*(2+4) + 4*log(exp(log_precision)-1)
}

lprec = seq(0.001, 6, length.out = 1000)
prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
                    sep = "", collapse = " ")


res.spde.points.pr <- inla(f.spde,
                           control.compute=list(dic=TRUE),
                           data=inla.stack.data(stk.all),
                           control.family=list(
                             hyper=list(
                               prec=list(
                                 prior = prior.table,
                                 initial=2))),
                           control.predictor=list(
                             A=inla.stack.A(stk.all), compute=TRUE),
                           scale=N,  # scale is on precision scale
                           control.fixed=list(mean=list(alpha=0),
                             prec=list(alpha=1/100)))

res.spde.points <- inla(f.spde,
                        control.compute=list(dic=TRUE),
                        data=inla.stack.data(stk.dat.coords),
                        control.family=list(
                          hyper=list(
                            prec=list(
                              prior = prior.table,
                              initial=2
                            ))),
                        control.predictor=list(
                          A=inla.stack.A(stk.dat.coords), compute=TRUE),
                        scale=N,  # scale is on precision scale
                        control.fixed=list(mean=list(alpha=0),
                                           prec=list(alpha=1/100)))
