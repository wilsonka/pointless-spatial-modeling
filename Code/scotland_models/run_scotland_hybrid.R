## Hybrid Approach

## Step 1: Get EB estimates for theta
Obj.hybrid <- MakeADFun(data=Data.tmb, parameters=Params.tmb,
                        DLL="U_hybrid", random=c("alpha", "S_j"), silent=T)

Opt.hybrid <- optim(par=Obj.hybrid$par, fn=Obj.hybrid$fn,
                    gr=Obj.hybrid$gr, method="BFGS", hessian=T)

#Opt.hybrid$par # 2.39, -3.35

cov.hybrid <- solve(Opt.hybrid$hessian)


## Step 2: Run HMC for fixed theta
init.q <- c(Obj.hybrid$env$parList()[3]$alpha,
            Obj.hybrid$env$parList()[2]$S + Obj.hybrid$env$parList()[3]$alpha)
Q <- inla.spde2.precision(spde, theta=Opt.hybrid$par)

set.seed(98115)
system.time(run.initial <- runManyHybrid(250, c(0.01, 0.02), c(50, 50), init.q,
                                   scot.dat$Observed, scot.dat$Expected, D, Q,
                                   var.p=rep(1,length(init.q)), thin=1,
                                   burnin=50))
posterior.vars <- apply(run.initial, 2, var)


set.seed(2987)
system.time(runc1 <- runManyHybrid(1500, c(0.02, 0.07), c(50, 75), init.q,
                             scot.dat$Observed, scot.dat$Expected, D, Q,
                             var.p=1/posterior.vars, thin=1, burnin=500))

system.time(runc2 <- runManyHybrid(1500, c(0.02, 0.072), c(50, 80), init.q,
                             scot.dat$Observed, scot.dat$Expected, D, Q,
                             var.p=1/posterior.vars, thin=1, burnin=500))

system.time(runc3 <- runManyHybrid(1500, c(0.02, 0.072), c(50, 90), init.q,
                             scot.dat$Observed, scot.dat$Expected, D, Q,
                             var.p=1/posterior.vars, thin=1, burnin=500))

system.time(runc4 <- runManyHybrid(1500, c(0.02, 0.071), c(50, 90), init.q,
                             scot.dat$Observed, scot.dat$Expected, D, Q,
                             var.p=1/posterior.vars, thin=1, burnin=500))

run.fixed.theta <- rbind(runc1[-1, ], runc2[-1, ], runc3[-1, ], runc4[-1, ])

