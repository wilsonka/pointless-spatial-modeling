## Fully Bayesian Approach

initial.q <- c(Opt.eb$par[1], Obj.eb$env$parList()[3]$S)
initial.theta <- Opt.eb$par[2:3]
Y <- scot.dat$Observed
E <- scot.dat$Expected
A <- D
spde <- spde
c.Sig <- 0.05 * cov.eb[2:3, 2:3]

if(run_fb) {
  run.update.theta1 <- runManyFB(M=1010000, stepsize=0.01, L=10, initial.q,
                                 initial.theta , Y, E, A, spde, (2/5*c.Sig), maxL=NULL,
                                 thin=1000, burnin=10000)
  
  run.update.theta2 <- runManyFB(M=560000, stepsize=0.008, L=20, initial.q,
                                 initial.theta , Y, E, A, spde, (2.5/5*c.Sig), maxL=NULL,
                                 thin=550, burnin=10000)
  
  run.update.theta3 <- runManyFB(M=510000, stepsize=0.01, L=13, initial.q,
                                 initial.theta , Y, E, A, spde, (2/5*c.Sig), maxL=NULL,
                                 thin=500, burnin=10000)
  
  run.update.theta4 <- runManyFB(M=1010000, stepsize=0.01, L=10, initial.q,
                                 initial.theta , Y, E, A, spde, (2/5*c.Sig), maxL=NULL,
                                 thin=1000, burnin=10000)

} else {
  load("Results/scotland_fb_results.RData")
}

run.update.theta <- rbind(run.update.theta1[-1, ], run.update.theta2[-1, ],
                          run.update.theta3[-1, ], run.update.theta4[-1, ])
