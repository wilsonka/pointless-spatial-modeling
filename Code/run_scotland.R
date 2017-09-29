# Run Scotland Example

#########################
########## ICAR #########
#########################
source("Code/scotland_models/run_scotland_icar.R")

#########################
#### Empirical Bayes ####
#########################
source("Code/scotland_models/run_scotland_eb.R")

#########################
######### Hybrid ########
#########################
source("Code/scotland_models/run_scotland_hybrid.R")

#########################
##### Fully Bayesian ####
#########################
source("Code/scotland_models/run_scotland_fully_bayesian.R")


#########################
# Convergence Assessment
#########################

## Rhat values

# Hybrid
x1 <- mcmc(runc1[-1, ])
x2 <- mcmc(runc2[-1, ])
x3 <- mcmc(runc3[-1, ])
x4 <- mcmc(runc4[-1, ])
r.coda <- gelman.diag(mcmc.list(x1, x2, x3, x4), autoburnin=F, multivariate=F)
r.coda$psrf[1,]
apply(r.coda$psrf[2:nrow(r.coda$psrf), ], 2, summary)

# Fully Bayesian
x1 <- mcmc(run.update.theta1[-1, ])
x2 <- mcmc(run.update.theta2[-1, ])
x3 <- mcmc(run.update.theta3[-1, ])
x4 <- mcmc(run.update.theta4[-1, ])
r.coda <- gelman.diag(mcmc.list(x1, x2, x3, x4), autoburnin=F, multivariate=F)
r.coda$psrf[1:3, ]
apply(r.coda$psrf[4:nrow(r.coda$psrf), ], 2, summary)

#########################
##### Summary stats #####
#########################

## ICAR

# exp(beta0)
inla.zmarginal(inla.tmarginal(function(x) exp(x),  
                              res.icar$marginals.fixed$`(Intercept)`))
# tau
signif(res.icar$summary.hyperpar, 3)

## EB

# delta method

# estimate for exp(beta0)
signif(exp(Opt.eb$par[1]), 3)
signif(exp(Opt.eb$par[1] + qnorm(c(0.025, 0.975))* sqrt(cov.eb[1, 1])), 3)

# on log scale for phi and lambda2
cov.eb.log.sp.parms <- matrix(c(-2, 0, -2, -1), nrow=2) %*% 
  cov.eb[2:3, 2:3] %*% t(matrix(c(-2, 0, -2, -1), nrow=2))

# estimate for phi
sqrt(8)/exp(Opt.eb$par[3])
exp(log(sqrt(8)) - Opt.eb$par[3] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.eb.log.sp.parms[2,2]))

# estimate for lambda^2
1/(4*pi*exp(2*Opt.eb$par[3] + 2*Opt.eb$par[2]))
exp(-log(4*pi) - 2 * Opt.eb$par[3] - 2 * Opt.eb$par[2] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.eb.log.sp.parms[1,1]))


## Hybrid

# exp(beta0)
signif(quantile(exp(run.fixed.theta[, 1]), probs=c(0.025, 0.5, 0.975)), 3)
# theta1 = log(tau)
Opt.hybrid$par[1]
# theta2 = log(kappa)
Opt.hybrid$par[2]

# on log scale for phi and lambda2
cov.hybrid.log.sp.parms <- matrix(c(-2, 0, -2, -1), nrow=2) %*% 
  cov.hybrid %*% t(matrix(c(-2, 0, -2, -1), nrow=2))


# phi
kappa.fixed.theta <- exp(Opt.hybrid$par[2])
signif(sqrt(8)/kappa.fixed.theta, 3)

exp(log(sqrt(8)) - Opt.hybrid$par[2] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.hybrid.log.sp.parms[2,2]))


# lambda2 = 1/4pik^2t^2
tau.fixed.theta <- exp(Opt.hybrid$par[1])
signif(1/(4*pi*kappa.fixed.theta^2*tau.fixed.theta^2), 3)

exp(-log(4*pi) - 2 * Opt.hybrid$par[2] - 2 * Opt.hybrid$par[1] + qnorm(c(0.025, 0.975)) * 
      sqrt(cov.hybrid.log.sp.parms[1,1]))

## Fully Bayesian

# exp(beta0)
signif(quantile(exp(run.update.theta[, 1]), probs=c(0.025, 0.5, 0.975)), 3)
# theta1 = log(tau)
signif(quantile(run.update.theta[, 2], probs=c(0.025, 0.5, 0.975)), 3)
# theta2 = log(kappa)
signif(quantile(run.update.theta[, 3], probs=c(0.025, 0.5, 0.975)), 3)
# phi
kappas <- exp(run.update.theta[, 3])
signif(quantile(sqrt(8)/kappas, probs=c(0.025, 0.5, 0.975)), 3)
# lambda2 = 1/4pik^2t^2
taus <- exp(run.update.theta[, 2])
signif(quantile(1/(4*pi*kappas^2*taus^2), probs=c(0.025, 0.5, 0.975)), 3)


