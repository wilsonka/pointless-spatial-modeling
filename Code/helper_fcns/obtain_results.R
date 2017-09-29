ObtainResults <- function(res.spde, field.true, mesh.in.kenya) {
  # Feed in INLA results and get nice table as output
  #
  # Args:
  #   res.spde: INLA summary object
  #   field.true: true field to be compared to
  #   mesh.in.kenya: T/F vector of which mesh points are within Kenya
  
  # Intercept: beta0
  int <- paste0(signif(res.spde$summary.fixed$`0.5quant`, 3), " & (", 
                signif(res.spde$summary.fixed$`0.025quant`, 3), ", ",
                signif(res.spde$summary.fixed$`0.975quant`, 3), ")")
  
  # sigma^2
  sigma2.tmp <- inla.zmarginal(inla.tmarginal(function(x) 1/x,
                                           res.spde$marginals.hyperpar$`Precision for the Gaussian observations`),
                            silent=T)
  sigma2 <- paste0(signif(sigma2.tmp$quant0.5, 3), " & (", 
                signif(sigma2.tmp$quant0.025, 3), ", ",
                signif(sigma2.tmp$quant0.975, 3), ")")
  
  # Phi and lambda^2
  res.spde.field <- inla.spde2.result(res.spde, "i", spde, do.transf=TRUE)
  phi <- inla.qmarginal(p=c(0.5, 0.025, 0.975),
                        res.spde.field$marginals.range.nominal[[1]])
  phi <- paste0(signif(phi[1], 3), " & (", 
                signif(phi[2], 3), ", ", 
                signif(phi[3], 3), ")")
  
  lambda2 <- inla.qmarginal(p=c(0.5, 0.025, 0.975), 
                            res.spde.field$marginals.variance.nominal[[1]])
  lambda2 <- paste0(signif(lambda2[1], 3), " & (",
                    signif(lambda2[2], 3), ", ",
                    signif(lambda2[3], 3), ")")
  
  # MSE
  MSE <- signif(mean((res.spde$summary.random$i$mean[mesh.in.kenya] - 
                        field.true[mesh.in.kenya])^2), 3)
  # MAE
  MAE <- signif(mean(abs(res.spde$summary.random$i$mean[mesh.in.kenya] - 
                           field.true[mesh.in.kenya])), 3)
  
  rbind(int, sigma2, phi, lambda2, MSE, MAE)
}

