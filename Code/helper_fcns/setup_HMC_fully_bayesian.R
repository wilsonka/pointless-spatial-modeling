# Fully Bayesian HMC

# Y|S, beta0 ~ Pois(E A exp(S + beta0))
#   S ~ N(0, Q^{-1})

library(mvtnorm)

## U function
U <- function(q, Y, E, A, Q) {
  ## Inputs: values for the latent variables (q),
  ##         obs (Y), expected countes (E), projector matrix (A),
  ##         Q -- the precision matrix
  ## Outputs: value of -log posterior evaluated at q
  alpha <- q[1]
  S <- q[2:length(q)]
  
  -alpha*sum(Y) - t(Y) %*% log(A %*% exp(S)) + t(E) %*% A %*% exp(alpha + S) + 1/2*t(S) %*% Q %*% S + 1/200 * alpha^2
}

## grad_U function
grad_U <- function(q, Y, E, A, Q) {
  ## Inputs: values for the latent variables (q),
  ##         obs (Y), expected countes (E), projector matrix (A),
  ##         Q -- the precision matrix
  ## Outputs: vector that is of length q, dU/dalpha, dU/S
  alpha <- q[1]
  S <- q[2:length(q)]
  
  dU_dalpha <-  as.vector(-sum(Y) + exp(alpha) * t(E) %*% A %*% exp(S) + 
                            alpha/100)
  Astar <- t(t(A) * exp(S)) # this gives a matrix dim(A) where each entry is Astar[i,j] = A[i,j]*exp(S)
  dU_dS <- as.vector(t(Astar) %*% (E*exp(alpha)-Y/(A%*%exp(S)))) + 
    as.vector(t(S) %*% Q)
  c(dU_dalpha, dU_dS)
}

HMC <- function(stepsize, L, current_q, Y, E, A, Q) {
  ## Inputs: step size, number of steps in each iteration (L), 
  ##    current values of latent variables (current_q), 
  ##    precision matrix for S (Q)
  
  eps <- stepsize*runif(1, .9, 1.1)
  
  q <- current_q
  p <- rnorm(length(q), 0, sd=1) # ind normals
  
  current_p <- p
  
  p <- p - eps/2 * grad_U(q, Y, E, A, Q)
  
  for (i in 1:L) {
    q <- q + eps * p
    if (i != L) {
      p <- p - eps * grad_U(q, Y, E, A, Q)
    }
  }
  
  p <- p - eps/2 * grad_U(q, Y, E, A, Q)
  p <- -p
  
  current_U <- U(current_q, Y, E, A, Q)
  current_K <- sum(current_p^2)/2
  proposed_U <- U(q, Y, E, A, Q)
  proposed_K <- sum(p^2)/2
  
  if (log(runif(1)) < (current_U - proposed_U + current_K - proposed_K)[1]) {
    return (list(q=q,accept=1))
  } else {
    return (list(q=current_q,accept=0))
  }
}

steptheta <- function(q, theta, spde, c.Sig) {
  accept <- 0
  
  m <- spde$n.spde
  thetamu <- spde$param.inla$theta.mu
  S <- q[2:length(q)]
  theta1 <- theta[1]
  theta2 <- theta[2]
  
  thetastar <- as.vector(rmvnorm(1, mean=theta, sigma=c.Sig))
  theta1star <- thetastar[1]
  theta2star <- thetastar[2]
  Q <- inla.spde.precision(spde, theta=theta)
  Qstar <- inla.spde.precision(spde, theta=thetastar)
  
  posteriorstar <- inla.qsample(Q=Qstar, mu=rep(0, m),
                                sample=S, logdens=T)$logdens -
    1/20 * (theta1star - thetamu[1])^2 - 1/20 * (theta2star - thetamu[2])^2
  posterior <- inla.qsample(Q=Q, mu=rep(0, m),
                            sample=S, logdens=T)$logdens -
    1/20 * (theta1 - thetamu[1])^2 - 1/20 * (theta2 - thetamu[2])^2
  
  if(log(runif(1)) < (posteriorstar - posterior) ) {
    accept <- 1
    theta <- thetastar
  }
  
  return(list(accept=accept, theta=theta))
}

runManyFB <- function(M, stepsize, L, initial.q, initial.theta , Y, E, A, spde,
                      c.Sig, maxL=NULL, thin=1, burnin=1000) {
  qs <- matrix(0, ncol=length(initial.q), nrow=((M-burnin)/thin+1))
  qs[1,] <- initial.q
  q.prev <- initial.q

  thetas <- matrix(0, ncol=2, nrow=((M-burnin)/thin + 1))
  thetas[1,] <- initial.theta
  theta.prev <- initial.theta

  accepts <- matrix(0, nrow=((M-burnin)/thin), ncol=2)

  j <- 1
  for (i in 1:M) {
    if(!is.null(maxL)){
      L <- sample(1:maxL,1)
    }
    if(i %% 100 == 0) {
      print(i)
    }

    outstep1 <- steptheta(q.prev, theta.prev, spde, c.Sig)
    theta.prev <- outstep1$theta
    Q <- inla.spde.precision(spde, theta=theta.prev)
    outstep2 <- HMC(stepsize, L, q.prev, Y, E, A, Q)
    q.prev <- outstep2$q
    if(i %% thin == 0 && i > burnin) {
      qs[j+1, ] <- q.prev
      thetas[j+1, ] <- theta.prev
      accepts[j, ] <- c(outstep1$accept, outstep2$accept)
      j <- j+1
    }
  }
  print(paste("Acceptance theta:",mean(accepts[,1])))
  print(paste("Acceptance others:",mean(accepts[,2])))

  #list(qs=qs, thetas=thetas)
  cbind(qs[,1],thetas,qs[,2:ncol(qs)])
}
