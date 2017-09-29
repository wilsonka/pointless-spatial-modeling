# Hybrid HMC

# Y|S.c ~ Pois(E A exp(S.c))
#   S.c ~ N(alpha, Q^{-1})


## U function
U <- function(q, Y, E, A, Q) {
  ## Inputs: values for the latent variables (q),
  ##         obs (Y), expected countes (E), projector matrix (A),
  ##         Q -- the precision matrix
  ## Outputs: value of -log posterior evaluated at q
  alpha <- q[1]
  S.c <- q[2:length(q)]
  
  - t(Y) %*% log(A %*% exp(S.c)) + t(E) %*% A %*% exp(S.c) + 
    1/2*t(S.c - alpha) %*% Q %*% (S.c - alpha) + 1/200 * alpha^2
}

## grad_U function
grad_U <- function(q, Y, E, A, Q) {
  ## Inputs: values for the latent variables (q),
  ##         obs (Y), expected countes (E), projector matrix (A),
  ##         Q -- the precision matrix
  ## Outputs: vector that is of length q, dU/dalpha, dU/S
  alpha <- q[1]
  S.c <- q[2:length(q)]
  
  dU_dalpha <- as.vector(alpha*sum(Q) - sum(Q %*% S.c) + alpha/100)
  
  # this gives a matrix dim(A) where each entry is Astar[i,j] = A[i,j]*exp(S)
  Astar <- t(t(A) * exp(S.c)) 
  
  dU_dS.c <- as.vector(t(Astar) %*% (E-Y/(A%*%exp(S.c)))) + 
    as.vector(t(S.c) %*% Q) - as.vector(alpha * colSums(Q))
  
  c(dU_dalpha, dU_dS.c)
}

HMC <- function (stepsize, L, current_q, Y, E, A, Q, var.p) {
  ## Inputs: step size, number of steps in each iteration (L), 
  ##    current values of latent variables (current_q), 
  ##    precision matrix for S (Q), and variance of momentum variables (var.p)
  
  eps <- stepsize*runif(1,.9,1.1)
  
  q <- current_q
  p <- rnorm(length(q),0,sd=sqrt(var.p)) # ind normals
  
  current_p <- p
  
  #
  p <- p - eps/2 * grad_U(q, Y, E, A, Q)
  
  for (i in 1:L) {
    q <- q + eps * 1/var.p * p
    if (i != L) {
      p <- p - eps * grad_U(q, Y, E, A, Q)
    }
  }
  

  p <- p - eps/2 * grad_U(q, Y, E, A, Q)
  p <- -p
  
  current_U <- U(current_q, Y, E, A, Q)
  current_K <- sum(current_p^2/var.p)/2
  proposed_U <- U(q, Y, E, A, Q)
  proposed_K <- sum(p^2/var.p)/2
  
  if (log(runif(1)) < (current_U - proposed_U + current_K - proposed_K)[1]) {
    return (list(q=q,accept=1))
  } else {
    return (list(q=current_q,accept=0))
  }
}

runManyHybrid <- function (M, stepsize, L, initial_q , Y, E, A, Q, var.p,
                     maxL=NULL, thin=1, burnin=1000) {
  qs <- matrix(0,ncol=length(initial_q),nrow=((M-burnin)/thin+1))
  qs[1,] <- initial_q
  q.prev <- initial_q
  accepts <- rep(0,((M-burnin)/thin))
  
  j <- 1
  for (i in 1:M) {
    #print(i)
    if(!is.null(maxL)){
      L <- sample(1:maxL,1)
    }
    
    if( i <= burnin) {
      out <- HMC(stepsize[1], L[1], q.prev, Y, E, A, Q, var.p)
      q.prev <- out$q
      
    } else {
      out <- HMC(stepsize[2], L[2], q.prev, Y, E, A, Q, var.p)
      q.prev <- out$q
      if(i %% thin == 0) {
        qs[j+1,] <- out$q
        accepts[j] <- out$accept
        j <- j+1
      }
    }
    
  }   
  print(paste("Acceptance:",mean(accepts)))
  
  qs=qs
}
