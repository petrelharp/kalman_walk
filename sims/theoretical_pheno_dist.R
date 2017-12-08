A <- matrix(c(0,-1,1,0), nrow=2);
B <- matrix(c(1,1), ncol=1);
C <- matrix(c(1,0), nrow=1);

V <- function (tau) {
  matrix(c(1,tau,0,1-tau), nrow=2) 
}

S <- function(tau) {
  V(tau)%*%A%*%solve(V(tau))
}

gamete_gen <- function(tau, eps) {
  gametes <- list()
  count <- 0
  for (i in 0:1)
  {
    for (j in 0:1)
    {
      for (k in 0:1)
      {
        for (l in 0:1)
        {
          count <- count + 1
          M <- matrix(c(i,k,j,l), ncol=2, nrow=2)
          gametes[[count]] <- S(tau)*M + S(tau + eps)*(1 - M)
        }
      }
    }
  }
  return(gametes)
}

progeny_gen <- function(tau, eps){
  gametes <- gamete_gen(tau, eps)
  progeny <- list()
  count2 <- 0
  for (i in 1:16)
  {
    for (j in 1:16)
    {
      count2 <- count2 + 1
      progeny[[count2]] <- (gametes[[i]] + gametes[[j]])/2
    }
  }
  return(progeny)
}

U = matrix(c(sqrt(1/2), 0+1i/sqrt(2), sqrt(1/2), 0-1i/sqrt(2)), ncol=2)
invU = matrix(c(sqrt(1/2), sqrt(1/2),0-1i/sqrt(2), 0+1i/sqrt(2)), nrow=2)

theory_comp <- function(eps) {
  progeny <- progeny_gen(0, eps)
  mean_dist <- 0
  for (k in 1:256) 
  {
    X0 <- invU%*%progeny[[k]]%*%U
    f <- function (t) {
        sapply(t, function (tt) {
              X <- UZU * matrix( c( sin(tt), tt*exp(tt*(0+1i)), sin(tt), tt*exp(tt*(0-1i)) ), nrow=2)
              return( (C%*%U%*%X%*%invU%*%B)^2 )
           })
    }
    mean_dist <- mean_dist + 0.5*(eps^2)*integrate(f, upper = Inf, lower = 0)
  }
  return(mean_dist)
}
