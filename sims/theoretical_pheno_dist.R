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
          name <- paste("G", count, sep = ""); 
          assign(name, S(tau)*M + S(tau + eps)*(1 - M))
          gametes[[name]] <- get(name)
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
  for (k in k:256) 
  {
    X <- invU%*%progeny[[k]]%*%U * matrix( c( sin(t), t*exp(t*(0+1i)), sin(t), t*exp(t*(0-1i)) ), nrow=2)
    ingd <- C%*%U%*%X%*%invU%*%B%*%t(B)%*%t(invU)%*%t(X)%*%t(U)%*%t(C)
    mean_dist <- mean_dist + 0.5*(eps^2)*integrate(ingd, upper = 10, lower = 0)
  }
  return(mean_dist)
}