####################

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

########
h <- function (t,tau=0,A=A_fn(tau)) { sapply(t, function (tt) C %*% expm(tt*A) %*% B) }
optimal_h <- function (t) { sin(t)+cos(t) }
Here is the impulse response function:
  
  tt <- seq(0,2*pi, length.out=100)
plot(tt, h(tt), type='l',
     xlab='t', ylab=expression(h(t)), main='impulse response function')
#lines(tt, h(tt, tau=0.1), col='blue', lty=3)
#lines(tt, optimal_h(tt), col='red', lty=2)
######

pdyn <- function(t, A) { sapply(t, function (tt) C %*% expm(tt*A) %*% B) }
tt <- seq(0,10*pi, length.out=1000)
plot(tt, pdyn(tt, A=A), type='l',
     xlab='t', ylab=expression(h(t)), main='Phenotype', xlim=c(0,30), ylim=c(-6,6))
for (i in 1:256){
lines(tt, pdyn(tt, A=progeny[[i]]), col=sample(256,1), lty=1)
}
