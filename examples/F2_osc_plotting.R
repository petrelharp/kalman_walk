####################
Df <- function (A,t, BK, CK) {
  exp(-t/(4*pi)) * ( h(t, M=A, BK=B, CK=C) - optimal_h(t) )^2
}
####################################
D <- function (A, BK, CK, upper=10, ...) {
  f <- function (t) { Df(A,t, B, C) }
  integrate(f, lower=0, upper=upper, ...)$value
}
#############################
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


##################

f2_divergence <- function(tau, eps) {
  progeny <- progeny_gen(tau,eps)
  avg_div <- 0
  for (i in 1:256){
    avg_div <- avg_div + sqrt(D(progeny[[i]]))/256
  }
  return(avg_div)
}

d <- seq(0,-0.5,-0.01)
df2 <- sapply(d, function(d) f2_divergence(0,d))

##########

f1_divergence <- function(tau, eps){
  f1div <- sqrt(D((S(tau) + S(tau+eps))/2))
  return(f1div)
}

df1 <- sapply(d, function(d) f1_divergence(0,d))

###############
gen_d <- sapply(d, function(d) norm(S(0)-S(0+d),"F"))

pdf("~/kalman_walk/examples/F2_vs_F1_divergence_tau0.pdf")
plot(df1~gen_d, ylim=c(0,2), type = "l", lwd=3, col="red", xlab="Genetic distance", ylab="Phenotypic divergence")
lines(df2~gen_d, lwd=3, col="blue")
legend(0.2,1.5, legend=c("F1 hybrids", "F2 hybrids"), col=c("red", "blue"), lty=1, lwd=3)
dev.off()
