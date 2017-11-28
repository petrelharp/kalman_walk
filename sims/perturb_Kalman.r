library(MASS)
library(Matrix)
library(expm)

#####################################
h <- function (t, M, BK, CK) { sapply(t, function (tt) CK %*% expm(tt*M) %*% BK) }
optimal_h <- function (t) { sin(t)+cos(t) }
###################################
Df <- function (A,t, BK, CK) {
  exp(-t/(4*pi)) * ( h(t, M=A, BK=BK, CK=CK) - optimal_h(t) )^2
}
####################################
D <- function (A, BK, CK, upper=10, ...) {
  f <- function (t) { Df(A,t, BK, CK) }
  integrate(f, lower=0, upper=upper, ...)$value
}

##################################
recombine<- function (X,Y) 
{
  n <- dim(X)[1];
  r <- matrix(c(rbinom(n^2,1,0.5)),nrow=n);
  G <- X * r + Y * (1-r)
  return(G)
}

make_P <- function (p, perturb=FALSE) {
  # Make the matrix P from a vector of entries
  stopifnot(length(p)==25)
  x <- if(perturb){0} else {1}
  P <- matrix(c(x,x,x,0,0,0,
                0,x+p[1],p[2],p[3],p[4],p[5],
                0,0,x+p[7],p[8],p[9],p[10],
                0,0,p[12],x+p[13],p[14],p[15],
                0,0,p[17],p[18],x+p[19],p[20],
                0,0,p[22],p[23],p[24],x+p[25]), nrow=6);
  P[2:6,1] <- P[2:6,1] - P[2:6,2]
  return(P)
}

# Generate a starting matrix by adding noise 
# to the free parameters of the Kalman decomposition
# for the oscillator.
random_oscillator <- function(s) {
  A <- matrix(c(0,-1,1,0), nrow=2);
  B <- matrix(c(1,1), ncol=1);
  C <- matrix(c(1,0), nrow=1);
  #A14 <- matrix(c(rnorm(4,mean=0, sd=1)),nrow=2);
  A12 <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  Arno <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  A13 <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  #A24 <- matrix(c(rnorm(4,mean=0, sd=1)),nrow=2);
  Anrno <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  p <- rnorm(25,mean=0, sd=s);
  P <- make_P(p)
  return(list(A=A, A12=A12, Arno=Arno, Anrno=Anrno, A13=A13, B=B, C=C, P=P))
}

# assemble the Kalman system from a list of components
assemble_oscillator <- function (AL) {
    # AL should be a list with components:
    # A
    # A12
    # Arno
    # Anrno
    # A13
    # B
    # C
    # P
  A00 <- matrix(c(0,0,0,0), nrow=2);
  B00 <- matrix(c(0,0), ncol=1);
  C00 <- matrix(c(0,0), nrow=1);
  AK <- (rbind(cbind(AL$A, A00, A00),
               cbind(AL$A12, AL$Arno, AL$A13),
               cbind(A00, A00, AL$Anrno)));
  # pad B and C with zeros
  BK <- matrix(c(AL$B, rep(0, 4)), ncol=1)
  CK <- matrix(c(AL$C, rep(0, 4)), nrow=1)
  invP <- solve(AL$P);
  KalC <- CK%*%AL$P;
  KalB <- AL$P%*%BK;
  KalA <- AL$P%*%AK%*%invP
  return(list(A=KalA, B=BK, C=CK))
}

# Perturb an oscillator of the form above
perturb_oscillator <- function(AL, eps) {
    mA12 <- matrix(c(rnorm(4, mean=0, sd=eps)),nrow=2,ncol=2)
    mArno <- matrix(c(rnorm(4, mean=0, sd=eps)),nrow=2,ncol=2)
    mA13 <- matrix(c(rnorm(4, mean=0, sd=eps)),nrow=2,ncol=2)
    mAnrno <- matrix(c(rnorm(4, mean=0, sd=eps)),nrow=2,ncol=2)
    mp <- rnorm(25,mean=0,sd=eps)
    
    new_A12 <- AL$A12 * (1 + mA12)
    new_Arno <- AL$Arno * (1 + mArno)
    new_A13 <- AL$A13 * (1 + mA13)
    new_Anrno <- AL$Anrno * (1 + mAnrno)
    mP <- AL$P + make_P(mp, perturb=TRUE)
    return(list(A=AL$A, A12=new_A12, Arno=new_Arno, Anrno=new_Anrno, 
                A13=new_A13, B=AL$B, C=AL$C, P=mP))
}

#####################################
sys_gen2 <- function (s=sigma, eps=0.05, nreps=100) 
{
  KalA_list <- random_oscillator(s=s)
  KalA <- assemble_oscillator(KalA_list)
  
  output <- matrix(0,nrow=nreps,ncol=2)
  colnames(output) <- c("sys_dist", "pheno_dist")
  for (i in 1:nreps)
  {
    
    new_KalA_list <- perturb_oscillator(KalA_list, eps=eps)
    new_KalA <- assemble_oscillator(new_KalA_list)
    dist <- norm(KalA$A-new_KalA$A,"F")
    
    avg_div <- 0
    jj <- 0
    for (j in 1:100) 
    { 
      Z1 <- recombine(KalA$A, new_KalA$A)
      Z2 <- recombine(KalA$A, new_KalA$A)
      Z3 <- (Z1 + Z2)/2
      DZ <- tryCatch(D(Z3, BK=KalA$B, CK=KalA$C), error=function(e) NULL)
      if (is.null(DZ)==1)
      {
        j = j-1
        jj = jj+1
        if(jj==10)
        {
          break
        }
      }
      if (is.null(DZ)==0) 
        {
          avg_div <- avg_div + DZ/100
        }
    }
    output[i,] <- c(dist,avg_div)
  }
  return(output)
}

#########################

pert_sim <- function(rgens){
  pheno_divergence <- sys_gen2(s=0.1, eps=0.05, nreps=100)
  for(i in 2:rgens)
  {
    ph0 <- sys_gen2(s=0.01,eps=0.05,nreps=100)
    pheno_divergence <- rbind(pheno_divergence, ph0)
  }
}
