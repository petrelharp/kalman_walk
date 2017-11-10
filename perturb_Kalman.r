sys_gen2 <- function (s=sigma) 
{
  A <- matrix(c(0,-1,1,0), nrow=2);
  B <- matrix(c(1,1), ncol=1);
  C <- matrix(c(1,0), nrow=1);
  A00 <- matrix(c(0,0,0,0), nrow=2);
  B00 <- matrix(c(0,0), ncol=1);
  C00 <- matrix(c(0,0), nrow=1);
  #A14 <- matrix(c(rnorm(4,mean=0, sd=1)),nrow=2);
  A12 <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  Arno <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  A13 <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  #A24 <- matrix(c(rnorm(4,mean=0, sd=1)),nrow=2);
  Anrno <- matrix(c(rnorm(4,mean=0, sd=s)),nrow=2);
  
  p <- rnorm(25,mean=0, sd=s);
  P <- matrix(c(1,1-p[1],1-p[2],p[3],p[4],p[5],
                0,p[1],p[2],-p[3],-p[4],-p[5],
                0,0,p[7],p[8],p[9],p[10],
                0,0,p[12],p[13],p[14],p[15],
                0,0,p[17],p[18],p[19],p[20],
                0,0,p[22],p[23],p[24],p[25]), nrow=6);
  invP <- solve(P);
  
  AK <- (rbind(cbind(A,A00,A00),
               cbind(A12,Arno,A13),
               cbind(A00,A00, Anrno)));
  BK <- matrix(c(1,1,0,0,0,0),ncol=1);
  CK <- matrix(c(1,0,0,0,0,0),nrow=1);
  KalC <- CK%*%P;
  KalB <- P%*%BK;
  KalA <- P%*%AK%*%invP

output <- matrix(0,nrow=1000,ncol=2)
for (i in 1:1000)
{
nam0 <- paste("new_KalA", i, sep = ""); 

mA12 <- matrix(c(rnorm(4, mean=0, sd=0.01)),nrow=2,ncol=2)
mArno <- matrix(c(rnorm(4, mean=0, sd=0.01)),nrow=2,ncol=2)
mA13 <- matrix(c(rnorm(4, mean=0, sd=0.01)),nrow=2,ncol=2)
mAnrno <- matrix(c(rnorm(4, mean=0, sd=0.01)),nrow=2,ncol=2)
mp <- rnorm(25,mean=0,sd=0.01)

new_A12 <- A12 + A12*mA12
new_Arno <- Arno + Arno*mArno
new_A13 <- A13 + A13*mA13
new_Anrno <- Anrno #+ Anrno*mAnrno
new_p <- p + p*mp
mP <- matrix(c(1,1-new_p[1],1-new_p[2],new_p[3],new_p[4],new_p[5],
               0,new_p[1],new_p[2],-new_p[3],-new_p[4],-new_p[5],
               0,0,new_p[7],new_p[8],new_p[9],new_p[10],
               0,0,new_p[12],new_p[13],new_p[14],new_p[15],
               0,0,new_p[17],new_p[18],new_p[19],new_p[20],
               0,0,new_p[22],new_p[23],new_p[24],new_p[25]), nrow=6);
invmP <- solve(mP);
new_AK <- (rbind(cbind(A,A00,A00),
                 cbind(new_A12,new_Arno,new_A13),
                 cbind(A00,A00, new_Anrno)))
temp_KalA <- mP%*%new_AK%*%invmP

assign(nam0, temp_KalA)

dist <- norm(KalA-get(nam0),"F")

avg_div <- 0
for(j in 1:100) 
{ 
  nam <- paste("Z", j, sep = "") 
  assign(nam, recombine(KalA,get(nam0)))
  
  nam2 <- paste("DZ", j, sep = "")
  assign(nam2, D(get(nam)))
  
  avg_div <- avg_div + get(nam2)/100
}
output[i,] <- c(dist,avg_div)
}
output <- output
}