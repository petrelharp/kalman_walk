source("perturb_Kalman.r")

#####################
# COMPARE minimal vs non-minimal system #

min0 <- min_perturb(tau=0)
nonmin0 <- sys_gen2(s=0.01, eps=0.1, nreps=100)

pdf("~/kalman_walk/sims/2d_vs_6d_oscillator_tau0.pdf")
plot(min0, col="blue")
lines(lowess(min0), col="blue")
points(nonmin0, col="black")
lines(lowess(nonmin0), col="black")
legend(0.1,4,legend=c("2D minimal system", "6D system"), col=c("blue", "black"), pch=1)
dev.off()

##############

