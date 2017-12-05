
f <- function (t, a, b, c_mu=1/20) {
    # a = sigma_N^2 / gamma^2
    # b = sigma_S^2 / gamma^2
    # omega = 1
    f1(t, a, b, c_mu) * f2(t, a, b, c_mu)
}

f2 <- function (t, a, b, c_mu=1/20) {
    1/sqrt(1 + 4 * a * t / (1+b))
}

f1 <- function (t, a, b, c_mu=1/20) {
    exp(-(c_mu * t)^2 / (1 + b + 4 * b * t))
}

plotf <- function (main, Ne=1.0, tmax=5,
        ...) {
    tt <- seq(0, tmax/Ne, length.out=200)
    plot(tt*Ne, f(tt, ...), type='l', ylim=c(0,1), main=main,
            xlab=if (Ne==1) { expression(T/N[e])} else { expression(T) }, 
            ylab="F2 fitness")
    lines(tt*Ne, f1(tt, ...), col='red', lty=3)
    # lines(tt*Ne, f2(tt, ...), col='green', lty=3)
}

pdf(file="speciation_rates.pdf", width=7, height=3, pointsize=10)

layout(t(1:3))
plotf(a=1, b=1, main="large population")
plotf(a=.1, b=.1, main="small population")
plotf(a=10, b=1, main="metapopulation")

legend("topright", lty=c(1,3), col=c('black', 'red', 'green'), bg='white',
        legend=c("F2 fitness", "F1 fitness"))

# legend("topright", lty=c(1,3,3), col=c('black', 'red', 'green'),
#         legend=c("total fitness", "shift of mean", "increase of variance"))

dev.off()



pdf(file="speciation_rates_2.pdf", width=7, height=6, pointsize=10)

layout(matrix(1:6, nrow=2))
plotf(a=1, b=1, main="large population")
plotf(a=1, b=1, main="large population, Ne=1e5", Ne=1e5, tmax=1e4)

plotf(a=.1, b=.1, main="small population")
plotf(a=.1, b=.1, main="small population, Ne=1e3", Ne=1e3, tmax=1e4)

plotf(a=10, b=1, main="metapopulation")
plotf(a=10, b=1, main="metapopulation, Ne=1e7", Ne=1e7, tmax=1e4)

legend("topright", lty=c(1,3,3), col=c('black', 'red', 'green'),
        legend=c("total fitness", "shift of mean", "increase of variance"))

dev.off()


