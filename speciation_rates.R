
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

plotf <- function (main, ...) {
    tt <- seq(0, 10, length.out=200)
    plot(tt, f(tt, ...), type='l', ylim=c(0,1), main=main,
            xlab=expression(T/N[e]), ylab="F2 fitness")
    lines(tt, f1(tt, ...), col='red', lty=3)
    lines(tt, f2(tt, ...), col='green', lty=3)
}

pdf(file="speciation_rates.pdf", width=7, height=3, pointsize=10)

layout(t(1:3))
plotf(a=1, b=1, main="large population")
plotf(a=.1, b=.1, main="small population")
plotf(a=10, b=1, main="metapopulation")

legend("topright", lty=c(1,3,3), col=c('black', 'red', 'green'),
        legend=c("total fitness", "mean", "variance"))

dev.off()


