# make a conceptual figure of a fitness landscape

dx <- 0.4
xx <- yy <- seq(1/sqrt(2) - dx,1/sqrt(2) + 0.7 * dx,length.out=100)
z <- abs(sqrt(outer(xx^2, yy^2, "+")) - 1.0)

plot_pop <- function (x, sigma, nreps=100, pch=20, cex=1, ...) {
    y <- matrix(rnorm(length(x)*nreps, mean=rep(x,each=nreps), sd=sigma), nrow=nreps)
    points(y, pch=pch, cex=cex, ...)
    return(invisible(y))
}

do.plot <- function (p1, p2, sigma, add=FALSE) {
    F1 <- (p1+p2)/2
    F2 <- rbind( c(p1[1],p2[2]), c(p2[1],p1[2]) )
    proj <- function (xy) { th <- atan2(x=xy[1],y=xy[2]); c(cos(th), sin(th)) }
    if (!add) contour(xx, yy, z, lty=3, 
                      xaxt='n', yaxt='n', xlab='', ylab='')
    lines(xx, sqrt(1-xx^2), lty=1, lwd=2)
    segments(x0=c(p1[1],p2[1]), x1=c(F1[1],F1[1]), 
             y0=c(p1[2],p2[2]), y1=c(F1[2],F1[2]), 
             lty=2, col='red')
    segments(x0=c(p1[1],p2[1]), x1=c(F2[2,1],F2[2,1]), 
             y0=c(p1[2],p2[2]), y1=c(F2[2,2],F2[2,2]), 
             lty=2, col='purple')
    segments(x0=F1[1], y0=F1[2], x1=proj(F1)[1], y1=proj(F1)[2], col='red')
    segments(x0=F2[2,1], y0=F2[2,2], x1=proj(F2[2,])[1], y1=proj(F2[2,])[2], col='purple')
    plot_pop(p1, sigma, col=adjustcolor('black', 0.2))
    plot_pop(p2, sigma, col=adjustcolor('black', 0.2))
    plot_pop(F1, sigma, col=adjustcolor('red', 0.2))
    plot_pop(F2[2,], sigma, col=adjustcolor('purple', 0.2))
    points(rbind(p1,p2,F1,F2[2,]), cex=2, pch=20,
           col=c("black","black","red","purple") )
    if (!add) legend("topright", pch=20, pt.cex=2, bg='white',
                   col=c("black","red","purple"),
                   legend=c("parental", "F1", "F2"))
}

p1 <- c(0.90, sqrt(1-0.90^2))
p2 <- c(0.5, sqrt(1-0.50^2))


pdf(file="conceptual_fig.pdf", width=3, height=3, pointsize=10)
par(mar=c(1,1,1,1))

do.plot(p1, p2, sigma=.02)

dev.off()


dx <- 0.25
xx <- yy <- seq(1/sqrt(2)-dx,1/sqrt(2)+dx,length.out=100)
z <- abs(sqrt(outer(xx^2, yy^2, "+")) - 1.0)


p3 <- c(0.75, sqrt(1-0.75^2))
p4 <- c(0.50, sqrt(1-0.50^2))


p5 <- c(0.78, sqrt(1-0.78^2))
p6 <- c(0.85, sqrt(1-0.85^2))



pdf(file="conceptual_fig_two_pairs.pdf", width=3, height=3, pointsize=10)
par(mar=c(1,1,1,1))

do.plot(p3, p4, sigma=.012)

do.plot(p5, p6, sigma=.012, add=TRUE)

dev.off()
