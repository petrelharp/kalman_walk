library(expm)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)

pdf(file="blue_hybrids.pdf", width=6, height=4, pointsize=10)

A0 <- matrix(c(0,-1,1,0), nrow=2)
B <- matrix(c(1,1), ncol=1)
C <- matrix(c(1,0), nrow=1)

A_fn <- function (tau) { 
    matrix(c(tau/(tau+1), -(2*tau*(tau+1)+1)/(tau+1), 
             1/(tau+1), -tau/(tau+1)), nrow=2) 
}

tt <- seq(0,6*2*pi, length.out=500)

h <- function (t,tau=0,A=A_fn(tau)) { sapply(t, function (tt) C %*% expm(tt*A) %*% B) }


F1 <- function (tau1,tau2) { (A_fn(tau1)+A_fn(tau2))/2 }

get_F2s <- function(tau1, tau2) 
{
    vc <- function(tau1, tau2, x, y) 
    {
        out <- c(A_fn(tau1)[x,y], A_fn(tau2)[x,y])
        return(out)
    }
    mlist <- list()
    count <- 1
    for (i in 1:2) 
    {
        for (j in 1:2) 
        {
            for (k in 1:2) 
            {
                for (l in 1:2)
                {
                    mlist[[count]] <- matrix(c(vc(tau1, tau2, 1, 1)[i], vc(tau1, tau2, 2, 1)[j], vc(tau1, tau2, 1, 2)[k], vc(tau1, tau2, 2, 2)[l]),2,2)
                    count <- count + 1
                }
            }
        }
    }
    f2list <- list()
    count2 <- 1
    for(i in 1:length(mlist)) 
    {
        for (j in 1:length(mlist))
        {
            f2list[[count2]] <- (mlist[[i]] + mlist[[j]])/2
            count2 <- count2 + 1
        }
    }
    return(unique(f2list))
}


hhf0 <- sapply( list(A_fn(0), A_fn(0)), function (A) { sapply(tt, h, A=A) } )

hhf1_1 <- sapply( list(A_fn(0), A_fn(0.01), F1(0, 0.01)), function (A) { sapply(tt, h, A=A) } )
hhf1_2 <- sapply( list(A_fn(0), A_fn(0.1), F1(0, 0.1)), function (A) { sapply(tt, h, A=A) } )
hhf1_3 <- sapply( list(A_fn(0), A_fn(0.5), F1(0, 0.5)), function (A) { sapply(tt, h, A=A) } )

hhf2_1 <- sapply( get_F2s(0, 0.01), function (A) { sapply(tt, h, A=A) } )
hhf2_2 <- sapply( get_F2s(0, 0.1), function (A) { sapply(tt, h, A=A) } )
hhf2_3 <- sapply( get_F2s(0, 0.5), function (A) { sapply(tt, h, A=A) } )



cf <- colorRampPalette(colors = c("grey", "steelblue"))
#cf <- colorRampPalette(colors = palette())
hx <- data.frame("value"=hhf0[,1], "time"=tt, "variable"="parent")

#F1s 

k <- data.frame(hhf1_1)
k$time <- tt
pl11 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(3)) + ylim(c(-6,6)) + xlab("") + ylab("") +
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.01")), x=0, y=6, hjust=0, vjust=1, size=6)

k <- data.frame(hhf1_2)
k$time <- tt
pl12 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(3)) + ylim(c(-6,6)) + xlab("") + ylab("") + 
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.1")), x=0, y=6, hjust=0, vjust=1, size=6)

k <- data.frame(hhf1_3)
k$time <- tt
pl13 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(3)) + ylim(c(-6,6)) + xlab("") + ylab("") +
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.5")), x=0, y=6, hjust=0, vjust=1, size=6)

#F2s 
k <- data.frame(hhf2_1)
k$time <- tt
pl21 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(81)) + ylim(c(-6,6)) + xlab("") + ylab("") +
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.01")), x=0, y=6, hjust=0, vjust=1, size=6)

k <- data.frame(hhf2_2)
k$time <- tt
pl22 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(81)) + ylim(c(-6,6)) + xlab("") + ylab("") +
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.1")), x=0, y=6, hjust=0, vjust=1, size=6)

k <- data.frame(hhf2_3)
k$time <- tt
pl23 <- (melt(k, id.vars = "time")) %>% ggplot(aes(x=time, y=value, group=variable, col=variable)) + 
    geom_line() + theme_classic(12) + theme(legend.position = "none") + 
    scale_color_manual(values = cf(81)) + ylim(c(-6,6)) + xlab("") + ylab("") +
    geom_line(data=hx, aes(x=time, y=value), color="black", size=1) + 
    annotate(geom = "text", label=expression(paste(epsilon, " = 0.5")), x=0, y=6, hjust=0, vjust=1, size=6)

tr.grob <- textGrob(bquote(F[2]), 
                    gp=gpar(fontface="bold", col="black", fontsize=24))
tl.grob <- textGrob(bquote(F[1]), 
                    gp=gpar(fontface="bold", col="black", fontsize=24))
x.grob <- textGrob("time", 
                   gp=gpar(col="black", fontsize=18))
y.grob <- textGrob("phenotype", 
                   gp=gpar(col="black", fontsize=18), rot = 90)

k1 <- plot_grid(pl11, pl12, pl13, ncol=1)
k2 <- plot_grid(pl21, pl22, pl23, ncol=1)

k11 <- grid.arrange(k1, top=tl.grob)
k22 <- grid.arrange(k2, top=tr.grob)

k0 <- plot_grid(k11, k22, nrow = 1)
k00 <- grid.arrange(k0, bottom=x.grob, left=y.grob)

dev.off()
