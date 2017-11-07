library(expm)

system_diff <- function (A0, B0, C0, A1, B1, C1, g = 1/(4*pi)) {
    # Evaluate
    # \int_0^\infty e^{-tg} |C0 e^{t A0} B0 - C1 e^{t A1} B1|^2 dt
    n0 <- nrow(A0)
    n1 <- nrow(A1)
    eig0 <- eigen(A0, symmetric=TRUE)
    eig1 <- eigen(A1, symmetric=TRUE)
    # this is the list of C v v^T B, where v is an eigenvector
    z0 <- lapply(1:ncol(eig0$vectors), function (k) { C0 %*% eig0$vectors[,k] %*% crossprod(eig0$vectors[,k], B0) })
    z1 <- lapply(1:ncol(eig1$vectors), function (k) { C1 %*% eig1$vectors[,k] %*% crossprod(eig1$vectors[,k], B1) })
    ZZ <- get_zmat(z0, z0)
    ZY <- get_zmat(z0, z1)
    YY <- get_zmat(z1, z1)
    LL <- 1/(g - outer(eig0$values, eig0$values, "+"))
    LN <- 1/(g - outer(eig0$values, eig1$values, "+"))
    NN <- 1/(g - outer(eig1$values, eig1$values, "+"))
    return(sum(ZZ*LL) + sum(ZY*LN) + sum(YY*NN))
}

get_zmat <- function (z0, z1) {
    # this is the matrix of tr(z[i]^T z[j])
    n0 <- length(z0)
    n1 <- length(z1)
    ZZ <- matrix(0, nrow=n0, ncol=n1)
    for (i in 1:n0) {
        for (j in 1:n1) {
            ZZ[i,j] <- sum(diag(crossprod(z0[[i]], z1[[j]])^2))
        }
    }
    return(ZZ)
}

system_diff_slow <- function (A0, B0, C0, A1, B1, C1, g = 1/(4*pi), upper=10) {

    h <- function (t,A,B,C) { sapply(t, function (tt) C %*% expm(tt*A) %*% B) }

    f <- function (t) {
        exp(-t/(4*pi)) * sum(( h(t,A0,B0,C0) - h(t,A1,B1,C1) )^2)
    }

    integrate(f, lower=0, upper=upper)$value
}

if (TRUE) {

    input_dim <- 1
    output_dim <- 1

    gen_A <- function (n) {
        A <- matrix(rnorm(n^2), nrow=n)
        A <- A + t(A)
        eigA <- eigen(A)
        v <- c(0, (-1) * rexp(n-1))
        return(eigA$vectors %*% (v * t(eigA$vectors)))
    }

    n0 <- 4
    A0 <- gen_A(n0)
    B0 <- matrix(rnorm(n0*input_dim), ncol=input_dim)
    C0 <- matrix(rnorm(n0*output_dim), nrow=output_dim)

    n1 <- 6
    A1 <- gen_A(n1)
    B1 <- matrix(rnorm(n1*input_dim), ncol=input_dim)
    C1 <- matrix(rnorm(n1*output_dim), nrow=output_dim)

    system_diff(A0,B0,C0,A1,B1,C1)
    system_diff_slow(A0,B0,C0,A1,B1,C1)

}
