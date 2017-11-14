source("perturb_Kalman.r")

set.seed(23)
num_errors <- 0

check_system <- function (sys) {
    obs_D <- abs(D(sys$A, sys$B, sys$C))
    if (obs_D > 1e-7) {
        cat("Error: D=", obs_D, "\n")
        return(1)
    }
    return(0)
}

cat("check that D(optimal) = 0\n")
optimal_system <- list(
  A = matrix(c(0,-1,1,0), nrow=2),
  B = matrix(c(1,1), ncol=1),
  C = matrix(c(1,0), nrow=1) )
num_errors <- num_errors + check_system(optimal_system)

cat("check perturbation functions\n")
stopifnot(all(make_P(rep(0,25), perturb=TRUE)==0))


for (sigma in c(0.001, 0.1, 1.0)) {
    for (eps in c(0.001, 0.1, 1.0)) {
        cat(sprintf("Testing with sigma=%f and eps=%f\n", sigma, eps))
        AL <- random_oscillator(sigma)
        random_system <- assemble_oscillator(AL)

        cat("check that random_oscillator returns an oscillator\n")
        num_errors <- num_errors + check_system(random_system)

        cat("check that random_oscillator does not change B and C\n")
        num_errors <- num_errors + !(all(abs(random_system$B - c(1,1,rep(0,length(random_system$B)-2))) < 1e-16))
        num_errors <- num_errors + !(all(abs(random_system$C - c(1,0,rep(0,length(random_system$B)-2))) < 1e-16))

        pAL <- perturb_oscillator(AL, eps=eps)
        perturbed_system <- assemble_oscillator(pAL)
        num_errors <- num_errors + !(all(sapply(AL, dim) == sapply(pAL, dim)))

        cat("check that perturb_oscillator returns an oscillator\n")
        num_errors <- num_errors + check_system(perturbed_system)

        cat("check that perturb_oscillator does not change B and C\n")
        num_errors <- num_errors + !(all(abs(perturbed_system$B - c(1,1,rep(0,length(perturbed_system$B)-2))) < 1e-16))
        num_errors <- num_errors + !(all(abs(perturbed_system$C - c(1,0,rep(0,length(perturbed_system$B)-2))) < 1e-16))
    }
}

if (num_errors>0) {
    stop("Failed with ", num_errors, " errors.")
} else {
    cat("All good.\n")
}
