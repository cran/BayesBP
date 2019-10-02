#'Gelman Rubin statistics.
#'@description Check Markov chains for convergence.
#'@param M An n x m numeric matrix of Markov Chains.
#'@param burn.in The default value 0.5 means that the second halves of chains will be used to compute.
#'@return Gelman Rubin statistics.
#'@references Gelman A.,Carlin J.B.,Stern H.S.,and Rubin D.B.(2004),Bayesian Data Analysis,Boca Raton,FL:Chapman&Hall/CRC.
#'@export Rhat
#'
Rhat <- function(M, burn.in = 0.5) {
    m <- ncol(M)
    x <- M[round(((burn.in * nrow(M)) + 1), 0):nrow(M), ]
    n <- nrow(x)
    phibar <- mean(x)
    phi <- colMeans(x)
    B <- n/(m - 1) * sum((phi - phibar)^2)
    W <- mean(apply(x, 2, var))
    post.var <- (n - 1)/n * W + B/n
    r <- sqrt(post.var/W)
    r <- ifelse(r < 1, 1, r)
    return(r)
}
