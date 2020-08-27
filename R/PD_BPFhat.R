#'Two dimensional Bernstein polynomial
#'@description Given Bernstein polynomial coeffients to compute Fhat.
#'@param coef  Bernstein polynomial coefficients.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@param pdbasis Partial differential Bernstein polynomial basis.
#'@param by 1: partial differential by ages; 2: partial differential by years.
#'@return Partial differential Bernstein polynomial given coefficients.
#'@examples
#'coef <- runif(9)
#'ages <- 35:85
#'years <- 1988:2007
#'pdbasis <- PD_BPbasis(ages,years,10,N=1,by=1)
#'PD_BPFhat(coef,ages,years,pdbasis,by=1)
#'@family outer Bernstein polynomial
#'@export PD_BPFhat
PD_BPFhat <- function(coef, ages, years, pdbasis, by = 1) {
    coef <- as.vector(coef)
    lx <- length(ages)
    ly <- length(years)
    n <- sqrt(length(coef))
    if (round(n) != n) {
        stop("length(coef) should be square number")
    }
    coef <- matrix(coef, n, n)
    d <- sapply(pdbasis, dim)
    chs <- which(d[1, ] == (n * n - n))
    if (by == 1) {
        coef <- coef[2:n, ] - coef[1:(n - 1), ]
        coef <- as.vector(coef)
        est <- colSums(coef * pdbasis[[chs]])
    } else if (by == 2) {
        coef <- coef[, 2:n] - coef[, 1:(n - 1)]
        coef <- as.vector(coef)
        est <- colSums(coef * pdbasis[[chs]])
    }
    return(matrix(est, lx, ly))
}
