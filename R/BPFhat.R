#'Two dimensional Bernstein polynomial
#'@description Given Bernstein polynomial coeffients to compute Fhat.
#'@param coef  Bernstein polynomial coefficients.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@param basis Bernstein polynomial basis.
#'@return This function return outer Bernstein polynomial using coefficients.
#'@examples
#'coef <- runif(9)
#'ages <- 35:85
#'years <- 1988:2007
#'list.basis <- BPbasis(ages,years,10)
#'BPFhat(coef,ages,years,list.basis)
#'@family outer Bernstein polynomial
#'@export BPFhat

BPFhat <- function(coef, ages, years, basis) {
    coef <- as.vector(coef)
    lx <- length(ages)
    ly <- length(years)
    n <- sqrt(length(coef))
    if (round(n) != n) {
        stop("length(coef) should be square number")
    }
    d <- sapply(basis, dim)
    chs <- which(d[1, ] == (n * n))
    est <- colSums(coef * basis[[chs]])
    return(matrix(est, lx, ly))
}


