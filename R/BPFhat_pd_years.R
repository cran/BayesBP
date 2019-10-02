#'Bernstein polynomial
#'@description Given Bernstein polynomial coeffients to compute partial differential by years Fhat.
#'@param coef  Bernstein polynimial coefficients.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@return This function return outer Bernstein polynomial using coefficients.
#'@examples
#'coef<-runif(9)
#'ages<-35:85
#'years<-1988:2007
#'BPFhat_pd_years(coef,ages,years)
#'@family outer Bernstein polynomial
#'@export BPFhat_pd_years

BPFhat_pd_years <- function(coef, ages, years) {
    n <- sqrt(length(coef))
    if (round(n) != n) {
        stop("length(coef) should be square number")
    }
    a.basis <- BPbasis_pd_years(n - 1, ages, years, n - 1)
    coef <- matrix(coef, n, n)
    coef <- coef[, 2:n] - coef[, 1:(n - 1)]
    aa <- matrix(coef, nrow = n * (n - 1), ncol = length(ages) * length(years))
    m <- matrix(colSums(aa * a.basis[[1]]), length(ages), length(years))
    colnames(m) <- years
    row.names(m) <- ages
    return(m)
}
