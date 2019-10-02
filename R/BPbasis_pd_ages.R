#'Partial differential Bernstein polynomial basis.
#'@description This function build two dimensional partial differential Bernstein polynomial basis by ages.
#'@param n0 Upper bound of possion random variable.
#'@param nn Lower bound of possion random variable.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@return Partial differential Bernstein basis by ages.
#'@examples
#'ages<-35:85
#'years<-1988:2007
#'list.basis<-BPbasis_pd_ages(10,ages,years)
#'list.basis
#'@family Bernstein basis
#'@export BPbasis_pd_ages



BPbasis_pd_ages <- function(n0, ages, years, nn = 1) {
    x <- mapping_to_01(ages)
    y <- mapping_to_01(years)
    lx <- length(x)
    ly <- length(y)
    x1 <- rep(x, length(y))
    y1 <- rep(y, each = length(x))
    list.basis <- list()
    bin <- function(n, i, x) {
        choose(n, i) * x^i * (1 - x)^(n - i)
    }
    basis <- function(i, j) {
        n * bin(n - 1, i, x1[k]) * bin(n, j, y1[k])
    }
    for (n in nn:n0) {
        i <- 0:(n - 1)
        j <- 0:n
        g <- array(0, dim = c(n, n + 1, lx * ly))
        for (k in 1:(lx * ly)) {
            g[, , k] <- outer(i, j, basis)
        }
        parameter <- matrix(as.vector(g), nrow = (n + 1) * n, lx * ly)
        list.basis[[length(list.basis) + 1]] <- parameter
        rm(g, parameter)
    }
    return(list.basis)
}

