#'Bernstein polynomial basis.
#'@description This function build two dimensional Bernstein polynomial basis.
#'@param n0 Upper bound of possion random variable.
#'@param nn Lower bound of possion random variable.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@return Bernstein basis.
#'@examples
#'ages<-35:85
#'years<-1988:2007
#'list.basis<-BPbasis(10,ages,years)
#'list.basis
#'@family Bernstein basis
#'@export BPbasis



BPbasis <- function(n0, ages, years, nn = 1) {
    x <- mapping_to_01(ages)
    y <- mapping_to_01(years)
    lx <- length(x)
    ly <- length(y)
    x1 <- rep(x, ly)
    y1 <- rep(y, each = lx)
    list.basis <- list()
    bin <- function(n, i, x) {
        choose(n, i) * x^i * (1 - x)^(n - i)
    }
    basis <- function(i, j) {
        bin(n, i, x1[k]) * bin(n, j, y1[k])
    }
    for (n in nn:n0) {
        i <- 0:n
        j <- 0:n
        g <- array(0, dim = c(n + 1, n + 1, lx * ly))
        for (k in 1:(lx * ly)) {
            g[, , k] <- outer(i, j, basis)
        }
        parameter <- matrix(as.vector(g), nrow = (n + 1)^2, lx * ly)
        list.basis[[length(list.basis) + 1]] <- parameter
        rm(g, parameter)
    }
    return(list.basis)
}
