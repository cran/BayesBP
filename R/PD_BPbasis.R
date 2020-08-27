#'Partial differential Bernstein polynomial basis.
#'@description This function build two dimensional Bernstein polynomial basis.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@param n0 Upper bound of possion random variable.
#'@param N Lower bound of possion random variable.
#'@param by 1: partial differential by ages; 2: partial differential by years.
#'@return Partial differential Bernstein basis.
#'@examples
#'ages <- 35:85
#'years <- 1988:2007
#'pdbasis <- PD_BPbasis(ages,years,10,by = 1)
#'pdbasis
#'@family Bernstein basis
#'@export PD_BPbasis

PD_BPbasis <- function(ages, years, n0, N = 1, by = 1) {
    x <- scale_to_01(ages)
    y <- scale_to_01(years)
    xy <- expand.grid(x, y)
    basis_list <- list()
    for (n in N:n0) {
        if (by == 1) {
            i <- 0:(n - 1)
            j <- 0:n
            g <- array(0, dim = c(n, n + 1, nrow(xy)))
            for (k in seq_len(nrow(xy))) {
                g[, , k] <- outer(i, j, function(i, j) {
                  n * bin(n - 1, i, xy[k, 1]) * bin(n, j, xy[k, 2])
                })
            }
        } else if (by == 2) {
            i <- 0:n
            j <- 0:(n - 1)
            g <- array(0, dim = c(n + 1, n, nrow(xy)))
            for (k in seq_len(nrow(xy))) {
                g[, , k] <- outer(i, j, function(i, j) {
                  n * bin(n, i, xy[k, 1]) * bin(n, j - 1, xy[k, 2])
                })
            }
        } else {
            stop("Error")
        }
        parameter <- matrix(as.vector(g), nrow = (n + 1) * n, nrow(xy))
        basis_list[[length(basis_list) + 1]] <- parameter
        rm(g, parameter)
    }
    class(basis_list) <- "BPbasis"
    return(basis_list)
}
