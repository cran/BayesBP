#'Bernstein polynomial basis.
#'@description This function build two dimensional Bernstein polynomial basis.
#'@param ages  Range of ages.
#'@param years Range of years.
#'@param n0 Upper bound of possion random variable.
#'@param N Lower bound of possion random variable.
#'@return Bernstein basis.
#'@examples
#'ages <- 35:85
#'years <- 1988:2007
#'list.basis <- BPbasis(ages,years,10)
#'list.basis
#'@family Bernstein basis
#'@export BPbasis

BPbasis <- function(ages, years, n0, N = 1) {
    x <- scale_to_01(ages)
    y <- scale_to_01(years)
    xy <- expand.grid(x, y)
    basis_list <- list()
    for (n in N:n0) {
        i <- j <- 0:n
        g <- array(0, dim = c(n + 1, n + 1, nrow(xy)))
        for (k in seq_len(nrow(xy))) {
            g[, , k] <- outer(i, j, function(i, j) {
                bin(n, i, xy[k, 1]) * bin(n, j, xy[k, 2])
            })
        }
        parameter <- matrix(as.vector(g), nrow = (n + 1)^2, ncol = nrow(xy))
        basis_list[[length(basis_list) + 1]] <- parameter
    }
    class(basis_list) <- "BPbasis"
    return(basis_list)
}
