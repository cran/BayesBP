#'Scale to [0,1]
#'@examples
#'scale_to_01(35:85)
#'(35:85-35)/(85-35)
#'scale_to_01(runif(10))
#'@param x Vector.
#'@export scale_to_01
scale_to_01 <- function(x) {
    a <- min(x)
    b <- max(x)
    return((x - a)/(b - a))
}
NULL
#'Binomial function
#'@param n Integer.
#'@param i Integer(i < n).
#'@param x Numeric(0<= x <=1).
#'@examples
#'bin(5,3,.5)
#'@export bin
bin <- function(n, i, x) {
    choose(n, i) * x^i * (1 - x)^(n - i)
}
NULL
