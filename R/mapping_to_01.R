#'Mapping to [0,1]
#'@examples
#'mapping_to_01(35:85)
#'(35:85)/(85-35)
#'mapping_to_01(runif(10))
#'@param x Vector.
#'@export mapping_to_01


mapping_to_01 <- function(x) {
    if (length(x) < 2) 
        stop("length(x) must bigger then 2")
    m <- max(x)
    n <- min(x)
    ans <- (x - n)/(m - n)
    return(ans)
}
