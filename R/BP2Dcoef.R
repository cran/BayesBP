#'@title Getting coefficeint from BP2D result.
#'@description This function will return coefficient and length of each set of coefficeint.
#'@param result This is output of BP2D.
#'@return Coefficients table.
#'@family Bayesain estimate
#'@export BP2D_coef


BP2D_coef <- function(result) {
    if (class(result) == "BP2D_result") {
        len_coef <- list()
        for (j in seq_len(length(result$store_coefficient))) {
            len_coef[[j]] <- sapply(result$store_coefficient[[j]], length)
        }
        len_coef <- do.call(cbind, len_coef)
        store_coef <- lapply(result$store_coef, unlist)
        max_len <- max(sapply(store_coef, length))
        store_coef <- sapply(store_coef, function(x) c(x, rep(NA, max_len - length(x))))
        colnames(len_coef) <- colnames(store_coef) <- names(result$store_coefficient)
        row.names(len_coef) <- seq_len(nrow(len_coef))
        row.names(store_coef) <- seq_len(nrow(store_coef))
        return(list(store_coefficient = store_coef, len_coefficient = len_coef))
    } else {
        stop("Please check your input and try again")
    }
}
