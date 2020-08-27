#'@title Table and Criterion.
#'@description  If you give more groups of prior,you can use this function to get the table and T criterion.
#'@param results_list A vector of characters.
#'@return Table and criterion T.
#'@family Bayesain estimate
#'@export BP2D_table

#'
BP2D_table <- function(results_list) {
    if (class(results_list) == "character") {
        l <- length(results_list)
        if (l > 1) {
            allFhat <- lapply(1:l, function(x) get(results_list[x])$Fhat)
            T. <- list()
            for (i in 1:l) {
                for (j in i:l) {
                  if (i != j) {
                    T.[[length(T.) + 1]] <- data.frame(paste0("(", i, ",", j, ")"), 10^5 * 
                      max(abs(allFhat[[i]] - allFhat[[j]])))
                  }
                }
            }
            T. <- do.call(rbind, T.)
            output <- lapply(1:l, function(x) {
                tmp <- cbind(x, get(results_list[x])$output)
                colnames(tmp)[1] <- "set of priors"
                tmp
            })
            output <- do.call(rbind, output)
            colnames(T.) <- c("set of priors", "T")
            rownames(T.) <- seq_len(nrow(T.))
            r <- list(output = output, T. = T.)
            class(r) <- "BPtable"
            return(r)
        } else if (l == 1) {
            r <- get(results_list)$output
            class(r) <- "BPtable"
            return(r)
        }
    } else {
        stop("Please check resultlist again")
    }
}
