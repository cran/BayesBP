#'Credible interval.
#'@description Builing partial differential Bernstein polynomial credible interval.
#'@param result This is output of BP2D.
#'@param n_cluster Muticores is remmended.(default:n_cluster=1)
#'@param alpha     Level of significance.
#'@return This function return Bayesian credible interval with level of significance.
#'@family Credible interval
#'@export Credible_interval_pd_ages
#'
Credible_interval_pd_ages <- function(result, n_cluster = 1, alpha = 0.05) {
    t1 <- Sys.time()
    if (class(result) == "BP2D_result") {
        if (n_cluster > detectCores()) {
            n_cluster <- max(detectCores() - 1, 1)
        }
        ages <- as.numeric(rownames(result$Fhat))
        years <- as.numeric(colnames(result$Fhat))
        cl <- makeCluster(n_cluster)
        clusterExport(cl, c("BPFhat_pd_ages", "BPbasis_pd_ages", "mapping_to_01",
            "ages", "years"))
        a <- parSapply(cl, unlist(result$store_coefficient, recursive = F), function(x) BPFhat_pd_ages(x,
            ages, years))
        stopCluster(cl)
        gc()
        temp <- apply(a, 1, function(x) quantile(x, c(alpha/2, 1 - alpha/2)))
        lowerbound <- matrix(temp[1, ], nrow = length(ages), ncol = length(years))
        upperbound <- matrix(temp[2, ], nrow = length(ages), ncol = length(years))
        row.names(upperbound) <- ages
        colnames(upperbound) <- years
        row.names(lowerbound) <- ages
        colnames(lowerbound) <- years
        r <- list(lowerbound, upperbound)
        names(r) <- paste0(c(alpha/2, 1 - alpha/2) * 100, "%C.I.")
        t2 <- Sys.time()
        message("It costs ", round(difftime(t2, t1, units = "secs"), digits = 3),
            " secs\n")
        return(r)
    } else {
        stop("Input result is not right!")
    }
}
