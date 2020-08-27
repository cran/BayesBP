#'Credible interval.
#'@description Builing two dimensional Bernstein polynomial credible interval.
#'@param result This is output of BP2D.
#'@param n_cluster Muticores is remmended.(default:n_cluster=1)
#'@param alpha     Level of significance.
#'@param by 1: partial differential by ages; 2: partial differential by years.
#'@references L.H. Chien, T.J. Tseng, C.H. Chen, H.F. Jiang, F.Y. Tsai, T.W. Liu, C.A. Hsiung, I.S. Chang Comparison of annual percentage change in breast cancer incidence rate between Taiwan and the United States-A smoothed Lexis diagram approach.
#'@return Bayesian credible interval with level of significance.
#'@family Credible interval
#'@export PD_Credible_interval
#'
PD_Credible_interval <- function(result, n_cluster = 1, alpha = 0.05, by = 1) {
    if (class(result) != "BP2D_result") {
        stop("'result' is not correct!")
    } else {
        coef_list <- unlist(result$store_coefficient, recursive = F)
        ages <- as.numeric(rownames(result$Fhat))
        years <- as.numeric(colnames(result$Fhat))
        lx <- length(ages)
        ly <- length(years)
        n0 <- max(sqrt(unique(sapply(coef_list, length))))
        pdbasis <- PD_BPbasis(ages, years, n0)
        LU <- c(alpha/2, 1 - alpha/2)
        cl <- makeCluster(n_cluster)
        clusterExport(cl, c("PD_BPFhat", "scale_to_01"))
        ret <- parSapply(cl, coef_list, function(x) PD_BPFhat(x, ages, years, pdbasis, 
            by))
        ret <- t(parApply(cl, ret, 1, function(x) c(quantile(x, LU), mean(x))))
        stopCluster(cl)
        lowerCI <- matrix(ret[, 1], lx, ly)
        upperCI <- matrix(ret[, 2], lx, ly)
        M <- matrix(ret[, 3], lx, ly)
        row.names(lowerCI) <- row.names(M) <- row.names(upperCI) <- ages
        colnames(lowerCI) <- colnames(M) <- colnames(upperCI) <- years
        ret <- list(lowerCI, AVG = M, upperCI)
        names(ret)[c(1, 3)] <- paste0(LU * 100, "%CI")
        return(ret)
    }
}


