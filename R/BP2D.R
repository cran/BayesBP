#'@title Bayesian estimation using two dimensions Bernstein polynomial
#'@description This function runs Metropolis-Hasting algorithm which is given setting prior and data.This algorithm starts storing coefficients when it runs halfway,so we use second halves of coefficients compute Rhat to check convergence.
#'@import parallel iterators stats utils
#'@references Li-Chu Chien,Yuh-Jenn Wu,Chao A. Hsiung,Lu-Hai Wang,I-Shou Chang(2015).Smoothed Lexis Diagrams With Applications to Lung and Breast Cancer Trends in Taiwan,Journal of the American Statistical Association, Taylor & Francis Journals, vol. 110(511), pages 1000-1012, September.
#'@param prior prior=(n0,alpha,L) where alpha is a Poisson parameter,n0 is upper bound of alpha
#'L can be every number which is bigger than one.
#'@param input_data  It contain disease and population(ex: simulated_data_1).
#'@param ages        Range of ages.
#'@param years       Range of years.
#'@param Iterations  Iterations of chain.
#'@param n_chain     Number of Markov chain.
#'@param RJC         Control parameter for transfer dimension.
#'@param nn          The parameter nn is lower bound of alpha.
#'@param seed        Set seed yes or not.
#'@param set         Choose seed.(defaults:set=1)
#'@param n_cluster   This parameter means number of cores, five cores is recommended.(default: n_cluster=1).
#'@param interval    Each hundreds save one coefficient.
#'@param double      If R.hat >1.1 then double the iterations of times.
#'@return
#'  This function will return Bayesian estimate of incidence,Stored parameters,posterior mean,posterior max and table.
#'\item{Fhat}{Bayesian estimate of incidence.}
#'\item{chain}{Bayesian estimate of posterior p-value mean.}
#'\item{maxchain}{Bayesian estimate of posterior p-value max.}
#'\item{store_coefficients}{Two dimensional Bernstein coefficients.}
#'\item{output}{When M-H algorithm ends,contruct the table which contains norm,mean of Fhat,maximum of Fhat,R.hat,iterations,P-value and elasped time.}
#'
#'@examples
#'\donttest{
#'library(BayesBP)
#'#simulated_data_1,simulated_data_2
#'#Ages 1~85,years 1988~2007
#'#Data are zero from 0 to 34
#'#Given one prior and simulated_data_1
#'data('simulated_data_1')
#'ages<-35:85
#'years<-1988:2007
#'prior<-c(10,5,2)
#'result<-BP2D(prior,simulated_data_1,ages,years,n_cluster=1)
#'result$Fhat
#'result$chain
#'result$maxchain
#'result$output
#'result$store_coefficients$chain_1
#'matplot(result$chain,type='l',main='Posterior mean trace plot')
#'matplot(result$maxchain,type='l',main='Posterior max trace plot')
#'BP2D_coef(result)
#'write.BP(result,filename = 'result.xlsx')
#'write.BP('result',filename = 'result.xlsx')
#'
#'#Given four prior and simulated_data_2
#'data('simulated_data_2')
#'n0<-c(10,20,10,20)
#'alpha<-c(5,10,5,10)
#'L<-c(2,2,4,4)
#'prior<-cbind(n0,alpha,L)
#'ages<-35:85
#'years<-1988:2007
#'results_list<-paste0('result_',letters[1:4])
#'for(i in 1:4){
#'  assign(results_list[i],BP2D(prior[i,],simulated_data_2,ages,years,n_cluster=1))
#'}
#'BPtable<-BP2D_table(results_list)
#'write.BPtable(BPtable,filename = 'BPtable.xlsx')
#'mapply(write.BP,results_list,paste0(results_list,'.xlsx'))
#'#Credible interval
#'CI<-Credible_interval(result,n_cluster = 1)
#'CI_pda<-Credible_interval_pd_ages(result,n_cluster = 1)
#'CI_pdy<-Credible_interval_pd_years(result,n_cluster = 1)
#'CI
#'CI_pda
#'CI_pdy
#'}
#'@export BP2D
#'@family Bayesain estimate



BP2D <- function(prior, input_data = input_data, ages = ages, years = years, Iterations = 2e+05,
    n_cluster = 1, n_chain = 5, RJC = 0.35, nn = 2, seed = TRUE, set = 1, interval = 100,
    double = 4) {
    if ((prior[1] < prior[2]) || (length(prior) != 3)) {
        stop("The prior is not correct!")
    } else if (RJC > 0.5 || RJC < 0) {
        stop("'RJC' is not correct!")
    } else if (nn < 1 && (!round(nn) == nn)) {
        stop("Can not compute in this method!")
    } else if (n_cluster > detectCores()) {
        n_cluster <- max(detectCores() - 1, 1)
    } else if (!round(double) == double) {
        stop("'double' should be integer!")
    } else if (n_chain %in% c(0, 1)) {
        stop("'n_chain' should bigger then one")
    }
    n0 <- prior[1]
    alpha <- prior[2]
    LL <- max(prior[3], 1)
    tao11 <- min(ages)
    tao12 <- max(ages)
    tao21 <- min(years)
    tao22 <- max(years)
    input_data <- as.matrix(input_data)
    nc_d <- ncol(input_data)/2
    Ndata <- input_data[, 1:nc_d]
    Mdata <- input_data[, 1:nc_d + nc_d]
    disease <- Ndata[tao11:tao12, tao21:tao22 - tao21 + 1]
    population <- Mdata[tao11:tao12, tao21:tao22 - tao21 + 1]
    rdm <- disease/population
    if ((sum(population == 0) > 0) || (!dim(Ndata) == dim(Mdata))) {
        stop("Please check your data again")
    }
    disease <- as.matrix(disease)
    population <- as.matrix(population)
    x <- (ages - tao11)/(tao12 - tao11)
    y <- (years - tao21)/(tao22 - tao21)
    lx <- length(x)
    ly <- length(y)
    list.basis <- BPbasis(n0, ages, years)
    Bernstein <- function(a) {
        n <- sqrt(length(a))
        aa <- matrix(a, nrow = n^2, ncol = lx * ly)
        temp <- colSums(aa * list.basis[[n - 1]])
        Mean <- mean(temp)
        Max <- max(temp)
        return(list(Max = Max, Mean = Mean))
    }
    Pvalue <- function(a) {
        n <- sqrt(length(a))
        aa <- matrix(a, nrow = n^2, ncol = lx * ly)
        M <- matrix(colSums(aa * list.basis[[n - 1]]), lx, ly)
        M1 <- matrix(sapply(M * population, function(x) rpois(1, x)), lx, ly)/population
        EZR1 <- max(abs(M1 - M))
        EZR2 <- abs(mean(M1) - mean(M))
        ZR1 <- max(abs(M - rdm))
        ZR2 <- abs(mean(M) - mean(rdm))
        PMax <- EZR1 > ZR1
        PMean <- EZR2 > ZR2
        return(list(Max = PMax, Mean = PMean))
    }
    Fhat <- function(a) {
        n <- sqrt(length(a))
        aa <- matrix(a, nrow = n^2, ncol = lx * ly)
        matrix(colSums(aa * list.basis[[n - 1]]), lx, ly)
    }
    logLF <- function(a) {
        n <- sqrt(length(a))
        aa <- matrix(a, nrow = n^2, ncol = lx * ly)
        M <- colSums(aa * list.basis[[n - 1]]) * population
        sum(disease * log(M) - M - lfactorial(disease))
    }

    r1 <- LL * max(rdm[1, ])
    r2 <- LL * max(rdm[nrow(rdm), ])
    c1 <- LL * max(rdm[, 1])
    c2 <- LL * max(rdm[, ncol(rdm)])
    ld <- c(r1, r2, c1, c2)
    M2 <- LL * max(rdm)
    M1 <- 0

    if (seed == TRUE) {
        set.seed(set)
    } else {
        rm(set)
    }
    P1 <- sum(dpois(0:nn, alpha))
    Pn <- sapply((nn + 1):(n0 - 1), function(x) dpois(x, alpha))
    P <- c(P1, Pn, 1 - sum(Pn) - P1)
    n1 <- sample(nn:n0, n_chain, P, replace = T)
    P <- c(rep(0, nn - 1), P, 0)
    initial <- function(n1) {
        h1 <- c(runif(n1 + 1, M1, r1))
        h2 <- matrix(c(runif(n1 - 1, M1, c1), runif((n1 - 1)^2, M1, M2), runif(n1 -
            1, M1, c2)), nrow = n1 - 1, ncol = n1 + 1)
        h3 <- c(runif(n1 + 1, M1, r2))
        Rprm <- rbind(h1, h2, h3)
        return(Rprm)
    }
    initialvalue <- lapply(n1, initial)
    MHRJ_Agorithm <- function(initialvalue, kk = 0.5) {
        gg <- kk * Iterations
        store <- list()
        storeF <- 0
        chain <- c()
        maxchain <- c()
        Pmax <- c()
        Pmean <- c()
        n1 <- sqrt(length(unlist(initialvalue))) - 1
        Rprm <- matrix(unlist(initialvalue), n1 + 1, n1 + 1)
        R <- as.vector(Rprm)
        si <- (ifelse(kk == 0.5, 1, kk * Iterations + 1))
        ei <- (2 * kk * Iterations)
        for (i in si:ei) {
            RJP <- c(min(1, P[n1 - 1]/P[n1]), min(1, P[n1 + 1]/P[n1]))
            RJP <- c(RJC * RJP[1], 1 - sum(RJC * RJP), RJC * RJP[2])
            if (RJP[1] == 0) {
                H <- sample(c(1, 2), 1, prob = RJP[2:3])
            } else if (RJP[3] == 0) {
                H <- sample(c(0, 1), 1, prob = RJP[1:2])
            } else {
                H <- sample(c(0, 1, 2), 1, prob = RJP)
            }
            Sprm <- Rprm
            if (H == 0) {
                n2 <- n1 - 1
                delete <- sample(0:n1 + 1, 2, replace = T)
                Sprm <- Rprm[-delete[1], -delete[2]]
                S <- as.vector(Sprm)
                Jcb <- sum(log(ld - M1)) - log(M2 - M1) * (2 * n2 - 3)
                rio <- logLF(S) - logLF(R) + log(P[n2]/P[n1]) - Jcb
            } else if (H == 2) {
                n2 <- n1 + 1
                Sprm <- Rprm
                add <- sample(0:n1, 2, replace = T)
                m1 <- Rprm[1:add[1], ]
                m2 <- Rprm[add[1]:nrow(Rprm), ]
                m <- rbind(m1, m2)
                m1 <- m[, 1:add[2]]
                m2 <- m[, add[2]:ncol(m)]
                Sprm <- cbind(m1, m2)
                ls <- sqrt(length(Sprm)) - 2
                ww <- c(runif(1, M1, r1), runif(ls, M1, M2), runif(1, M1, r2))
                vv <- c(runif(1, M1, c1), runif(ls, M1, M2), runif(1, M1, c2))
                Sprm[add[1] + 1, ] <- ww
                Sprm[, add[2] + 1] <- vv
                S <- as.vector(Sprm)
                Jcb <- sum(log(ld - M1)) + log(M2 - M1) * (2 * n2 - 1)
                rio <- logLF(S) - logLF(R) - log(P[n2]/P[n1]) + Jcb
            } else {
                Sprm <- Rprm
                n2 <- n1
                if (n2 == 1) {
                  stay <- sample(1:4, 1)
                } else {
                  stay <- sample(1:9, 1)
                }
                if (stay == 1) {
                  Sprm[1, 1] <- runif(1, M1, r1)
                } else if (stay == 2) {
                  Sprm[1, ncol(Sprm)] <- runif(1, M1, r1)
                } else if (stay == 3) {
                  Sprm[nrow(Sprm), 1] <- runif(1, M1, r2)
                } else if (stay == 4) {
                  Sprm[nrow(Sprm), ncol(Sprm)] <- runif(1, M1, r2)
                } else if (stay == 5) {
                  Sprm[1, 2:(ncol(Sprm) - 1)] <- runif(ncol(Sprm) - 2, M1, r1)
                } else if (stay == 6) {
                  Sprm[nrow(Sprm), 2:(ncol(Sprm) - 1)] <- runif(ncol(Sprm) - 2,
                    M1, r2)
                } else if (stay == 7) {
                  Sprm[2:(nrow(Sprm) - 1), 1] <- runif(ncol(Sprm) - 2, M1, c1)
                } else if (stay == 8) {
                  Sprm[2:(nrow(Sprm) - 1), ncol(Sprm)] <- runif(ncol(Sprm) - 2,
                    M1, c2)
                } else {
                  ch <- sample(2:n1, 1)
                  Sprm[ch, 2:(nrow(Sprm) - 1)] <- runif(ncol(Sprm) - 2, M1, M2)
                }
                S <- as.vector(Sprm)
                rio <- logLF(S) - logLF(R)
            }
            if (rio > 0) {
                nnext <- n2
                Xnext <- Sprm
            } else {
                if (log(runif(1)) < rio) {
                  nnext <- n2
                  Xnext <- Sprm
                } else {
                  nnext <- n1
                  Xnext <- Rprm
                }
            }
            Rprm <- Xnext
            R <- as.vector(Xnext)
            n1 <- nnext
            condition <- ifelse(interval == 1, i > gg, i > gg && (i%%gg%%interval ==
                1))
            if (condition) {
                chain[length(chain) + 1] <- Bernstein(Xnext)$Mean
                maxchain[length(maxchain) + 1] <- Bernstein(Xnext)$Max
                tmp <- Pvalue(Xnext)
                Pmax[length(Pmax) + 1] <- tmp$Max
                Pmean[length(Pmean) + 1] <- tmp$Mean
                store[[length(store) + 1]] <- Xnext
                storeF <- storeF + Fhat(Xnext)
                rm(tmp)
            }
        }
        return(list(chain = chain, store = store, storeF = storeF/gg * interval,
            maxchain = maxchain, Pmax = Pmax, Pmean = Pmean))
    }
    norm_fun <- function(x, y) {
        c(L1norm = mean(abs(x - y)), supnorm = max(abs(x - y)))
    }
    kk <- 0.5
    out.put <- list()
    cl <- makeCluster(n_cluster)
    t1 <- Sys.time()
    prior <- paste0("n0=", prior[1], ",alpha=", prior[2], ",max=", prior[3])
    message("Prior: ", prior)
    repeat {
        result <- parLapply(cl, 1:n_chain, function(i) {
            if (seed == T) {
                set.seed(set)
                a <- runif(n_chain)
                set.seed(a[i])
            }
            MHRJ_Agorithm(initialvalue[[i]], kk = kk)
        })
        gc()
        chain <- rbind(sapply(1:n_chain, function(i) result[[i]]$chain))
        maxchain <- rbind(sapply(1:n_chain, function(i) result[[i]]$maxchain))
        R.hat <- Rhat(chain, burn.in = 0)
        Fprime <- 0
        for (j in 1:n_chain) {
            Fprime <- Fprime + (result[[j]]$storeF)/n_chain
        }
        Px <- rbind(sapply(1:n_chain, function(i) result[[i]]$Pmax))
        Pn <- rbind(sapply(1:n_chain, function(i) result[[i]]$Pmean))
        Fhat.mean <- mean(Fprime) * 10^5
        Fhat.max <- max(Fprime) * 10^5
        norm <- norm_fun(rdm, Fprime) * 10^5
        Pvalue.mean <- mean(Pn)
        Pvalue.max <- mean(Px)
        t2 <- Sys.time()
        out.put[[length(out.put) + 1]] <- cbind(L1_norm = norm[1], sup_norm = norm[2],
            Fhat.mean, Fhat.max, R.hat, Iterations = Iterations * kk * 2, Pvalue.mean,
            Pvalue.max, time.mins = difftime(t2, t1, units = "mins"))
        message("Iterations: ", Iterations * kk * 2)
        message("Time: ", round(difftime(t2, t1, units = "mins"), digits = 3),
            " mins")
        message("R.hat: ", round(R.hat, digits = 6))
        if (R.hat < 1.1) {
            message("Markov chains is convergence.\n")
            break
        } else if (kk == 2^double/2) {
            message("Markov chains is not convergence.\n")
            break
        } else {
            message("Markov chains is not convergence.")
            message("Double the number of iterations.")
            kk <- kk * 2
            initialvalue <- lapply(1:n_chain, function(i) tail(result[[i]]$store,
                1))
        }
    }
    stopCluster(cl)
    storeparameter <- lapply(1:n_chain, function(i) result[[i]]$store)
    names(storeparameter) <- paste0("chain_", 1:n_chain)
    names(storeparameter) <- colnames(chain) <- colnames(maxchain) <- paste0("chain_",
        1:n_chain)
    output <- do.call(rbind, out.put)
    rownames(output) <- rep(prior, dim(output)[1])
    row.names(Fprime) <- ages
    colnames(Fprime) <- years
    r <- list(Fhat = Fprime, chain = chain, maxchain = maxchain, store_coefficients = storeparameter,
        output = output)
    class(r) <- "BP2D_result"
    return(r)
}
