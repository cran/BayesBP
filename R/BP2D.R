#'@title Bayesian estimation using two dimensions Bernstein polynomial
#'@description This function runs Metropolis-Hasting algorithm which is given setting prior and data.This algorithm starts storing coefficients when it runs halfway,so we use second halves of coefficients compute Rhat to check convergence.
#'@import parallel iterators stats utils
#'@references Li-Chu Chien,Yuh-Jenn Wu,Chao A. Hsiung,Lu-Hai Wang,I-Shou Chang(2015).Smoothed Lexis Diagrams With Applications to Lung and Breast Cancer Trends in Taiwan,Journal of the American Statistical Association, Taylor & Francis Journals, vol. 110(511), pages 1000-1012, September.
#'@param prior prior=(n0,alpha,L) where alpha is a Poisson parameter,n0 is upper bound of alpha
#'L can be every number which is bigger than one.
#'@param disease     Disease matrix.
#'@param population  Population matrix.
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
#'# ---------------------------------------- #
#'library(BayesBP)
#'ages<-35:85
#'years<-1988:2007
#'prior<-c(10,5,2)
#'data(simulated_data_1)
#'disease<-simulated_data_1$disease
#'population<-simulated_data_1$population
#'result<-BP2D(prior,ages,years,disease,population)
#'# ---------------------------------------- #
#'# Bernstein basis
#'basis<-BPbasis(ages,years,10)
#'pdbasis1<-PD_BPbasis(ages,years,10,by = 1)
#'pdbasis2<-PD_BPbasis(ages,years,10,by = 2)
#'# Bernstein polynomial
#'coef<-result$store_coefficients$chain_1[[1]]
#'BPFhat(coef,ages,years,basis)
#'PD_BPFhat(coef,ages,years,pdbasis1,by = 1)
#'PD_BPFhat(coef,ages,years,pdbasis2,by = 2)
#'# Credible interval
#'Credible_interval(result)
#'PD_Credible_interval(result,by = 1)
#'PD_Credible_interval(result,by = 2)
#'# ---------------------------------------- #
#'# Given four prior set
#'ages<-35:85
#'years<-1988:2007
#'data(simulated_data_2)
#'disease<-simulated_data_2$disease
#'population<-simulated_data_2$population
#'p<-expand.grid(n0=c(10,20),alpha=c(5,10),LL=c(2,4))
#'prior_set<-p[p$n0==p$alpha*2,]
#'result_list<-paste0('result',1:nrow(prior_set))
#'for (i in seq_len(nrow(prior_set))) {
#'    prior<-prior_set[i,]
#'    assign(result_list[i],BP2D(prior,ages,years,disease,population))
#'    write.BP(get(result_list[i]),sprintf('%s.xlsx',result_list[i]))
#'}
#'tab<-BP2D_table(result_list)
#'write.BPtable(tab,'result_table.xlsx')
#'# ---------------------------------------- #
#'}
#'@export BP2D
#'@family Bayesain estimate

BP2D <- function(prior, ages, years, disease, population, Iterations = 2e+05, n_chain = 5,
    n_cluster = 1, nn = 2, interval = 100, RJC = 0.35, seed = TRUE, set = 1, double = 4) {
    if ((prior[1] < prior[2]) || (length(prior) != 3)) {
        stop("'prior' is not correct!")
    } else if (RJC > 0.5 || RJC < 0) {
        stop("'RJC' is not correct!")
    } else if (nn < 1 && (!round(nn) == nn)) {
        stop("'nn' is not correct!")
    } else if (n_cluster > detectCores()) {
        n_cluster <- max(detectCores() - 1, 1)
    } else if (!round(double) == double) {
        stop("'double' should be integer!")
    } else if (n_chain %in% c(0, 1)) {
        stop("'n_chain' should bigger then one")
    }

    prior <- unlist(prior)
    n0 <- prior[1]
    alpha <- prior[2]
    LL <- max(prior[3], 1)
    basis <- BPbasis(ages, years, n0)
    disease <- as.matrix(disease)
    population <- as.matrix(population)
    rate <- disease/population

    lx <- length(ages)
    ly <- length(years)
    l1 <- LL * max(rate[1, ], na.rm = T)
    l2 <- LL * max(rate[nrow(rate), ], na.rm = T)
    l3 <- LL * max(rate[, 1], na.rm = T)
    l4 <- LL * max(rate[, ncol(rate)], na.rm = T)
    ld <- c(l1, l2, l3, l4)
    M2 <- LL * max(rate, na.rm = T)
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

    initialvalue <- lapply(n1, function(n) {
        h1 <- runif(n + 1, M1, l1)
        h2 <- matrix(c(runif(n - 1, M1, l3), runif((n - 1)^2, M1, M2), runif(n - 1, M1,
            l4)), nrow = n - 1, ncol = n + 1)
        h3 <- runif(n + 1, M1, l2)
        ret <- rbind(h1, h2, h3)
        return(ret)
    })

    BP <- function(a) {
        a <- as.vector(a)
        n <- sqrt(length(a))
        est <- colSums(a * basis[[n - 1]])
        return(matrix(est, lx, ly))
    }

    logLF <- function(a) {
        M <- BP(a) * population
        sum(disease * log(M) - M, na.rm = T)
    }

    MHRJ_Algorithm <- function(initialvalue, kk = 0.5) {
        iter <- kk * Iterations
        store_coef <- list()
        storeF <- 0
        chain <- c()
        maxchain <- c()
        Pmax <- c()
        Pmean <- c()
        nx <- nrow(initialvalue) - 1
        X <- initialvalue
        start <- (ifelse(kk == 0.5, 1, kk * Iterations + 1))
        end <- (2 * kk * Iterations)
        for (i in start:end) {
            RJP <- c(min(1, P[nx - 1]/P[nx]), min(1, P[nx + 1]/P[nx]))
            RJP <- c(RJC * RJP[1], 1 - sum(RJC * RJP), RJC * RJP[2])
            if (RJP[1] == 0) {
                H <- sample(c(1, 2), 1, prob = RJP[2:3])
            } else if (RJP[3] == 0) {
                H <- sample(c(0, 1), 1, prob = RJP[1:2])
            } else {
                H <- sample(c(0, 1, 2), 1, prob = RJP)
            }
            Y <- X
            if (H == 0) {
                ny <- nx - 1
                delete <- sample(0:nx + 1, 2, replace = T)
                Y <- X[-delete[1], -delete[2]]
                Jcb <- sum(log(ld - M1)) - log(M2 - M1) * (2 * ny - 3)
                rio <- logLF(Y) - logLF(X) + log(P[ny]/P[nx]) - Jcb
            } else if (H == 2) {
                ny <- nx + 1
                Y <- matrix(0, ny + 1, ny + 1)
                s <- sample(2:ny, 2, replace = T)
                Y[-s[1], -s[2]] <- X
                Y[s[1], ] <- c(runif(1, M1, l1), runif(nx, M1, M2), runif(1, M1, l2))
                Y[, s[2]] <- c(runif(1, M1, l3), runif(nx, M1, M2), runif(1, M1, l4))
                Jcb <- sum(log(ld - M1)) + log(M2 - M1) * (2 * ny - 1)
                rio <- logLF(Y) - logLF(X) - log(P[ny]/P[nx]) + Jcb
            } else {
                Y <- X
                ny <- nx
                n <- nrow(Y)
                if (ny == 1) {
                  state <- sample(1:4, 1)
                } else {
                  state <- sample(1:9, 1)
                }
                if (state == 1) {
                  Y[1, 1] <- runif(1, M1, l1)
                } else if (state == 2) {
                  Y[1, n] <- runif(1, M1, l1)
                } else if (state == 3) {
                  Y[n, 1] <- runif(1, M1, l2)
                } else if (state == 4) {
                  Y[n, n] <- runif(1, M1, l2)
                } else if (state == 5) {
                  Y[1, 2:(n - 1)] <- runif(n - 2, M1, l1)
                } else if (state == 6) {
                  Y[n, 2:(n - 1)] <- runif(n - 2, M1, l2)
                } else if (state == 7) {
                  Y[2:(n - 1), 1] <- runif(n - 2, M1, l3)
                } else if (state == 8) {
                  Y[2:(n - 1), n] <- runif(n - 2, M1, l4)
                } else {
                  ch <- sample(2:ny, 1)
                  Y[ch, 2:(n - 1)] <- runif(n - 2, M1, M2)
                }
                rio <- logLF(Y) - logLF(X)
            }
            if (rio > 0) {
                nnext <- ny
                Xnext <- Y
            } else {
                if (log(runif(1)) < rio) {
                  nnext <- ny
                  Xnext <- Y
                } else {
                  nnext <- nx
                  Xnext <- X
                }
            }
            X <- Xnext
            nx <- nnext
            condition <- ifelse(interval == 1, i > iter,
                                i > iter && (i%%iter%%interval == 1))
            if (condition) {
                est <- BP(Xnext)
                chain[length(chain) + 1] <- mean(est)
                maxchain[length(maxchain) + 1] <- max(est)
                Rp <- rpois(lx * ly, est * population)/population
                EZR1 <- max(abs(est - Rp), na.rm = T)
                EZR2 <- abs(mean(est, na.rm = T) - mean(Rp, na.rm = T))
                ZR1 <- max(abs(est - rate), na.rm = T)
                ZR2 <- abs(mean(est, na.rm = T) - mean(rate, na.rm = T))
                Pmax[length(Pmax) + 1] <- EZR1 > ZR1
                Pmean[length(Pmean) + 1] <- EZR2 > ZR2
                store_coef[[length(store_coef) + 1]] <- Xnext
                storeF <- storeF + est
            }
        }
        return(list(chain = chain, maxchain = maxchain, store_coef = store_coef, storeF = storeF/iter *
            interval, Pmax = Pmax, Pmean = Pmean))
    }

    kk <- 0.5
    out.put <- list()
    cl <- makeCluster(n_cluster)
    t1 <- Sys.time()
    pr <- paste0("n0=", prior[1], ",alpha=", prior[2], ",max=", prior[3])
    message("Prior: ", pr)
    repeat {
        result <- parLapply(cl, 1:n_chain, function(i) {
            if (seed == T) {
                set.seed(set)
                a <- sample(-10000:10000,n_chain)
                set.seed(a[i])
            }
            MHRJ_Algorithm(initialvalue[[i]], kk = kk)
        })
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
        L1norm <- mean(abs(rate - Fprime), na.rm = T) * 10^5
        supnorm <- max(abs(rate - Fprime), na.rm = T) * 10^5
        Pvalue.mean <- mean(Pn)
        Pvalue.max <- mean(Px)
        t2 <- Sys.time()
        out.put[[length(out.put) + 1]] <- cbind(L1_norm = L1norm, sup_norm = supnorm, Fhat.mean,
            Fhat.max, R.hat, Iterations = Iterations * kk * 2, Pvalue.mean, Pvalue.max,
            time.mins = difftime(t2, t1, units = "mins"))
        message("Iterations: ", Iterations * kk * 2)
        message("Time: ", round(difftime(t2, t1, units = "mins"), digits = 3), " mins")
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
            initialvalue <- sapply(1:n_chain, function(i)
                tail(result[[i]]$store_coef, 1))
        }
    }
    stopCluster(cl)
    storeparameter <- lapply(1:n_chain, function(i) result[[i]]$store_coef)
    names(storeparameter) <- colnames(chain) <- colnames(maxchain) <-
        paste0("chain_", 1:n_chain)
    output <- do.call(rbind, out.put)
    rownames(output) <- rep(pr, dim(output)[1])
    row.names(Fprime) <- ages
    colnames(Fprime) <- years
    r <- list(Fhat = Fprime, chain = chain, maxchain = maxchain,
              store_coefficients = storeparameter, output = output)
    class(r) <- "BP2D_result"
    return(r)
}

