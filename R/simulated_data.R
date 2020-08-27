#'Risky population function
#'@param x Numeric.
#'@param y Numeric.
#'@keywords datasets
#'@export M
M <- function(x, y) {
    r0 = 152040
    r1 = -285270
    r2 = 110410
    r3 = 173900
    r4 = -49950
    r5 = -33630
    r6 = -19530
    r7 = -110330
    r8 = 88840
    r9 = -7990
    population <- r0 + r1 * x + r2 * y +
        r3 * x^2 + r4 * x * y + r5 * y^2 +
        r6 * x^3 +  r7 * x^2 * y +
        r8 * x * y^2 + r9 * y^3
    return(population)
}
NULL

#'Generated data
#'@param ages Ages.
#'@param years Years.
#'@param FT Rate function.
#'@param M  Population function.
#'@keywords datasets
#'@export gen_data
gen_data <- function(ages, years, FT, M) {
    x <- scale_to_01(ages)
    y <- scale_to_01(years)
    disease <- outer(x, y, M) * outer(x, y, FT)
    population <- outer(x, y, M)
    row.names(disease) <- row.names(population) <- ages
    colnames(disease) <- colnames(population) <- years
    return(list(disease = disease, population = population))
}
NULL

#'Generate simulated data 1
#'@description Given rate function 1 generated data.
#'@docType data
#'@keywords datasets
#'@name simulated_data_1
#'@usage data(simulated_data_1)
#'@format list of matrix
#'@examples
#'ages <- 35:85
#'years <- 1988:2007
#'FT1 <- function(x,y){0.00148*sin(0.5*pi*x*y)+0.00002}
#'simulated_data_1 <- gen_data(ages,years,FT1,M)
NULL

#'Generate simulated data 2
#'@description Given rate function 2 generated data.
#'@docType data
#'@keywords datasets
#'@name simulated_data_2
#'@usage data(simulated_data_2)
#'@format list of matrix
#'@examples
#'ages <- 35:85
#'years <- 1988:2007
#'FT2 <- function(x,y){0.00148*sin(0.5*pi*x*(y+0.2))+0.00002}
#'simulated_data_2 <- gen_data(ages,years,FT2,M)
NULL
