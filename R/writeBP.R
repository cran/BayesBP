#'Write xlsx file
#'@description This function will write result of BP2D to xlsx file.
#'@import openxlsx
#'@param writedata result of BP2D(character or list).
#'@param filename  xlsx file name.
#'@export write.BP
write.BP <- function(writedata, filename) {
    wb <- createWorkbook()
    if (is.character(writedata)) {
        writedata <- get(writedata)
    } else if (class(writedata) != "BP2D_result") {
        stop("writedata is not correct!")
    }
    a1 <- writedata$Fhat
    temp <- BP2D_coef(writedata)
    a2 <- temp$store_coefficient
    a3 <- temp$len_coefficient
    rm(temp)
    a4 <- writedata$chain
    a5 <- writedata$maxchain
    a6 <- writedata$output
    addWorksheet(wb, "Fhat")
    writeData(wb, sheet = 1, a1, rowNames = T)
    addWorksheet(wb, "Store_coefficient")
    writeData(wb, sheet = 2, a2)
    addWorksheet(wb, "len_coefficient")
    writeData(wb, sheet = 3, a3)
    addWorksheet(wb, "chain")
    writeData(wb, sheet = 4, a4)
    addWorksheet(wb, "maxchain")
    writeData(wb, sheet = 5, a5)
    addWorksheet(wb, "output")
    writeData(wb, sheet = 6, a6, rowNames = T)
    setColWidths(wb, sheet = 6, cols = 1:10, widths = c(19, rep(12, 9)))
    saveWorkbook(wb, filename, overwrite = TRUE)
    return(invisible())
}
