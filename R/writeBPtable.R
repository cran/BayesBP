#'Write BPtalbe as xlsx file
#'@description If your environment has some result of BP2D,then you can use this function to store BPTable.
#'@import openxlsx
#'@param BPtable output of BP2D_table.
#'@param filename  xlsx file name.
#'@export write.BPtable
write.BPtable <- function(BPtable, filename) {
    wb <- createWorkbook()
    if (class(BPtable) != "BPtable") {
        stop("BPtable is not correct!")
    } else {
        a1 <- BPtable$output
        addWorksheet(wb, "Output")
        writeData(wb, sheet = 1, a1, rowNames = T)
        setColWidths(wb, sheet = 1, cols = 1:10, widths = c(22, rep(13, 9)))
        if ("T." %in% names(BPtable)) {
            a2 <- BPtable$T.
            addWorksheet(wb, "Criterion")
            setColWidths(wb, sheet = 2, cols = 1:3, widths = c(5, 13, 9))
            writeData(wb, sheet = 2, a2, rowNames = T)
        }
        saveWorkbook(wb, filename, overwrite = TRUE)
        return(invisible())
    }
}
