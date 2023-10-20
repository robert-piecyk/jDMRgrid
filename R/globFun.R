# Take 5 digit of decimal value and cut the numbers
floorDec <- function(valParm ,x){
    y <- function(x, level=1) round(x - 5*10^(-level-1), level)
    res <-y(as.numeric(valParm),x)
    return(res)
}

#' @param file_A
#' @importFrom utils
#' @importFrom stringr str_replace_all
#'
statusStringCheck <-  function(file_A){
    list_status <- c("Unmethylated", "Intermediate", "Methylated")
    strTocheckFileA <- head(file_A$status[1])
    if (strTocheckFileA %in% list_status) {
        file_A$status <- str_replace_all(
            file_A$status,pattern = "Unmethylated",replacement = "U")
        file_A$status <- str_replace_all(
            file_A$status,pattern = "Intermediate",replacement = "I")
        file_A$status <- str_replace_all(
            file_A$status,pattern = "Methylated",replacement = "M")
    }
    return(file_A)
}

