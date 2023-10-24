#------------------------------------------------------------------------------
#' Take x digit of decimal value and cut the numbers
#' @param valParm Numeric input. (num)
#' @param x Number of digits to be taken into account. (num)
#' @return Numeric input cut by x digits
floorDec <- function(valParm, x){
    y <- function(x, level=1) round(x - 5*10^(-level-1), level)
    res <-y(as.numeric(valParm),x)
    return(res)
}

#------------------------------------------------------------------------------
#' Replace methylation status by U, I and M states
#' @param file_A DataFrame with state-calls object; must contain status with
#'               'Intermediate', 'Unmethylated' or 'Methylated' calls 
#'               (DataFrame object)
#' @importFrom utils
#' @importFrom stringr str_replace_all
#' @return DataFrame with state-calls object; status calls replaced by U, I or
#'         M state calls
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