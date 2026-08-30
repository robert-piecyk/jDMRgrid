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

#' Sample name derived from a methylome file path
#'
#' \code{runjDMRgrid()} names its per-sample state-call files from the input path rather than
#' from the sheet, and \code{prepareFlist()} has to find those files again. Both now call this
#' helper so the two cannot drift apart.
#'
#' @param file Path(s) to methylome input file(s). (character)
#' @param if.Bismark Whether the inputs are Bismark CX reports. (logical)
#' @return The derived name(s). (character)
#' @keywords internal
deriveSampleName <- function(file, if.Bismark = FALSE) {
    file <- as.character(file)
    if (isTRUE(if.Bismark)) {
        ## FIX: the Bismark branch used gsub("|\\.txt|.CX_report.txt$", "", x). The pattern
        ## begins with an EMPTY alternative, which matches the empty string at every
        ## position and therefore strips nothing -- the "name" stayed a full path.
        sub("\\.CX_report\\.txt(\\.gz)?$", "", basename(file))
    } else {
        ## Accept gzipped input. data.table::fread() reads .gz transparently and the toy
        ## data ships compressed, but the pattern below only strips ".txt", so
        ## "methimpute_p1.txt.gz" became the sample name "methimpute_p1.gz".
        file <- sub("\\.gz$", "", file)
        basename(gsub(".*methylome_|\\.txt|_All.txt$", "", file))
    }
}

#-----------------------------------------------------------------------------
#' Coerce a GRanges to a data.frame without S3/S4 dispatch
#'
#' Bioconductor devel makes \code{GRanges} extend \code{List}, so
#' \code{as.data.frame()} on a GRanges selects \code{S4Vectors:::as.data.frame.List},
#' which calls \code{unlist()} on the object and fails with "GRanges objects don't
#' support [[, as.list(), lapply(), or unlist() at the moment". This builds the
#' frame from accessors, reproducing the column layout, names and factor levels of
#' \code{as.data.frame()} on every Bioconductor version.
#'
#' @param x A GRanges object. (GRanges)
#' @return A data.frame with columns seqnames, start, end, width, strand followed
#'         by the metadata columns; matrix metadata columns expand to
#'         \code{name.subname} exactly as \code{as.data.frame()} expands them.
#' @keywords internal
grToDF <- function(x) {
    df <- data.frame(
        seqnames = factor(as.character(GenomeInfoDb::seqnames(x)),
                          levels = GenomeInfoDb::seqlevels(x)),
        start  = BiocGenerics::start(x),
        end    = BiocGenerics::end(x),
        width  = BiocGenerics::width(x),
        strand = factor(as.character(BiocGenerics::strand(x)),
                        levels = c("+", "-", "*")),
        stringsAsFactors = FALSE)
    mc <- S4Vectors::mcols(x, use.names = FALSE)
    if (!is.null(mc) && NCOL(mc) > 0L) {
        cols <- lapply(seq_len(NCOL(mc)), function(i) mc[[i]])
        names(cols) <- colnames(mc)
        df <- cbind(df, do.call(data.frame,
                                c(cols, list(stringsAsFactors = FALSE))))
    }
    ## as.data.frame() on a GRanges leaves explicit integer row names, while data.frame()
    ## leaves them compact. The frames compare identical() either way, but cbind() warns
    ## "row names were found from a short variable and have been discarded" against a
    ## compact frame, so match as.data.frame() exactly.
    attr(df, "row.names") <- seq_len(NROW(df))
    df
}
