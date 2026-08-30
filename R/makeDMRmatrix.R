#------------------------------------------------------------------------------
#' Merge state calls (col6) and rc.meth.lvl  (col7) for all samples into  the 
#' single DataFrame; make a list of these two dataframes
#' @param filepath Paths to the methylome outputs for single samples. (char)
#' @param colm Column indices to be included in the final dataset. (num)
#' @param if.Bismark Whether the inputs were Bismark CX reports; controls how sample
#'   names are derived from file paths. (logical)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom dplyr inner_join
#' @return List of two data frames (for state calls and rc.meth.lvl) with
#'         merged columns
#' 
merge_cols <- function(filepath, colm, include.intermediate) {
    mylist <- list()
    for (l in seq_along(colm)){
        extract <- lapply(filepath, function(k){
            f <- fread(k,header=FALSE,skip=1,select=c(1, 2, 3, colm[l]))
            if (colm[l]==6) {
                if (include.intermediate==TRUE) {
                    f[,4] <- ifelse(f[,4] == "U", yes = 0, (
                        ifelse(f[,4] == "I", yes = 0.5, no = 1)))
                } else {
                    f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
                }
            }
            colnames(f)[4] <- basename(k)
            return(f)
            
        })
        df <- Reduce(function(x, y) {inner_join(x, y, by=c(
            "V1","V2","V3"))}, extract)
        
        mylist[[l]] <- df
    }
    return(mylist)
}

#------------------------------------------------------------------------------
#' Write DMR matrix data as txt file
#' @param out.df Output DataFrame with the dataset to be saved as txt file.
#'               (DataFrame object)
#' @param data.dir Path to the output directory. (char)
#' @param out.name Name of the file to be saved. (char)
#' @param contexts Methylation context; 'CG' or 'CHG' or 'CHH' . (char)
#' @importFrom data.table fwrite
#' @return DMR matrix saved in the output directory for the given context and
#'         data type given in the out.name
#' 
write.out <- function(out.df, data.dir, out.name, contexts) {
    fwrite(x=out.df, file=paste0(
        data.dir, "/", contexts,"_", out.name, ".txt"), quote=FALSE, 
        row.names=FALSE, col.names=TRUE, sep="\t")
}

#------------------------------------------------------------------------------
#' Builds a combined matrix for a given dataset (s-c, rcmthlvl or postMax) 
#' @inheritParams write.out
#' @param data DataFrame with the dataset to be saved as txt file.
#'             (DataFrame object)
#' @param data.type Type of the data; "postMax" or "StateCalls" or "rcMethlvl"
#'                  (char)
#' @param flist List of the full.path.MethReg to be saved as the columns; 
#'              gained from the prepareFlist function. (DataFrame object)
#' @param out.dir Path to the output directory. (char)
#' @param context Methylation context "CG" or "CHG" or "CHH". (char)
#' @return Saved DMR matrix in output directory for a given context and type
#' 
writeDMRmatrix <- function(data, data.type, flist, out.dir, context) {
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    for (a in 4:ncol(data)) {
        for (n in seq_along(flist$name)) {
            if (colnames(data)[a] ==  basename(flist$full.path.MethReg)[n]) {
                colnames(data)[a] <- flist$name[n]
            }
        }
    }
    
    names(data)[1] <- "seqnames"
    names(data)[2] <- "start"
    names(data)[3] <- "end"
    
    write.out(data, out.dir, data.type, context)
}

#------------------------------------------------------------------------------
#' Extract filenames and match them to the methylation contexts
#' @param context Methylation context; "CG" or "CHG" or "CHH". (char)
#' @param extractflist Vector with the files 
#'                     listed in the input directory. (char)
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @importFrom data.table rbindlist
#' @return Path to the methylome input files matching a given context
#' 
prepareFlist <- function(context, extractflist, samplelist,
                         if.Bismark = FALSE) {
    samplelist <- data.table::as.data.table(samplelist)
    if (is.null(samplelist$sample))
        stop("samplelist needs a 'sample' column")
    if (is.null(samplelist$file))
        stop("samplelist needs a 'file' column")

    ## FIX: this used to loop over the FILES and regex-grep each derived name against
    ## samplelist$file. That never consulted the 'sample' column the sheet advertises, the
    ## pattern contained unescaped regex metacharacters (a derived name such as
    ## "102F_All.gz" has a "." that matches anything), a zero-match silently produced an
    ## empty row that surfaced later as "subscript out of bounds", and the multiple-match
    ## branch did exactly what the single-match branch did.
    ##
    ## Now the sample sheet is authoritative: for each row we derive the name runjDMRgrid()
    ## would have written, match it exactly, and label the column with samplelist$sample.
    message("Matching state-call files to the samples in the sheet...")
    want  <- paste0(deriveSampleName(samplelist$file, if.Bismark), "_", context, ".txt")
    idx   <- match(want, basename(extractflist))
    gone  <- is.na(idx)
    if (any(gone))
        stop("no ", context, " state-call file for sample(s): ",
             paste(samplelist$sample[gone], collapse = ", "),
             ". Expected: ", paste(want[gone], collapse = ", "),
             ". runjDMRgrid() names its outputs from the methylome file path, so the sheet's ",
             "'file' column must point at the same inputs that produced these calls.")
    dup <- duplicated(idx)
    if (any(dup))
        stop("several samples map to the same state-call file: ",
             paste(samplelist$sample[dup], collapse = ", "),
             ". Sample names must be distinct in the sheet.")

    flist <- data.table::copy(samplelist)
    flist$full.path.MethReg <- extractflist[idx]
    message("Matched ", nrow(flist), " sample(s) for ", context, ".")
    flist
}

#------------------------------------------------------------------------------
#' Builds a DMR matrix for all samples
#' This function generates a binary matrix, a matrix of
#' recalibrated methylation levels and posterior probabilities for all samples.
#' @inheritParams prepareFlist
#' @inheritParams merge_cols
#' @inheritParams writeDMRmatrix
#' @param input.dir Path to the directory holding the region-call files written
#'                  by runjDMRgrid(). (char)
#' @param contexts Vector of cytosine contexts; default c('CG','CHG','CHH').
#'                 (char)
#' @param postMax.out Logical: also write the matrix of posterior probabilities
#'                     behind each state call. Needed by filterDMRmatrix(min.posterior=)
#'                     and by rankDMRs(). Costs roughly 35-40 percent more output.
#'                     (logical; default as TRUE)
#' @examples
#' sheet <- get(load(system.file("data", "listFiles2.RData", package = "jDMRgrid")))
#' sheet$file <- system.file("extdata", sheet$file, package = "jDMRgrid")
#'
#' grid <- file.path(tempdir(), "jDMRgrid_grid2")
#' mat  <- file.path(tempdir(), "jDMRgrid_matrix2")
#' runjDMRgrid(out.dir = grid, window = 200, step = 50, samplelist = sheet,
#'             contexts = "CG", min.C = 10, min.C.type = "percentile",
#'             mincov = 0, include.intermediate = FALSE, runName = "example")
#' makeDMRmatrix(contexts = "CG", samplelist = sheet,
#'               input.dir = grid, out.dir = mat, include.intermediate = FALSE)
#' list.files(mat)
makeDMRmatrix <- function(
        contexts=c("CG","CHG","CHH"), postMax.out=TRUE, samplelist,
        input.dir, out.dir, include.intermediate=FALSE, if.Bismark=FALSE) 
    {
    ## FIX: with an empty or wrong input.dir this returned quietly having written nothing,
    ## so a broken pipeline looked like a successful run.
    if (!dir.exists(input.dir)) stop("input.dir does not exist: ", input.dir)
    for (cx in contexts) {
        if (length(list.files(input.dir, pattern = paste0(cx, "\\.txt$"))) == 0L)
            stop("no per-sample region-call files for context ", cx, " in ", input.dir,
                 " -- did runjDMRgrid() write state calls? Check that paths in ",
                 "`samplefiles` resolve.")
    }
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    for (j in  seq_along(contexts)){
        # list all files in the input directory
        extractflist <- list.files(
            input.dir,pattern=paste0(contexts[j],".txt"),full.names=TRUE)
        if (length(extractflist) != 0){
            # Prepare flist to read files out of it
            flist <- prepareFlist(
                context = contexts[j], extractflist = extractflist, 
                samplelist = samplelist, if.Bismark = if.Bismark) 
            # Assign unique names for samples with or without replicate data
            if (!is.null(flist$replicate)) {
                message(
                    "Running context ", contexts[j], " for replicates...")
                flist$name <- paste0(flist$sample,"_", flist$replicate)
            } else {
                flist$name <- flist$sample
            }
            message("Now, constructing DMR matrix for ", contexts[j])
            # merge samples by Chr coordinates
            mydf <- merge_cols(
                flist$full.path.MethReg,include.intermediate,colm=c(5, 6, 7))
            # list containing state calls
            writeDMRmatrix(
                data = mydf[[2]], data.type = "StateCalls", flist = flist, 
                out.dir = out.dir, context = contexts[j])
            # list containing rcmethlvls
            writeDMRmatrix(
                data = mydf[[3]], data.type = "rcMethlvl", flist = flist, 
                out.dir = out.dir, context = contexts[j])
            # list containing postMax
            if (postMax.out == TRUE) {
                writeDMRmatrix(
                    data = mydf[[1]], data.type = "postMax", flist = flist, 
                    out.dir = out.dir, context = contexts[j])}
            message("Done!")
        } else{
            message("Files for context ", contexts[j],  " do not exist!")
        }
    }
}
