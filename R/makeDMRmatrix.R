#------------------------------------------------------------------------------
#' Merge state calls (col6) and rc.meth.lvl  (col7) for all samples into  the 
#' single DataFrame; make a list of these two dataframes
#' @param filepath Paths to the methylome outputs for single samples. (char)
#' @param colm Column indices to be included in the final dataset. (num)
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
prepareFlist <- function(context, extractflist, samplelist) {
    mynames <- gsub(
        paste0("_", context,".txt$"), "", basename(extractflist))
    selectlist <- list()
    message("Extracting filenames and matching them...")
    for (a1 in seq_along(mynames)){
        pat1 <- paste0(
            "_",mynames[a1],"_","|","_",mynames[a1],"|",mynames[a1])
        as <- samplelist[grepl(pat1, samplelist$file),]
        if (NROW(as)==1){
            as$full.path.MethReg <- grep(
                paste0(
                    "/", mynames[a1], "_",context, ".txt", sep=""
                    ), extractflist, value = TRUE)
            message(basename(as$full.path.MethReg)," found!")
            selectlist[[a1]] <- as
        } else {
            message(
                "Multiple files with  match ", mynames[a1]," found!")
            as$full.path.MethReg <- grep(
                paste0(
                    "/", mynames[a1], "_",context, ".txt", sep=""
                ), extractflist, value = TRUE)
            message(basename(as$full.path.MethReg)," found!")
            selectlist[[a1]] <- as
        }
    }
    flist <- rbindlist(selectlist)
    return(flist)
}

#------------------------------------------------------------------------------
#' Builds a DMR matrix for all samples
#' This function generates a binary matrix, a matrix of
#' recalibrated methylation levels and posterior probabilities for all samples.
#' @inheritParams prepareFlist
#' @inheritParams merge_cols
#' @inheritParams writeDMRmatrix
#' @param context Vector of cytosine contexts; default c('CG','CHG','CHH'). 
#'                (char)
#' @param postMax.out Logical if DMR matrix with postMax probabilities
#'                    should be output; default as FALSE. (logical)
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @param input.dir Path to the input directory with methylome calls
#'                  after jDMRgrid function. (char)
#' @param out.dir Path to the output directory. (char)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @import magrittr
#' @importFrom data.table fread fwrite rbindlist
#' @return Saved state-calls, rc-methylation and postMax DMR matrices
#'         as txt files.
#' @export
makeDMRmatrix <- function(
        contexts=c("CG","CHG","CHH"), postMax.out=FALSE, samplelist,
        input.dir, out.dir, include.intermediate=FALSE) 
    {
    for (j in  seq_along(contexts)){
        # list all files in the input directory
        extractflist <- list.files(
            input.dir,pattern=paste0(contexts[j],".txt"),full.names=TRUE)
        if (length(extractflist) != 0){
            # Prepare flist to read files out of it
            flist <- prepareFlist(
                context = contexts[j], extractflist = extractflist, 
                samplelist = samplelist) 
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
