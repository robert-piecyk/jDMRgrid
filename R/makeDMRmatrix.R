#'
#' @param filepath
#' @param colm
#' @param include.intermediate
#' @param mincov
#' @param nCytosines
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom dplyr inner_join
#'

# This function will merge (column 6) state calls and (column 7)
# rc.meth.lvl from all samples into one dataframe
# makes list of 2 dataframes
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

write.out <- function(out.df, data.dir, out.name, contexts){
    fwrite(x=out.df, file=paste0(
        data.dir,"/",contexts,"_",out.name,".txt"
        ),quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

#' Builds a combined matrix for a given dataset (s-c, rcmthlvl or postMax) 
#' @param data Dataset
#' @param data.type Character (postMax, StateCalls or rcMethlvl)
#' @param flist List of the output matrices for a given context
#' @param out.dir Output directory to save dataset
#' @param context Methylation context (CG, CHG or CHH)
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

#' Prepare flist out of names 
#' @param context Methylation context (CG, CHG or CHH)
#' @param extractflist All files in the input directory
#' @param samplelist Sample files in the input directory
#' @importFrom data.table rbindlist
prepareFlist <- function(context, extractflist, samplelist) {
    mynames <- gsub(
        paste0("_", context,".txt$"), "",basename(extractflist))
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
        }
    }
    flist <- rbindlist(selectlist)
    return(flist)
}

#' Builds a DMR matrix for all samples
#' This function generates a binary matrix, a matrix of
#' recalibrated methylation levels and posterior probabilities for all samples.
#' @param context Vector of cytosine contexts. (character vector)
#' @param samplefiles Path to the text file. (character)
#' @param include.intermediate A logical if intermediate calls should be used.
#' @param input.dir Path to the methylome calls. (character)
#' @param out.dir Path to output directory. (character)
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table rbindlist
#' @export
makeDMRmatrix <- function(
        contexts=c("CG","CHG","CHH"), postMax.out=FALSE, samplefiles,
        input.dir, out.dir, include.intermediate=FALSE) 
    {
    samplelist <- fread(samplefiles, header=TRUE) # Read the sample file 
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
            writeDMRmatrix(
                data = mydf[[3]], data.type = "rcMethlvl", flist = flist, 
                out.dir = out.dir, context = contexts[j])
            message("Done!")
        } else{
            message("Files for context ", contexts[j],  " do not exist!")
        }
    }
}
