
#' @param filepath
#' @param colm
#' @param include.intermediate
#' @param mincov
#' @param nCytosines
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom dplyr inner_join
#' @export
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

#' Builds a DMR matrix for all samples
#'
#' This function generates a binary matrix, a matrix of
#' recalibrated methylation levels and posterior probabilities for all samples.
#'
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
#'
makeDMRmatrix <- function(
        contexts=c("CG","CHG","CHH"), postMax.out=FALSE, samplefiles,
        input.dir, out.dir, include.intermediate=FALSE) 
    {
    # Read the sample file with filenames
    samplelist <- fread(samplefiles, header=TRUE)
    for (j in  seq_along(contexts)){
        
        # list all files in the input directory
        extractflist <- list.files(
            input.dir,pattern=paste0(contexts[j],".txt"),full.names=TRUE)
        
        #extract file subsets for construction of DMRmatrix
        if (length(extractflist) != 0){
            mynames <- gsub(
                paste0("_", contexts[j],".txt$"), "",basename(extractflist))
            selectlist <- list()
            message("Extracting filenames and matching them...")
            for (a1 in seq_along(mynames)){
                pat1 <- paste0(
                    "_",mynames[a1],"_","|","_",mynames[a1],"|",mynames[a1])
                #pat1 <- paste0("_",mynames[a1],"_")
                as <- samplelist[grepl(pat1, samplelist$file),]
                if (NROW(as)==1){
                    as$full.path.MethReg <- grep(
                        paste0(
                            "/", mynames[a1], "_",contexts[j], ".txt", sep=""
                            ),extractflist, value=TRUE)
                    message(basename(as$full.path.MethReg)," found!")
                    selectlist[[a1]] <- as
                } else {
                    message(
                        "Multiple files with  match ", mynames[a1]," found!")
                }
            }
            flist <- rbindlist(selectlist)
            #print(flist)
            
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
            #(column 6) state-calls and (column 7) rc.meth.lvl
            mydf <- merge_cols(
                flist$full.path.MethReg,include.intermediate,colm=c(5, 6, 7))
            
            # list containing state calls
            status.collect <- mydf[[2]]
            # renaming file names with sample names
            for (a in 4:length(colnames(status.collect))) {
                for (n in seq_along(flist$name)) {
                    if (colnames(status.collect)[a] ==
                        basename(flist$full.path.MethReg)[n]) {
                        colnames(status.collect)[a] <- flist$name[n]
                    }
                }
            }
            # list containing rcmethlvls
            rc.methlevel.collect <- mydf[[3]]
            # renaming file names with sample names
            for (a in 4:length(colnames(rc.methlevel.collect))) {
                for (n in seq_along(flist$name)) {
                    if (colnames(rc.methlevel.collect)[a] ==
                        basename(flist$full.path.MethReg)[n]) {
                        colnames(rc.methlevel.collect)[a] <- flist$name[n]
                    }
                }
            }
            if (postMax.out==TRUE){
                # list containing postmax
                postMax.collect <- mydf[[1]]
                # renaming file names with sample names
                for (a in 4:length(colnames(postMax.collect))) {
                    for (n in seq_along(flist$name)) {
                        if (colnames(postMax.collect)[a] ==
                            basename(flist$full.path.MethReg)[n]) {
                            colnames(postMax.collect)[a] <- flist$name[n]
                        }
                    }
                }
                names(postMax.collect)[1] <- "seqnames"
                names(postMax.collect)[2] <- "start"
                names(postMax.collect)[3] <- "end"
                write.out(
                    postMax.collect,out.dir,"postMax",contexts[j])
            }
            
            names(status.collect)[1] <- "seqnames"
            names(status.collect)[2] <- "start"
            names(status.collect)[3] <- "end"
            
            names(rc.methlevel.collect)[1] <- "seqnames"
            names(rc.methlevel.collect)[2] <- "start"
            names(rc.methlevel.collect)[3] <- "end"
            
            
            write.out(status.collect,out.dir,"StateCalls",contexts[j])
            write.out(rc.methlevel.collect,out.dir,"rcMethlvl",contexts[j])
            
            
            message("Done!")
        } else{
            message("Files for context ", contexts[j],  " do not exist!")
        }
    }
}
