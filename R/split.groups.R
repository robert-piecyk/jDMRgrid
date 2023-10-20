#------------------------------------------------------------------------------
#' Save subset dataset by control/treatment for the next analyses
#' @param input.dir
#' @param context
#' @param gp1
#' @param gp2
#' @param out.name
#' @param type.name
#' @import magrittr
#' @importFrom dplyr inner_join
#' @importFrom data.table fread fwrite
#'
saveSplitDataset <- function(
        input.dir, context, gp1, gp2, out.dir, out.name, type.name)
{
    fname <- paste0(input.dir, '/', context, type.name)
    if (file.exists(fname)) {
        data <- fread(fname)
        df1 <- subset(
            data,, which(colnames(data) %in% c(
                'seqnames','start','end', gp1, gp2)))
        fwrite(x=df1, file=paste0(
            out.dir, "/", context, "_", out.name, type.name), quote=FALSE, 
            row.names=FALSE, col.names=TRUE, sep="\t")
    }
}

#------------------------------------------------------------------------------
#' Split into control/treatment groups according to the sample and replicates
#' @param samplefiles Path to the text file. (character)
#' @param postMax.out A logical if postMax matrices should be included.
#' @param contexts Vector of cytosine contexts. (character vector)
#' @param input.dir Input directory with DMR matrices files. (character)
#' @param out.dir Output directory for DMR matrices files. (character)
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom dplyr inner_join
#' @export
#'
split.groups <- function(
        samplefiles,postMax.out=FALSE,contexts=c("CG","CHG","CHH"),input.dir,
        out.dir)
    {
    ft <- fread(samplefiles)
    ft$name <- paste0(ft$sample,"_", ft$replicate)
    gps <- ft$group[!ft$group %in% c('control')]
    gps <- unique(gps)
    for (m in seq_along(gps)){
        myvec <- c("control", gps[m])
        gp1 <- ft$name[which(ft$group==myvec[1])]
        gp2 <- ft$name[which(ft$group==myvec[2])]
        gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
        gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
        out.name <- paste0(gp1.sample, "_", gp2.sample)
        for (cn in seq_along(contexts)){
            saveSplitDataset(
                input.dir = input.dir, context = contexts[cn], gp1 = gp1, 
                gp2 = gp2, out.name = out.name, out.dir = out.dir,
                type.name = "_StateCalls.txt")
            saveSplitDataset(
                input.dir = input.dir, context = contexts[cn], gp1 = gp1, 
                gp2 = gp2, out.dir = out.dir, out.name = out.name, 
                type.name = "_rcMethlvl.txt")
            if (postMax.out==TRUE) {
                saveSplitDataset(
                    input.dir = input.dir, context = contexts[cn], gp1 = gp1, 
                    gp2 = gp2, out.dir = out.dir, out.name = out.name, 
                    type.name = "_postMax.txt")
            }
        }
    }
}