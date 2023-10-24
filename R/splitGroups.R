#------------------------------------------------------------------------------
#' Save subset dataset by control/treatment for the next analyses
#' @param input.dir Path to the input directory. (char)
#' @param context Methylation context; 'CG' or 'CHG' or 'CHH'. (char)
#' @param gp1 Filename for group 1. (char)
#' @param gp2 Filename for group 2. (char)
#' @param out.name Group name, as the part of the output filename
#'                 sample_treatment. (char)
#' @param type.name Type name, as the part of the output filename
#'                  (_StateCalls.txt or _rcMethlvl.txt or _postMax.txt). (char)
#' @import magrittr
#' @importFrom dplyr inner_join
#' @importFrom data.table fread fwrite
#' @return Dataset split by the sample/treatment pairs 
#'         saved in the output directory
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
#' @inheritParams saveSplitDataset
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @param postMax.out Logical if postMax matrices should be included. (logical)
#' @param contexts Vector of cytosine contexts; default as c('CG','CHG','CHH').
#'                 (char)
#' @param input.dir Path to the input directory with DMR matrices. (char)
#' @param out.dir Path to the output directory for DMR matrices. (char)
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom dplyr inner_join
#' @return DMR matrices split by sample/treatment groups, saved in the same
#'         directory with DMR matrices.
#' @export
#'
splitGroups <- function(
        samplelist, postMax.out=FALSE, contexts=c("CG","CHG","CHH"), input.dir,
        out.dir)
{
    samplelist$name <- paste0(samplelist$sample,"_", samplelist$replicate)
    gps <- samplelist$group[!samplelist$group %in% c('control')]
    gps <- unique(gps)
    for (m in seq_along(gps)){
        myvec <- c("control", gps[m])
        gp1 <- samplelist$name[which(samplelist$group==myvec[1])]
        gp2 <- samplelist$name[which(samplelist$group==myvec[2])]
        gp1.sample <- unique(samplelist$sample[which(samplelist$name==gp1)])
        gp2.sample <- unique(samplelist$sample[which(samplelist$name==gp2)])
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