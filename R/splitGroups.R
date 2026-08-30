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
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

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
#' @examples
#' ## a small DMR matrix, written to a temporary directory
#' d <- file.path(tempdir(), "jDMRgrid_example")
#' dir.create(d, showWarnings = FALSE)
#' co  <- data.frame(seqnames = 1L, start = c(1L, 201L, 401L), end = c(200L, 400L, 600L))
#' smp <- c("WT_rep1", "WT_rep2", "mutant1_rep1", "mutant1_rep2")
#' st  <- rbind(c(1, 1, 0, 0), c(0, 0, 1, 1), c(1, 1, 1, 1))
#' mk  <- function(nm, M) {
#'     x <- data.frame(co, M)
#'     names(x) <- c("seqnames", "start", "end", smp)
#'     data.table::fwrite(x, file.path(d, nm), sep = "\t")
#' }
#' mk("CG_StateCalls.txt", st)
#' mk("CG_rcMethlvl.txt", st * 0.9)
#' mk("CG_postMax.txt", matrix(0.99, 3, 4))
#'
#' sheet <- data.frame(
#'     file = NA_character_, sample = rep(c("WT", "mutant1"), each = 2),
#'     replicate = rep(c("rep1", "rep2"), 2),
#'     group = rep(c("control", "treatment1"), each = 2), stringsAsFactors = FALSE)
#'
#' ## one matrix per pairwise comparison
#' splitGroups(samplelist = sheet, contexts = "CG", input.dir = d, out.dir = d)
#' list.files(d, pattern = "WT_mutant1")
#' @export
#'
splitGroups <- function(
        samplelist, postMax.out=TRUE, contexts=c("CG","CHG","CHH"), input.dir,
        out.dir)
{
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

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