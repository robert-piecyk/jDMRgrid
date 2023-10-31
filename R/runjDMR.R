#-----------------------------------------------------------------------------
#' Export bins from the list and save it as RData file using dput
#' @param mylist List of two DataFrame objects : binned genome and statistics 
#'               meeting minC criterion. (list of DataFrame objects)
#' @param out.dir Path to the output directory. (char)
#' @param runName Character defining the name of the run; 
#'                default as GridGenome. (char)
#' @return Binned genome saved in output directory
#'
export.bins <- function(mylist, out.dir, runName)
{
    names_mylist <- names(mylist)
    lapply(seq_along(names_mylist), function(z1) {
        names_mylist.1 <- names_mylist[z1]
        window.size.1 <- strsplit(names_mylist.1, '_')[[1]][[2]]
        step.size.1 <- strsplit(names_mylist.1, '_')[[1]][[3]]
        mylist.1 <- mylist[[z1]]
        context.1 <- strsplit(names_mylist.1, '_')[[1]][[1]]
        out.name.1 <- paste0(
            out.dir, "/",runName,"_Win", window.size.1,"_Step", step.size.1,"_",
            context.1,".Rdata",sep="")
        dput(mylist.1, out.name.1)
    })
}

#-----------------------------------------------------------------------------
#' Main function to bin genome using non/sliding window approach
#' @param x1 Index of the set of window/step size to be checked; in single-mode
#'           it is set to 1. (num)
#' @param window Bin size in non/sliding window approach. (num)
#' @param step Step size in non/sliding window approach. (num)
#' @param gr GenomicRanges object of the start-end positions within genome.
#'           (GenomicRanges object)
#' @param cyt_gr GenomicRanges object of cytosines for all contexts.
#'               (GenomicRanges object)
#' @param contexts Vector of cytosine contexts; 
#'                 default as c('CG','CHG','CHH'). (char)
#' @param min.C Percentile threshold based on EDF. (num; between 0 and 100)
#' @importFrom data.table fread
#' @importFrom stats quantile ecdf
#' @importFrom IRanges countOverlaps
#' @return List of two DataFrame objects with binned genome and statistics 
#'         how many regions is meeting minC criterion.
#'
binGenomeLoop <- function(x1, window, step, gr, cyt_gr, contexts, min.C) {
    window.size <- window[x1]
    step.size <- step[x1]
    # Binning genome
    binned.g <- slidingWindows(gr, window.size, step.size)
    message("Binning genome with windows of: ", window.size,
            "bp and step-size of: ", step.size, "bp.")
    # Creating a data frame from the binned data
    dd <- data.frame(unlist(binned.g))
    names(dd) <- c("chr", "start", "end", "cluster.length", "strand")
    # Storing the data frame in a list
    new <- list(dd)
    names(new) <- as.numeric(as.character(window.size))
    # Creating a GRange object
    data_gr <- GRanges(
        seqnames = dd$chr, ranges = IRanges(
            start = dd$start, end = dd$end), 
        clusterlen = dd$cluster.length)
    # Creating a data frame for bin and step sizes
    mydf <- data.frame(bin.size = window.size, step.size = step.size)
    # Process each context using lapply and combine the results
    results <- lapply(contexts, function(cx) {
        message("Extracting cytosines for ", cx, ".")
        # Filtering cytosines
        ref_gr <- cyt_gr[which(cyt_gr$context == cx),]
        # Counting cytosines in GRanges
        dat.collect <- countOverlaps(data_gr, ref_gr)
        # Create a empirical distribution of cytosines within bins and
        # find a threshold based on its min.C percentile
        #new.dat.collect <- dat.collect[which(dat.collect >= min.C)]
        new.dat.collect <- dat.collect[which(dat.collect >= as.numeric(
            quantile(ecdf(dat.collect),min.C/100)))]
        non.empty.bins <- length(new.dat.collect) / length(dat.collect)
        # Create a filtrated data frame
        data.out <- dd[which(dat.collect >= as.numeric(quantile(
            ecdf(dat.collect),min.C/100))),]
        return(list(non.empty.bins, data.out))
    })
    new.one <- lapply(results, function(x) x[[2]])
    names(new.one) <- paste0(contexts, '_', window.size, '_', step.size)
    # Add the results to the mydf data frame
    mydf <- cbind(mydf,data.frame(context = contexts, ratio = unlist(
        lapply(results, function(x) x[[1]]))))
    return(list(mydf = mydf, collect.bins = new.one))
}

#-----------------------------------------------------------------------------
#' Bin genome using non/sliding window approach
#' @inheritParams binGenomeLoop
#' @inheritParams export.bins
#' @param methimputefiles Vector of paths to text files of the sample reads.
#'                        (char)
#' @param contexts Vector of cytosine contexts; 
#'                 default as c('CG','CHG','CHH'). (char)
#' @param window Bin size in non/sliding window approach. (num)
#' @param step Step size in non/sliding window approach. (num)
#' @param min.C Percentile threshold based on EDF. (num; between 0 and 100)
#' @param out.dir Path to the output directory. (char)
#' @param runName Character defining the name of the run; 
#'                default as GridGenome. (char)
#' @param if.Bismark Logical if Bismark inputs (CX reports in txt format)
#'                   are used. Default as FALSE. (logical)
#' @param FASTA.file Path to the FASTA file; required if Bismark outputs are
#'                   used. Default as NULL. (char)
#' @importFrom data.table fread rbindlist fwrite
#' @importFrom stats quantile ecdf
#' @importFrom IRanges countOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings readDNAStringSet
#' @return Binned genome saved in the output directory using non/sliding
#'         window approach.
#'
binGenome <- function(
        methimputefiles, contexts, window, step, min.C, out.dir, runName,
        if.Bismark, FASTA.file)
{
    # message about creating grid
    message('Creating grid...')
    cyt.collect <- list()
    # from one of the methIMPUTE file extract all cytosines positions
    meth.out <- fread(methimputefiles[1],showProgress=FALSE)
    all.cyt.pos <- meth.out[,c('seqnames','start','strand','context')]
    colnames(all.cyt.pos) <- c('chr','pos','strand','context')
    # create a GRanges object from the cytosines positions
    cyt_gr <- GRanges(seqnames=all.cyt.pos$chr,ranges=IRanges(
        start=all.cyt.pos$pos, width=1),context=all.cyt.pos$context,
        strand=all.cyt.pos$strand)
    # get length of chromosomesś
    if (if.Bismark == FALSE & FASTA.file == NULL) {
        chr.names <- unique(meth.out$seqnames)
        chr.lengths <- unlist(
            lapply(chr.names, function(x) max(
                meth.out$start[meth.out$seqnames == x])))
        names(chr.lengths) <- chr.names
    } else {
        fasta.out <- readDNAStringSet(FASTA.file)
        chr.lengths <- unlist(lapply(fasta.out, length))
    }
    # create a GRanges object from the chromosome (fasta) start-end positions
    gr <- GRanges(seqnames=names(chr.lengths),ranges=IRanges(
        start=1, end=chr.lengths))
    # main loop function
    if (length(window) != length(step)) {
        message('Window and step vectors sizes must have same length.')
    } else {
        results <- lapply(
            seq_along(window), binGenomeLoop, window = window, step = step, 
            gr = gr, cyt_gr = cyt_gr, contexts = contexts, min.C = min.C)
        out <- rbindlist(lapply(results, function(x) x$mydf))
        out <- out[order(out$context, out$bin.size, out$step.size),]
        fwrite(out,file = paste0(
            out.dir, '/', runName, '_optimal_minC_threshold.csv'))
        mybins <- split(out, f=out$context)
        #mybins <- lapply(out, function(x) x[which.min(x$ratio),])
        collect.bins <- lapply(results, function(x) x$collect.bins)
        message("Exporting regions...")
        lapply(collect.bins, function(x) export.bins(
            mylist=x, out.dir=out.dir,runName=runName))
        return(list.files(out.dir, pattern=paste0(
                ".*", runName, ".*\\.Rdata$"), full.names=TRUE))
    }
    message("Done!")
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using future_apply parallel method
#' @inheritParams makeMethimpute
#' @param out.samplelist DataFrame consisting of: contexts, methImpute files,
#'                     sample names (from methImpute files) and IDs 
#'                     corresponding to the binned genome. (DataFrame object)
#' @param merge_list List of DataFrames with binned regions; for each context
#'                   separately (list of DataFrame objects)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @param out.dir Path to the output directory. (char)
#' @param mincov Minimum read coverage; default as 0. (num; between 0 and 1)
#' @param if.Bismark Logical if Bismark inputs (CX reports in txt format)
#'                   are used. Default as FALSE. (logical)
#' @param FASTA.file Path to the FASTA file; required if Bismark outputs are
#'                   used. Default as NULL. (char)
#' @import magrittr
#' @import future
#' @import future.apply
#' @import doParallel
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#' 
makeMethimpute_future <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov,
        if.Bismark, FASTA.file)
{
    plan(multisession)
    info_lapply <- future_lapply(
        seq_along(out.samplelist$context), function(j) {
            refRegion <- list(reg.obs = merge_list[[out.samplelist$id[j]]])
            message("Running file: ",out.samplelist$methfn[j]," for context: ",
                    out.samplelist$context[j], "\n")
            grid.out <- makeMethimpute(
                df = as.character(out.samplelist$file[j]),
                context = out.samplelist$context[j],
                refRegion = refRegion, fit.plot = FALSE,
                include.intermediate = include.intermediate,
                probability = "constrained", out.dir = out.dir,
                fit.name = paste0(
                    basename(out.samplelist$methfn[j]), "_",
                    out.samplelist$context[j]),
                name = basename(out.samplelist$methfn[j]), mincov = mincov,
                if.Bismark = if.Bismark, FASTA.file = FASTA.file)}, 
        future.seed = NULL)
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using future_apply parallel method
#' @inheritParams makeMethimpute
#' @param out.samplelist DataFrame consisting of: contexts, methImpute files,
#'                     sample names (from methImpute files) and IDs 
#'                     corresponding to the binned genome. (DataFrame object)
#' @param merge_list List of DataFrames with binned regions; for each context
#'                   separately (list of DataFrame objects)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @param out.dir Path to the output directory. (char)
#' @param mincov Minimum read coverage; default as 0. (num; between 0 and 1)
#' @param numCores Number of cores to be used if for_each method should be 
#'                 performed; default as NULL. (num)
#' @param if.Bismark Logical if Bismark inputs (CX reports in txt format)
#'                   are used. Default as FALSE. (logical)
#' @param FASTA.file Path to the FASTA file; required if Bismark outputs are
#'                   used. Default as NULL. (char)
#' @import magrittr
#' @import foreach
#' @import doParallel
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster stopCluster                
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#'
makeMethimpute_foreach <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov, 
        numCores, if.Bismark, FASTA.file)
{
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    runMethimputeJ <- function(jk) {
        refRegion <- list(reg.obs = merge_list[[out.samplelist$id[jk]]])
        message("Running file: ", out.samplelist$methfn[jk],
                " for context: ", out.samplelist$context[jk], "\n")
        grid.out <- makeMethimpute(
            df = as.character(out.samplelist$file[jk]),
            context = out.samplelist$context[jk],
            refRegion = refRegion, fit.plot = FALSE,
            include.intermediate = include.intermediate,
            probability = "constrained",out.dir = out.dir,
            fit.name = paste0(
                basename(out.samplelist$methfn[jk]), "_",
                out.samplelist$context[jk]),
            name = basename(out.samplelist$methfn[jk]), mincov = mincov,
            if.Bismark = if.Bismark, FASTA.file = FASTA.file)
        return(grid.out)
    }
    jk.list <- seq_along(out.samplelist$context)
    info_lapply <- foreach(
        jk = 1:max(jk.list), .combine = "c", .packages = c(
            'methimpute'), .export = ".env") %dopar% 
        {
            runMethimputeJ(jk)
        }
    stopCluster(cl)
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using single-core method
#' @inheritParams makeMethimpute
#' @param out.samplelist DataFrame consisting of: contexts, methImpute files,
#'                     sample names (from methImpute files) and IDs 
#'                     corresponding to the binned genome. (DataFrame object)
#' @param merge_list List of DataFrames with binned regions; for each context
#'                   separately (list of DataFrame objects)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @param out.dir Path to the output directory. (char)
#' @param mincov Minimum read coverage; default as 0. (num; between 0 and 1)
#' @param if.Bismark Logical if Bismark inputs (CX reports in txt format)
#'                   are used. Default as FALSE. (logical)
#' @param FASTA.file Path to the FASTA file; required if Bismark outputs are
#'                   used. Default as NULL. (char)
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#' 
makeMethImpute_normal <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov,
        if.Bismark, FASTA.file)
{
    info_lapply <- lapply(seq_along(out.samplelist$context), function(jn) {
        refRegion <- list(reg.obs = merge_list[[out.samplelist$id[jn]]])
        message(
            "Running file: ", out.samplelist$methfn[jn], " for context: ",
            out.samplelist$context[jn], "\n")
        grid.out <- makeMethimpute(
            df = as.character(out.samplelist$file[jn]),
            context = out.samplelist$context[jn],
            refRegion = refRegion, fit.plot = FALSE,
            include.intermediate = include.intermediate,
            probability = "constrained", out.dir = out.dir,
            fit.name = paste0(
                basename(out.samplelist$methfn[jn]), "_",
                out.samplelist$context[jn]), name = basename(
                    out.samplelist$methfn[jn]), mincov = mincov, 
            if.Bismark = if.Bismark, FASTA.file = FASTA.file)
    })
}

#-----------------------------------------------------------------------------
#' Run jDMR on binned genome
#' this function runs a HMM model on a genome binned using 
#' a sliding/non-sliding window approach
#' @param out.dir Path to the output directory. (char)
#' @param window Bin size in non/sliding window approach. (num)
#' @param step Step size in non/sliding window approach. (num)
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @param contexts Vector of cytosine contexts; 
#'                 default as c('CG','CHG','CHH'). (char)
#' @param min.C Percentile threshold based on EDF. (num; between 0 and 100)
#' @param mincov Minimum read coverage; default as 0. (num; between 0 and 1)
#' @param include.intermediate Logical if intermediate calls should be used;
#'                             default as FALSE. (logical)
#' @param runName Character defining the name of the run; 
#'                default as GridGenome. (char)
#' @param numCores Number of cores to be used if for_each method should be 
#'                 performed; default as NULL. (num)
#' @param parallelApply Logical if futureapply method should be performed;
#'                      default as FALSE. (logical)
#' @param if.Bismark Logical if Bismark inputs (CX reports in txt format)
#'                   are used; default as FALSE. (logical)
#' @param FASTA.file Path to the FASTA file; required if Bismark outputs are
#'                   used; default as NULL. (char)
#' @import magrittr
#' @import future.apply
#' @import doParallel
#' @import future
#' @import foreach
#' @importFrom foreach %dopar%
#' @importFrom data.table fread fwrite
#' @importFrom stringr str_remove_all
#' @importFrom IRanges slidingWindows
#' @importFrom parallel makeCluster stopCluster
#' @importFrom rlang .data
#' @return Output methylome files for the regions using grid genome from 
#'         non/sliding window approach.
#' @export
#'
runjDMRgrid <- function(
        out.dir, window, step, samplelist, contexts=c('CG','CHG','CHH'), 
        min.C, mincov=0, include.intermediate=FALSE, runName='GridGenome',
        numCores = NULL, parallelApply = FALSE, if.Bismark = FALSE,
        FASTA.file = NULL)
{
    methimputefiles <- samplelist$file
    # Step 1: Bin genome according to the window and step size
    bin.genome.files <- binGenome(
        methimputefiles = methimputefiles, contexts = contexts, 
        window = window, step = step, min.C = min.C, out.dir = out.dir,
        runName = runName, if.Bismark = if.Bismark, FASTA.file = FASTA.file)
    # Step 2: Extract RData for binned regions and name them as contexts
    bin.select <- lapply(seq_along(contexts), function(x) {
        b <- bin.genome.files[grep(contexts[x], bin.genome.files)]
        if (length(b) != 0) {
            names(b) <- contexts[x]
            return(b) }})
    bin.select <- Filter(Negate(is.null), bin.select)
    names(bin.select) <- unlist(lapply(bin.select, function(x) names(x)[1]))
    # Step 3: Read all RData binned regions
    merge_list <- lapply(bin.select, function(x) dget(x))
    # Step 4: Create data frame consisting of: contexts, methImpute files,
    # sample names (from methImpute files)
    # and IDs corresponding to the binned genome
    out.samplelist <- expand.grid(file = samplelist$file, context = contexts)
    out.samplelist <- merge(out.samplelist, data.frame(
        context = names(bin.select), id = seq(1,length(names(bin.select)))))
    if (if.Bismark == FALSE) {
        out.samplelist$methfn<-unlist(lapply(out.samplelist$file,function(xi){
                gsub(".*methylome_|\\.txt|_All.txt$","",xi)}))} else {
        out.samplelist$methfn<-unlist(lapply(out.samplelist$file,function(xi){
                gsub("|\\.txt|.CX_report.txt$","",xi)}))}
    # Step 5: Run makeMethimpute for files in out.samplelist
    if (parallelApply == TRUE) {
        makeMethimpute_future(
            out.samplelist, merge_list, include.intermediate, out.dir, mincov,
            if.Bismark, FASTA.file)
    }
    if (is.numeric(numCores) == TRUE) {
        makeMethimpute_foreach(
            out.samplelist, merge_list, include.intermediate, out.dir, mincov,
            numCores, if.Bismark, FASTA.file)
    }
    if (is.null(numCores) & parallelApply == FALSE) {
        makeMethImpute_normal(
            out.samplelist, merge_list, include.intermediate, out.dir, mincov,
            if.Bismark, FASTA.file)
    }
}
