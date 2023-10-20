#-----------------------------------------------------------------------------
#' 
#' @param mylist
#' @param out.dir
#' @param runName
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
#' f
#' @param x1
#' @param window
#' @param step
#' @param gr
#' @param cyt_gr
#' @param contexts
#' @param min.C
#' @importFrom data.table fread
#' @importFrom stats quantile ecdf
#' @importFrom IRanges countOverlaps
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
#' Bin genome 
#' @param methimputefiles
#' @param contexts Vector of cytosine contexts. (character vector)
#' @param window Bin size. (numeric)
#' @param step Step size. (numeric)
#' @param min.C Percentile threshold based on EDF. (numeric between 0 and 100)
#' @param out.dir Output directory. (character)
#' @param runName
#' @importFrom data.table fread rbindlist fwrite
#' @importFrom stats quantile ecdf
#' @importFrom IRanges countOverlaps
#' @importFrom GenomicRanges GRanges
#'
binGenome <- function(
        methimputefiles, contexts, window, step, min.C, out.dir, runName)
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
    chr.names <- unique(meth.out$seqnames)
    chr.lengths <- unlist(
        lapply(chr.names, function(x) max(
            meth.out$start[meth.out$seqnames == x])))
    names(chr.lengths) <- chr.names
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
        lapply(
            collect.bins, function(x) export.bins(
                mylist=x, out.dir=out.dir,runName=runName))
        return(list.files(out.dir, pattern=paste0(
                ".*", runName, ".*\\.Rdata$"),full.names=TRUE))
    }
    message("Done!")
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using future_apply parallel method
#' @import magrittr
#' @import future
#' @import future.apply
#' @import doParallel
#' 
makeMethimpute_future <- function(
        out.filelist, merge_list, include.intermediate, out.dir, mincov)
{
    plan(multisession)
    info_lapply <- future_lapply(
        seq_along(out.filelist$context), function(j) {
            refRegion <- list(reg.obs = merge_list[[out.filelist$id[j]]])
            message("Running file: ", out.filelist$methfn[j], " for context: ",
                    out.filelist$context[j], "\n")
            grid.out <- makeMethimpute(
                df = as.character(out.filelist$file[j]),
                context = out.filelist$context[j],
                refRegion = refRegion, fit.plot = FALSE,
                include.intermediate = include.intermediate,
                probability = "constrained", out.dir = out.dir,
                fit.name = paste0(
                    basename(out.filelist$methfn[j]), "_",
                    out.filelist$context[j]),
                name = basename(out.filelist$methfn[j]), mincov = mincov)}, 
        future.seed = NULL)
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using future_apply parallel method
#' @import magrittr
#' @import foreach
#' @import doParallel
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' 
makeMethimpute_foreach <- function(
        out.filelist, merge_list, include.intermediate, out.dir, mincov, 
        numCores)
{
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    runMethimputeJ <- function(jj) {
        refRegion <- list(reg.obs = merge_list[[out.filelist$id[jj]]])
        message("Running file: ", out.filelist$methfn[jj],
                " for context: ", out.filelist$context[jj], "\n")
        grid.out <- makeMethimpute(
            df = as.character(out.filelist$file[jj]),
            context = out.filelist$context[jj],
            refRegion = refRegion, fit.plot = FALSE,
            include.intermediate = include.intermediate,
            probability = "constrained",out.dir = out.dir,
            fit.name = paste0(
                basename(out.filelist$methfn[jj]), "_",
                out.filelist$context[jj]),
            name = basename(out.filelist$methfn[jj]), mincov = mincov)
        return(grid.out)
    }
    info_lapply <- foreach(
        jk = seq_along(out.filelist$context), .combine = "c", .packages = c(
            'methimpute')) %dopar% 
        {
            runMethimputeJ(jk)
        }
    stopCluster(cl)
}

#-----------------------------------------------------------------------------
#' Run makeMethimpute using single-core method
#' 
makeMethImpute_normal <- function(
        out.filelist, merge_list, include.intermediate, out.dir, mincov)
{
    info_lapply <- lapply(seq_along(out.filelist$context), function(jn) {
        refRegion <- list(reg.obs = merge_list[[out.filelist$id[jn]]])
        message(
            "Running file: ", out.filelist$methfn[jn], " for context: ",
            out.filelist$context[jn], "\n")
        grid.out <- makeMethimpute(
            df = as.character(out.filelist$file[jn]),
            context = out.filelist$context[jn],
            refRegion = refRegion, fit.plot = FALSE,
            include.intermediate = include.intermediate,
            probability = "constrained", out.dir = out.dir,
            fit.name = paste0(
                basename(out.filelist$methfn[jn]), "_",
                out.filelist$context[jn]), name = basename(
                    out.filelist$methfn[jn]), mincov = mincov)
    })
}

#-----------------------------------------------------------------------------
#' Run jDMR on binned genome
#' this function runs a HMM model on a genome binned using 
#' a sliding/non-sliding window approach
#' @param out.dir Output directory. (character)
#' @param window Bin size. (numeric)
#' @param step Step size. (numeric)
#' @param samplefiles Path to the text file. (character)
#' @param contexts Vector of cytosine contexts. (character vector)
#' @param min.C Percentile threshold based on EDF. (numeric between 0 and 100)
#' @param mincov Minimum read coverage. (numeric value between 0 and 1)
#' @param include.intermediate A logical if intermediate calls should be used.
#' @param runName 
#' @param numCores Number of cores. (numeric)
#' @param parallelApply A logical if futureapply method should be performed.
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
#' @export
#'
runjDMRgrid <- function(
        out.dir, window, step, samplefiles, contexts=c('CG','CHG','CHH'), 
        min.C, mincov=0, include.intermediate=FALSE, runName='GridGenome',
        numCores = NULL, parallelApply = FALSE)
{
    methimputefiles <- fread(samplefiles)$file
    # Step 1: Bin genome according to the window and step size
    bin.genome.files <- binGenome(
        methimputefiles = methimputefiles, contexts = contexts, 
        window = window, step = step, min.C = min.C, out.dir = out.dir,
        runName = runName)
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
    # Step 4: Read names of methImpute files which needs to be binned
    filelist <- fread(samplefiles, header=TRUE)
    # Step 5: Create data frame consisting of: contexts, methImpute files,
    # sample names (from methImpute files)
    # and IDs corresponding to the binned genome
    out.filelist <- expand.grid(file = filelist$file, context = contexts)
    out.filelist <- merge(out.filelist, data.frame(
        context = names(bin.select), id = seq(1,length(names(bin.select)))))
    out.filelist$methfn <- unlist(lapply(
        seq_along(out.filelist$file),function(xi) {
            gsub(".*methylome_|\\.txt|_All.txt$", "", out.filelist$file[xi])}))
    # Step 6: Run makeMethimpute for files in out.filelist
    if (parallelApply == TRUE) {
        makeMethimpute_future(
            out.filelist, merge_list, include.intermediate, out.dir, mincov)
    }
    if (is.numeric(numCores) == TRUE)
    {
        makeMethimpute_foreach(
            out.filelist, merge_list, include.intermediate, out.dir, mincov,
            numCores)
    }
    if (is.null(numCores) & parallelApply == FALSE)
    {
        makeMethImpute_normal(
            out.filelist, merge_list, include.intermediate, out.dir, mincov)
    }
}