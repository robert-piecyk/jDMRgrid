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
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

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
binGenomeLoop <- function(x1, window, step, gr, cyt_gr, contexts, min.C,
                          min.C.type = c("percentile", "count"),
                          region.mode = c("grid", "cluster-grid", "cluster"),
                          max.gap = 150L, min.sites = 5L, gap.quantile = 0.9) {
    min.C.type  <- match.arg(min.C.type)
    region.mode <- match.arg(region.mode)
    window.size <- window[x1]
    step.size   <- step[x1]
    mydf <- data.frame(bin.size = window.size, step.size = step.size)

    results <- lapply(contexts, function(cx) {
        message("Extracting cytosines for ", cx, ".")
        ref_gr <- cyt_gr[which(cyt_gr$context == cx), ]

        ## ---- candidate region space -------------------------------------------------
        if (region.mode == "grid") {
            ## uniform tiling of the whole genome (the historical behaviour)
            binned.g <- slidingWindows(gr, window.size, step.size)
        } else {
            ## data-defined candidate space: cluster the cytosines of this context first.
            gapw <- max.gap
            if (identical(max.gap, "adaptive")) {
                d <- diff(sort(GenomicRanges::start(ref_gr)))
                d <- d[d > 0]
                gapw <- max(1L, as.integer(stats::quantile(d, gap.quantile,
                                                           na.rm = TRUE)))
                message("  adaptive gap for ", cx, ": ", gapw, " bp (q",
                        gap.quantile, " of inter-cytosine distances)")
            }
            ## Cluster on an UNSTRANDED copy. cyt_gr carries the strand of each cytosine,
            ## and GenomicRanges::reduce() groups within strand, so clustering ref_gr
            ## directly yields a "+" and a "-" copy of every interval. The duplicate
            ## coordinates survive into the region table and make methimpute fail while
            ## compiling results ("each range must have an end that is greater or equal to
            ## its start minus one"). A cytosine cluster is a genomic interval, not a
            ## strand-specific one, and "grid" mode is likewise unstranded.
            ref_u <- ref_gr
            GenomicRanges::strand(ref_u) <- "*"
            clus <- GenomicRanges::reduce(ref_u, min.gapwidth = gapw)
            clus <- clus[countOverlaps(clus, ref_u) >= min.sites]
            if (length(clus) < 100L)
                warning("context ", cx, " collapsed to ", length(clus),
                        " cluster(s) at a gap of ", gapw, " bp. Dense contexts (CHH in ",
                        "particular) merge into chromosome-scale blocks at a fixed gap; ",
                        "consider max.gap = \"adaptive\".", call. = FALSE)
            binned.g <- if (region.mode == "cluster") {
                ## the clusters themselves are the regions
                GenomicRanges::GRangesList(clus)
            } else {
                ## tile fixed windows, but only inside clusters
                slidingWindows(clus, window.size, step.size)
            }
        }
        ## grToDF(), not as.data.frame(). On Bioc devel GRanges extends List, so
        ## as.data.frame() lands in as.data.frame.List, which unlist()s the GRanges and
        ## errors. grToDF() builds the frame from accessors and works everywhere.
        dd <- grToDF(unlist(binned.g))
        names(dd) <- c("chr", "start", "end", "cluster.length", "strand")
        data_gr <- GRanges(seqnames = dd$chr,
                           ranges = IRanges(start = dd$start, end = dd$end),
                           clusterlen = dd$cluster.length)

        ## ---- cytosine-count filter ---------------------------------------------------
        dat.collect <- countOverlaps(data_gr, ref_gr)
        thr <- if (min.C.type == "count") {
            min.C
        } else {
            as.numeric(quantile(ecdf(dat.collect), min.C / 100))
        }
        keep <- which(dat.collect >= thr)
        if (length(keep) == 0L)
            stop("no regions survived for context ", cx, " (region.mode = \"",
                 region.mode, "\", min.C = ", min.C, " as ", min.C.type,
                 "). Downstream HMM fitting cannot proceed on an empty region set; ",
                 "relax min.C, min.sites, or max.gap.", call. = FALSE)
        non.empty.bins <- length(keep) / length(dat.collect)
        list(non.empty.bins, dd[keep, ])
    })
    new.one <- lapply(results, function(x) x[[2]])
    names(new.one) <- paste0(contexts, '_', window.size, '_', step.size)
    mydf <- cbind(mydf, data.frame(context = contexts, ratio = unlist(
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
        if.Bismark, FASTA.file,
        min.C.type = c("percentile", "count"),
        region.mode = c("grid", "cluster-grid", "cluster"),
        max.gap = 150L, min.sites = 5L, gap.quantile = 0.9)
{
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    # message about creating grid
    message('Creating grid...'); cyt.collect <- list()
    # from one of the methIMPUTE file extract all cytosines positions
    meth.out <- fread(methimputefiles[1],showProgress=FALSE)
    if (if.Bismark == FALSE) {
        all.cyt.pos <- meth.out[,c('seqnames','start','strand','context')]
    } else {
        all.cyt.pos <- meth.out[,c(1,2,3,6)]
    }
    colnames(all.cyt.pos) <- c('chr','pos','strand','context')
    # create a GRanges object from the cytosines positions
    cyt_gr <- GRanges(seqnames=all.cyt.pos$chr,ranges=IRanges(
        start=all.cyt.pos$pos, width=1),context=all.cyt.pos$context,
        strand=all.cyt.pos$strand)
    # get length of chromosomes
    if (is.null(FASTA.file)) {
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
        min.C.type  <- match.arg(min.C.type)
        region.mode <- match.arg(region.mode)
        message("Region mode: ", region.mode,
                " | min.C interpreted as: ", min.C.type)
        results <- lapply(
            seq_along(window), binGenomeLoop, window = window, step = step, 
            gr = gr, cyt_gr = cyt_gr, contexts = contexts, min.C = min.C,
            min.C.type = min.C.type, region.mode = region.mode,
            max.gap = max.gap, min.sites = min.sites,
            gap.quantile = gap.quantile)
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
#' @param cyt.pos.all GRanges object from extractCytosinesFromFASTA. (GRanges)
#' @import magrittr
#' @import future
#' @import future.apply
#' @import doParallel
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#' 
makeMethimpute_future <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov,
        if.Bismark, cyt.pos.all)
{
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

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
                if.Bismark = if.Bismark, cyt.pos.all = cyt.pos.all)}, 
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
#' @param cyt.pos.all GRanges object from extractCytosinesFromFASTA. (GRanges)
#' @import magrittr
#' @import foreach
#' @import doParallel
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster stopCluster                
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#' @export
#'
makeMethimpute_foreach <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov, 
        numCores, if.Bismark, cyt.pos.all)
{
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    cl <- makeCluster(numCores)
    ## FIX: cluster workers start fresh R sessions and do NOT inherit the caller's
    ## .libPaths(). If jDMRgrid is installed outside the default library -- or if an older
    ## copy sits in it -- the workers silently load a different version, which surfaces as
    ## "unused arguments" from inside the foreach loop. Propagate the caller's paths.
    .libs <- .libPaths()
    parallel::clusterExport(cl, ".libs", envir = environment())
    parallel::clusterEvalQ(cl, .libPaths(.libs))
    registerDoParallel(cl)
    runMethimputeJ <- function(jj) {
        refRegion <- list(reg.obs = merge_list[[out.samplelist$id[jj]]])
        message("Running file: ", out.samplelist$methfn[jj],
                " for context: ", out.samplelist$context[jj], "\n")
        grid.out <- makeMethimpute(
            df = as.character(out.samplelist$file[jj]),
            context = out.samplelist$context[jj],
            refRegion = refRegion, fit.plot = FALSE,
            include.intermediate = include.intermediate,
            probability = "constrained",out.dir = out.dir,
            fit.name = paste0(
                basename(out.samplelist$methfn[jj]), "_",
                out.samplelist$context[jj]),
            name = basename(out.samplelist$methfn[jj]), mincov = mincov,
            if.Bismark = if.Bismark, cyt.pos.all = cyt.pos.all)
        return(grid.out)
    }
    jk <- NULL
    info_lapply <- foreach(
        jk = seq_along(out.samplelist$context), .combine = "c", .packages = c(
            'methimpute'), .export = "jk") %dopar% 
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
#' @param cyt.pos.all GRanges object from extractCytosinesFromFASTA. (GRanges)
#' @return Methylome for regions taken out grid genome from non/sliding window
#'         approach
#' 
makeMethImpute_normal <- function(
        out.samplelist, merge_list, include.intermediate, out.dir, mincov,
        if.Bismark, cyt.pos.all)
{
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

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
            if.Bismark = if.Bismark, cyt.pos.all = cyt.pos.all)
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
#' @importFrom methimpute extractCytosinesFromFASTA
#' @importFrom rlang .data
#' @return Output methylome files for the regions using grid genome from 
#'         non/sliding window approach.
#' @param min.C.type How \code{min.C} is interpreted. Has no default and must be given:
#'   the meaning of \code{min.C} changed between versions and the two readings differ by
#'   orders of magnitude. \code{"percentile"} (an ECDF
#'   percentile of the per-window cytosine counts, the behaviour since v0.2.x) or
#'   \code{"count"} (an absolute number of cytosines, the behaviour of earlier versions).
#'   The same value means very different things under the two rules -- on A. thaliana at
#'   \code{min.C = 10} the percentile rule keeps 2 278 082 CG bins and the absolute rule
#'   1 012 452 -- so state it explicitly when reproducing older analyses. (character)
#' @param region.mode How candidate regions are defined. \code{"grid"} tiles the whole
#'   genome uniformly (default, historical behaviour). \code{"cluster-grid"} clusters the
#'   cytosines of each context first and tiles fixed windows only inside those clusters.
#'   \code{"cluster"} uses the clusters themselves as regions. (character)
#' @param max.gap Maximum distance between cytosines within one cluster, in bp, or
#'   \code{"adaptive"} to derive it per context from \code{gap.quantile} of the observed
#'   inter-cytosine distances. Cluster modes only. (numeric or character)
#' @param min.sites Minimum cytosines for a cluster to be kept. Cluster modes only. (numeric)
#' @param gap.quantile Quantile of inter-cytosine distances used when
#'   \code{max.gap = "adaptive"}. Lower values give tighter, more numerous clusters. (numeric)
#' @examples
#' ## region calls for one context, using the packaged toy data
#' sheet <- get(load(system.file("data", "listFiles2.RData", package = "jDMRgrid")))
#' sheet$file <- system.file("extdata", sheet$file, package = "jDMRgrid")
#'
#' out <- file.path(tempdir(), "jDMRgrid_grid")
#' runjDMRgrid(out.dir = out, window = 200, step = 50, samplelist = sheet,
#'             contexts = "CG", min.C = 10, min.C.type = "percentile",
#'             mincov = 0, include.intermediate = FALSE, runName = "example")
#' list.files(out)
#' @export
#'
runjDMRgrid <- function(
        out.dir, window, step, samplelist, contexts=c('CG','CHG','CHH'), 
        min.C, mincov=0, include.intermediate=FALSE, runName='GridGenome',
        numCores = NULL, parallelApply = FALSE, if.Bismark = FALSE,
        FASTA.file = NULL,
        min.C.type = c("percentile", "count"),
        region.mode = c("grid", "cluster-grid", "cluster"),
        max.gap = 150L, min.sites = 5L, gap.quantile = 0.9)
{
    ## min.C changed meaning between versions of this package -- it was an absolute count of
    ## cytosines and became a percentile of the observed distribution. The two differ by
    ## orders of magnitude, and neither reading announces itself in the output, so a script
    ## carried over from an earlier version silently filters something quite different from
    ## what its author intended. There is no safe default, so the interpretation must be
    ## stated. A percentile is also actively dangerous where the grid extends well beyond the
    ## data: if most windows hold no cytosines, a low percentile evaluates to zero and every
    ## empty window is retained, after which the HMM has nothing to fit.
    if (missing(min.C.type))
        stop("min.C.type must be given explicitly.\n",
             "  \"count\"      - min.C is an absolute number of cytosines per region\n",
             "  \"percentile\" - min.C is a percentile of the observed distribution\n",
             "The two differ by orders of magnitude and the meaning of min.C changed between ",
             "versions, so no default is applied.", call. = FALSE)

    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    methimputefiles <- samplelist$file
    bin.genome.files <- binGenome(
        methimputefiles = methimputefiles, contexts = contexts, 
        window = window, step = step, min.C = min.C, out.dir = out.dir,
        runName = runName, if.Bismark = if.Bismark, FASTA.file = FASTA.file,
        min.C.type = match.arg(min.C.type), region.mode = match.arg(region.mode),
        max.gap = max.gap, min.sites = min.sites, gap.quantile = gap.quantile)
    bin.select <- lapply(seq_along(contexts), function(x) {
        b <- bin.genome.files[grep(contexts[x], bin.genome.files)]
        if (length(b) != 0) {
            names(b) <- contexts[x]
            return(b) }})
    bin.select <- Filter(Negate(is.null), bin.select)
    names(bin.select) <- unlist(lapply(bin.select, function(x) names(x)[1]))
    merge_list <- lapply(bin.select, function(x) dget(x))
    ## FIX: expand.grid() still defaults to stringsAsFactors = TRUE in R 4.x (unlike
    ## data.frame(), which changed in 4.0), so `file` came back as a factor. gsub() coerces
    ## factors silently, which hid this; basename() does not.
    out.samplelist <- expand.grid(file = samplelist$file, context = contexts,
                                  stringsAsFactors = FALSE)
    out.samplelist <- merge(out.samplelist, data.frame(
        context = names(bin.select), id = seq(1,length(names(bin.select)))))
    if (if.Bismark == FALSE) {
        out.samplelist$methfn <- deriveSampleName(out.samplelist$file)
        cyt.pos.all<-NULL} else {
        out.samplelist$methfn <- deriveSampleName(out.samplelist$file,
                                                  if.Bismark = TRUE)
        cyt.pos.all<-extractCytosinesFromFASTA(
            file = FASTA.file, contexts = contexts)}
    if (parallelApply == TRUE) {
        makeMethimpute_future(
            out.samplelist=out.samplelist,merge_list=merge_list,mincov=mincov,
            include.intermediate=include.intermediate,out.dir=out.dir,
            if.Bismark=if.Bismark,cyt.pos.all=cyt.pos.all) }
    if (is.numeric(numCores) == TRUE) {
        makeMethimpute_foreach(
            out.samplelist=out.samplelist,merge_list=merge_list,mincov=mincov,
            include.intermediate=include.intermediate,out.dir=out.dir,
            numCores=numCores,if.Bismark=if.Bismark,cyt.pos.all=cyt.pos.all) }
    if (is.null(numCores) & parallelApply == FALSE) {
        makeMethImpute_normal(
            out.samplelist=out.samplelist,merge_list=merge_list,mincov=mincov,
            include.intermediate=include.intermediate,out.dir=out.dir,
            if.Bismark=if.Bismark,cyt.pos.all=cyt.pos.all) }
}
