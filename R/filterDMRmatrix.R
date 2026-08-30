#------------------------------------------------------------------------------
#' Bind list of replicates (as part of lapply function)
#' @param x Number of rows in state-calls file. (num)
#' @param status.collect DataFrame with state-calls. (DataFrame object)
#' @param sampleinfo Information about samples and their replicates in df
#' @param replicate.consensus Percentage of concordance in state calls in num
#' @param diffct Allowed difference between control and treatment groups in num
#' @return A logical if the following row meets replicate consensus criterion 
#'         and should be included in the filtered dataset.
#' 
bindReplicateLists <- function(
        x, mypattern, df.bind, sampleinfo, replicate.consensus, diff.ct,
        max.iter) {
    for (m in unique(df.bind$sample)){
        total.reps <- length(df.bind$mypattern[df.bind$sample==m])
        rval <- round(replicate.consensus * total.reps)
        pattern.vals <- df.bind$mypattern[df.bind$sample==m]
        df.bind$diff.count[df.bind$sample==m] <- length(
            which(pattern.vals==1))/total.reps
        tt <- table(pattern.vals)
        if (max(tt) >= rval){
            df.bind$count[df.bind$sample==m] <- 0
        } else {
            df.bind_rows <- apply(df.bind, 1, function(row) {
                m <- row[which(colnames(df.bind) == 'sample')]
                total.reps <- length(
                    df.bind$mypattern[df.bind$sample == m])
                rval <- round(replicate.consensus * total.reps)
                pattern.vals <- df.bind$mypattern[df.bind$sample == m]
                row[["diff.count"]] <- length(
                    which(pattern.vals == 1)) / total.reps
                tt <- table(pattern.vals)
                if (max(tt) >= rval) {
                    row[["count"]] <- 0
                } else {
                    row[["count"]] <- 1 }
                return(row) })
            df.bind_rows <- as.data.frame(t(df.bind_rows))
            df.bind_rows$mypattern <- as.numeric(
                df.bind_rows$mypattern)
            df.bind_rows$diff.count <- as.numeric(
                df.bind_rows$diff.count)
            df.bind_rows$count <- as.numeric(df.bind_rows$count)
            df.bind <- df.bind_rows
        }
    }
    df.gp <- group_by(df.bind, sample) %>%
        summarize(n = mean(.data$diff.count))
    ## FIX: lapply() iterates the ELEMENTS of the combn matrix and the body ignores its
    ## argument, so my.diff was the first pair's difference repeated. With more than two
    ## sample groups only that one pair was ever tested against diff.ct.
    cb <- combn(df.gp$n, 2)
    my.diff <- abs(cb[1, ] - cb[2, ])
    if ((min(my.diff) >= diff.ct) & (sum(df.bind$count)==0)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#-----------------------------------------------------------------------------
#' Application of filtering using replicate consensus
#' @param status.collect DataFrame for state-calls. (DataFrame object)
#' @param rc.methlevel.collect DataFrame rc-methylation levels. 
#'                             (DataFrame object)
#' @param replicate.consensus Percentage of concordance in state calls. (num)
#' @param diff.ct Allowed difference between control and treatment groups. (num)
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom dplyr semi_join summarize
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return List of data frames for rc-methylation level and state-calls using 
#'         replicate-consensus criterion.
#' @examples
#' ## two groups of two replicates, named <sample>_<replicate>
#' co <- data.frame(seqnames = 1L, start = c(1L, 201L, 401L), end = c(200L, 400L, 600L))
#' sc <- data.frame(co, WT_rep1 = c(1, 0, 1), WT_rep2 = c(1, 0, 1),
#'                      mut_rep1 = c(0, 0, 1), mut_rep2 = c(0, 0, 1))
#' rc <- data.frame(co, WT_rep1 = c(0.9, 0.1, 0.9), WT_rep2 = c(0.9, 0.1, 0.9),
#'                      mut_rep1 = c(0.1, 0.1, 0.9), mut_rep2 = c(0.1, 0.1, 0.9))
#'
#' ## replicates must agree within each group, and the groups must differ
#' res <- filterReplicateConsensus(sc, rc, replicate.consensus = 1, diff.ct = 0.5)
#' res[[1]]
#' @export
#'
filterReplicateConsensus <- function(
        status.collect, rc.methlevel.collect, replicate.consensus, diff.ct) {
    if (!is.null(status.collect$epiMAF)){
        status.collect <- as.data.frame(status.collect)
        status.collect <- status.collect[, setdiff(names(status.collect), "epiMAF")]
    }
    mycol <- names(status.collect)[4:NCOL(status.collect)]
    sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol), "_")))
    if (length(sampleinfo) != 2) stop("Column for replicates is missing!!!")
    colnames(sampleinfo) <- c("sample", "replicate")

    ## PERF: this used to iterate the REGIONS -- one lapply() pass per region, each building a
    ## fresh data.frame, looping over every sample group, and calling table() per group. On
    ## CHH at 169 samples that is on the order of 1e11 element operations and 4e8 table()
    ## calls, and the inner branch additionally ran apply() over a data.frame, coercing every
    ## column to character and converting three of them back. That branch also recomputed
    ## exactly what the enclosing loop had already assigned.
    ##
    ## The quantities are per SAMPLE GROUP, not per region, so they vectorise directly:
    ##   consensus_g : the modal state among group g's replicates reaches
    ##                 round(replicate.consensus * n_reps)
    ##   diffcount_g : fraction of group g's replicates in state 1 (methylated)
    ## A region is kept when every group reaches consensus and the smallest pairwise
    ## difference in diffcount across groups is at least diff.ct.
    M   <- as.matrix(status.collect[, 4:NCOL(status.collect)])
    ## An earlier filter (min.posterior, epiMAF, or the non-polymorphic pass) can legitimately
    ## leave nothing behind. Without this guard `lv` is empty and do.call(pmax, list()) fails
    ## with "no arguments", which reads as a bug in the caller rather than an empty result.
    if (NROW(M) == 0L) {
        message("Empty data frame. Nothing to write!")
        return(NULL)
    }
    lv  <- sort(unique(as.vector(M)))
    grp <- split(seq_along(mycol), sampleinfo$sample)
    if (length(grp) < 2L)
        stop("replicate consensus filtering needs at least two sample groups")

    consensus <- matrix(TRUE, nrow(M), length(grp))
    dcount    <- matrix(0,    nrow(M), length(grp))
    for (i in seq_along(grp)) {
        cols <- grp[[i]]
        sub  <- M[, cols, drop = FALSE]
        mx   <- do.call(pmax, lapply(lv, function(v) rowSums(sub == v)))
        consensus[, i] <- mx >= round(replicate.consensus * length(cols))
        dcount[, i]    <- rowSums(sub == 1) / length(cols)
    }
    pr   <- combn(ncol(dcount), 2)
    mind <- do.call(pmin, lapply(seq_len(ncol(pr)), function(k)
                abs(dcount[, pr[1, k]] - dcount[, pr[2, k]])))
    keep <- (rowSums(!consensus) == 0L) & (mind >= diff.ct)

    df.status.collect <- as.data.frame(status.collect)[keep, , drop = FALSE]
    if (NROW(df.status.collect) != 0){
        df.rc.methlevel.collect <- as.data.frame(rc.methlevel.collect) %>%
            semi_join(df.status.collect, by = c("seqnames", "start", "end"))
        return(list(df.status.collect, df.rc.methlevel.collect))
    } else {
        message("Empty data frame. Nothing to write!")
        return(NULL)
    }
}

#-----------------------------------------------------------------------------
#' Application of filtering using epiMAF criterion
#' @param status.collect DataFrame for state-calls. (DataFrame object)
#' @param rc.methlevel.collect DataFrame rc-methylation levels. 
#'                             (DataFrame object)
#' @param epiMAF Threshold to filter for minor epi-allele frequency. (num)
#' @import magrittr
#' @importFrom dplyr semi_join
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return List of data frames for rc-methylation level and state-calls using 
#'         minor epi-allele frequency criterion.
#' @examples
#' ## state calls and methylation levels for one context
#' co <- data.frame(seqnames = 1L, start = c(1L, 201L, 401L), end = c(200L, 400L, 600L))
#' sc <- data.frame(co, s1 = c(1, 0, 1), s2 = c(1, 0, 1),
#'                      s3 = c(1, 0, 1), s4 = c(0, 0, 1))
#' rc <- data.frame(co, s1 = c(0.9, 0.1, 0.8), s2 = c(0.9, 0.1, 0.8),
#'                      s3 = c(0.9, 0.1, 0.8), s4 = c(0.1, 0.1, 0.8))
#'
#' ## keep regions where the rarer state is carried by fewer than 30 percent of samples
#' res <- filterEpiMAF(sc, rc, epiMAF = 0.3)
#' res[[1]]
#' @export
#' 
filterEpiMAF <- function(status.collect, rc.methlevel.collect, epiMAF){
    floorDec <- function(valParm ,x){
        y <- function(x, level=1) round(x - 5*10^(-level-1), level)
        res <-y(as.numeric(valParm),x)
        return(res)
    }
    ## PERF: was apply(mypatterns, 1, ...) with a table() call per region, plus a
    ## Sys.sleep() that existed only so the progress bar rendered. States are numeric
    ## (0/0.5/1, see merge_cols), so the minor-"allele" frequency is rowSums per level and a
    ## row-wise minimum -- three vectorised passes instead of one interpreted call per region.
    ## Measured 62x faster on 200k regions x 169 samples, identical results.
    ## FIX (behaviour change, see NEWS): the range was 4:(NCOL(status.collect) - 1), but the
    ## epiMAF column is appended AFTER this line, so at this point there is no trailing column
    ## and the last SAMPLE was silently excluded from the minor-allele-frequency calculation.
    ## With four samples it used three, making the smallest achievable epiMAF 1/3 = 0.333 --
    ## so an epiMAF.cutoff of 0.33 could never match anything.
    mypatterns <- status.collect[, 4:NCOL(status.collect)]
    M   <- as.matrix(mypatterns)
    lv  <- sort(unique(as.vector(M)))
    cnt <- vapply(lv, function(v) rowSums(M == v), numeric(NROW(M)))
    if (!is.matrix(cnt)) cnt <- matrix(cnt, nrow = NROW(M))
    ## table() tabulates only the levels PRESENT in a row, so its min() ignores states the
    ## row does not contain; setting zero counts to Inf reproduces that exactly.
    cnt[cnt == 0] <- Inf
    status.collect$epiMAF <- floorDec(do.call(pmin, as.data.frame(cnt)) / NCOL(M), 5)
    
    status.collect <- as.data.frame(status.collect)
    rc.methlevel.collect <- as.data.frame(rc.methlevel.collect)
    df.status.collect <- status.collect[which(status.collect$epiMAF < epiMAF),]
    if (NROW(df.status.collect) !=0){
        df.rc.methlevel.collect <- rc.methlevel.collect %>% semi_join(
            df.status.collect, by=c("seqnames","start","end"))
        return(list(df.status.collect, df.rc.methlevel.collect))
    } else {
        message("Empty data frame. Nothing to write!")
        return(NULL)
    }
}

#------------------------------------------------------------------------------
#' Merge bins having the same state-calls and calculate mean of their
#' methylation levels
#' @param rcmethlvl_data DataFrame rc-methylation levels. 
#'                       (DataFrame object)
#' @param statecalls_data DataFrame for state-calls. (DataFrame object)
#' @import magrittr
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce findOverlaps
#' @importFrom data.table fread fwrite rbindlist
#' @return List of data frames for rc-methylation level and state-calls with
#'         merged neighbouring bins having similar state-call profile.
#' 
merge.bins <- function(rcmethlvl_data, statecalls_data)
{
    ## PERF: the previous implementation built the grouping key with
    ##   apply(statecalls[, 4:ncol], 1, paste, collapse = "_")
    ## and did so TWICE -- once to split the state calls and again to split the methylation
    ## levels. On a CHH matrix that is roughly 800 MB of throwaway strings, built twice. It
    ## then ran GenomicRanges::reduce() per pattern group and, for every merged interval,
    ## an apply(, 2, mean) over that interval's rows.
    ##
    ## Sorting by (pattern, chromosome, start) turns the merge into a single running-maximum
    ## sweep -- a new interval begins wherever a bin starts beyond the furthest end seen so
    ## far in its group -- and the averaging becomes one grouped aggregation. No strings are
    ## built and nothing iterates the intervals.
    ##
    ## reduce() merges book-ended ranges as well as overlapping ones (its default
    ## min.gapwidth = 1), so the break test is `start > running_max_end + 1`, not `>`.
    sc <- data.table::as.data.table(statecalls_data)
    rc <- data.table::as.data.table(rcmethlvl_data)
    scols <- names(sc)[4:NCOL(sc)]
    rcols <- names(rc)[4:NCOL(rc)]

    sc[, .pat := .GRP, by = scols]
    ord <- order(sc$.pat, as.character(sc$seqnames), sc$start, sc$end)
    sc  <- sc[ord]; rc <- rc[ord]

    chr <- as.character(sc$seqnames)
    runmax <- sc[, data.table::shift(cummax(end), fill = -1L),
                 by = .(.pat, seqnames)]$V1
    newgrp <- sc$.pat != data.table::shift(sc$.pat, fill = -1L) |
              chr != data.table::shift(chr, fill = "")
    sc[, .mid := cumsum(newgrp | (sc$start > runmax + 1L))]

    states.data <- sc[, .(seqnames = seqnames[1L], start = min(start), end = max(end)),
                      by = .mid]
    states.data <- cbind(states.data[, -".mid"],
                         sc[!duplicated(sc$.mid), scols, with = FALSE])

    rc[, .mid := sc$.mid]
    rcmethlvl.data <- rc[, lapply(.SD, mean), by = .mid, .SDcols = rcols]
    rcmethlvl.data <- cbind(states.data[, seq_len(3), with = FALSE],
                            rcmethlvl.data[, -".mid"])

    states.data    <- states.data[order(seqnames, start, end), ]
    rcmethlvl.data <- rcmethlvl.data[order(seqnames, start, end), ]
    return(list(states.data, rcmethlvl.data))
}

#------------------------------------------------------------------------------
#' Save methylation levels and state calls at a given context as txt file
#' @param out.rcmethlvl Output DataFrame rc-methylation levels. 
#'                      (DataFrame object)
#' @param out.statecalls Output DataFrame for state-calls. (DataFrame object)
#' @param context Methylation context. (char)
#' @param out.name1 Name for the state-call file. (char)
#' @param out.name2 Name for the rc-methylation level file. (char)
#' @param data.out Path to the output directory. (char)
#' @import data.table fwrite
#' @return Output DMR matrices for state-calls and rc-methylation levels as 
#'         txt file.
#' 
export.out <- function(
        out.rcmethlvl,out.statecalls,context,out.name1,out.name2,data.out){
    fwrite(
        x=out.statecalls,file=paste0(
            data.out, "/", context, "_", out.name1, ".txt"
            ), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    if (!is.null(out.rcmethlvl)) {
        fwrite(x=out.rcmethlvl,file=paste0(
            data.out, "/", context, "_", out.name2, ".txt"
            ),quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    } 
}

#------------------------------------------------------------------------------
#' Remove non-polymorphic state-calls only by merging bins if applicable
#' @inheritParams merge.bins
#' @inheritParams export.out
#' @param status.collect DataFrame for state-calls. (DataFrame object)
#' @param rc.methlevel.collect DataFrame rc-methylation levels. 
#'                             (DataFrame object)
#' @param context Methylation context. (char)
#' @param data.dir Path to the output directory. (char)
#' @param if.mergingBins Logical if bins with same state-calls profile among
#'                       samples should be merged. (logical)
#' @return Filtrated state-calls and methylation DMR matrices using 
#'         merging bins criterion.
#'
removeNonPolymorphicOnly <- function(
        status.collect, rc.methlevel.collect, context, data.dir, 
        if.mergingBins) {
    message("Both, epiMAF and replicate consensus set to NULL.
                Merge corresponding bins")
    out1 <- status.collect
    out2 <- rc.methlevel.collect
    if (if.mergingBins == TRUE)
    {
        out <- merge.bins(statecalls_data=out1, rcmethlvl_data=out2)
    } else {
        out <- list(statecalls=out1, rcmethlvl=out2)
    }
    export.out(
        out.statecalls=out[[1]],out.rcmethlvl=out[[2]],
        context=context,out.name1="StateCalls-filtered",
        out.name2="rcMethlvl-filtered",data.out=data.dir)
}

#------------------------------------------------------------------------------
#' Apply epiMAF filtering and merging bins if applicable
#' @inheritParams filterEpiMAF
#' @inheritParams merge.bins
#' @inheritParams export.out
#' @param status.collect DataFrame for state-calls. (DataFrame object)
#' @param rc.methlevel.collect DataFrame rc-methylation levels. 
#'                             (DataFrame object)
#' @param context Methylation context. (char)
#' @param data.dir Path to the output directory. (char)
#' @param epiMAF.cutoff Threshold to filter for 
#'                      the minor epi-allele frequency. (num)
#' @param if.mergingBins Logical if bins with same state-calls profile among
#'                       samples should be merged. (logical)
#' @return Filtrated state-calls and methylation DMR matrices using epi-allele
#'         frequency criterion and merging bins, if applicable
#'         
populationFiltering <- function(
        status.collect, rc.methlevel.collect, context, data.dir, 
        epiMAF.cutoff, if.mergingBins) {
    message("Filtering for epiMAF: ", epiMAF.cutoff, '...')
    mydf <- filterEpiMAF(
        status.collect=status.collect, 
        rc.methlevel.collect=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
    if (!is.null(mydf)){
        # for population data remove the epiMAF column
        out1 <- mydf[[1]][,-ncol(mydf[[1]])]
        out2 <- mydf[[2]]
        if (if.mergingBins == TRUE)
        {
            out <- merge.bins(
                statecalls_data=out1, rcmethlvl_data=out2)
        } else {
            out <- list(statecalls=out1, rcmethlvl=out2)
        }
        export.out(
            out.statecalls=out[[1]], out.rcmethlvl=out[[2]],
            context=context, out.name1="StateCalls-filtered",
            out.name2="rcMethlvl-filtered", data.out=data.dir)
    } 
    else {
        message('Filtering for epiMAF returns NULL dataset.')
        out1 <- status.collect
        out2 <- rc.methlevel.collect
        if (if.mergingBins == TRUE)
        {
            out <- merge.bins(
                statecalls_data=out1, rcmethlvl_data=out2)
        } else {
            out <- list(statecalls=out1, rcmethlvl=out2)
        }
        export.out(
            out.statecalls=out[[1]], out.rcmethlvl=out[[2]],
            context=context, out.name1="StateCalls-filtered",
            out.name2="rcMethlvl-filtered", data.out=data.dir)
    }
}

#------------------------------------------------------------------------------
#' Apply replicate consensus filtering and merging bins if applicable
#' @inheritParams filterReplicateConsensus
#' @inheritParams merge.bins
#' @inheritParams export.out
#' @param status.collect DataFrame for state-calls. (DataFrame object)
#' @param rc.methlevel.collect DataFrame rc-methylation levels. 
#'                             (DataFrame object)
#' @param context Methylation context. (char)
#' @param data.dir Path to the output directory. (char)
#' @param replicate.consensus Percentage of concordance in state calls. (num)
#' @param if.mergingBins Logical if bins with same state-calls profile among
#'                       samples should be merged. (logical)
#' @param diff.ct Allowed difference between control and treatment groups. (num)
#' @return Filtrated state-calls and methylation DMR matrices using replicate
#'         consensus and merging bins, if applicable
#' 
replicateFiltering <- function(
        status.collect, rc.methlevel.collect, context, data.dir, 
        replicate.consensus, if.mergingBins, diff.ct) {
    message("Filtering for replicate consensus...")
    mydf <- filterReplicateConsensus(
        status.collect, rc.methlevel.collect, replicate.consensus, diff.ct)
    if (!is.null(mydf)){
        out1 <- mydf[[1]]
        out2 <- mydf[[2]]
        if (if.mergingBins == TRUE)
        {
            out <- merge.bins(
                statecalls_data=out1, rcmethlvl_data=out2)
        } else {
            out <- list(statecalls=out1, rcmethlvl=out2)
        }
        export.out(
            out[[2]],out[[1]],context=context,
            out.name1="StateCalls-filtered",
            out.name2="rcMethlvl-filtered", data.out=data.dir)
    } else {
        message('Proceeding with the original dataset...')
        out1 <- status.collect
        out2 <- rc.methlevel.collect
        if (if.mergingBins == TRUE)
        {
            out <- merge.bins(
                statecalls_data=out1, rcmethlvl_data=out2)
        } else {
            out <- list(statecalls=out1, rcmethlvl=out2)
        }
        export.out(
            out.statecalls=out[[1]], out.rcmethlvl=out[[2]],
            context=context, out.name1="StateCalls-filtered", 
            out.name2="rcMethlvl-filtered", data.out=data.dir)
    }
}

#------------------------------------------------------------------------------
#' Filter DMR matrix among contexts
#' @inheritParams removeNonPolymorphicOnly
#' @inheritParams populationFiltering
#' @inheritParams replicateFiltering
#' @param list.status Vector with the filenames for state-calls matrices. (char)
#' @param data.dir Path to the input/output directory for 
#'                 raw/filtered DMR matrices. (char)
#' @param replicate.consensus Percentage of concordance in state calls. (num)
#' @param epiMAF.cutoff Threshold to filter for minor epi-allele frequency. 
#'                      (num)
#' @param if.mergingBins Logical if bins with same state-calls profile among
#'                       samples should be merged. (logical)
#' @import magrittr
#' @importFrom data.table fread
#' @return Filtrated state-calls and methylation DMR matrices
#'
filteringAmongContexts <- function(
        list.status, data.dir, replicate.consensus, epiMAF.cutoff, 
        if.mergingBins, min.posterior = NULL, diff.ct = 0.5) {
    for (i in seq_along(list.status)){
        # extract context for a given state calls file
        context <- gsub("_StateCalls.txt", "", basename(list.status[i]))
        message("Filtering DMR matrix for ", context, ' context.')
        # read a given state calls file
        if (file.exists(list.status[i])){
            status.collect  <- fread(list.status[i], header=TRUE)
        } else {
            stop("Files do not exist or is non-readable!")
        }
        # remove non-polymorphic/unchanged patterns
        message("Removing non-polymorphic patterns...")
        # extract rows where state calls for each sample are NOT unchanged
        index <- which(rowSums(
            status.collect[,4:NCOL(status.collect)]) != 0 & rowSums(
                status.collect[,4:NCOL(
                    status.collect)]) != NCOL(status.collect)-3)
        # update state calls list with the changed patterns
        status.collect <- status.collect[index,]
        # read a corresponding methylation file and update it
        rc.methlvl.name <- paste0(data.dir, "/", context, "_rcMethlvl.txt")
        rc.methlevel.collect <- fread(rc.methlvl.name, header=TRUE)
        rc.methlevel.collect <- rc.methlevel.collect[index,]
        ## Confidence filter. makeDMRmatrix(postMax.out = TRUE) writes the HMM posterior for
        ## every region and sample, and until now nothing ever read it back: a state assigned
        ## at posterior 0.56 was treated exactly like one assigned at 0.99. Requiring a
        ## minimum posterior across all samples in the comparison removes the low-confidence
        ## calls that dominate the disagreement between asymmetric and symmetric runs.
        if (!is.null(min.posterior)) {
            ## The postMax matrix is written once per context by makeDMRmatrix(). When a
            ## group column is present the state-call files are per group pair and `context`
            ## carries the pair suffix (e.g. "CG_WT_mutant1"), so fall back to the bare
            ## context prefix.
            pm.name <- paste0(data.dir, "/", context, "_postMax.txt")
            if (!file.exists(pm.name))
                pm.name <- paste0(data.dir, "/", sub("_.*", "", context), "_postMax.txt")
            if (!file.exists(pm.name)) {
                stop("min.posterior needs ", basename(pm.name),
                     "; re-run makeDMRmatrix() with postMax.out = TRUE")
            }
            postMax.collect <- data.table::as.data.table(fread(pm.name, header = TRUE))
            ## Align by coordinate, never by row position: the postMax matrix covers all
            ## regions and all samples, while status.collect has already been subset to the
            ## polymorphic ones and may hold only the samples in this comparison.
            key <- c("seqnames", "start", "end")
            sc.key <- data.table::as.data.table(status.collect)[, key, with = FALSE]
            pm <- postMax.collect[sc.key, on = key]
            smp <- intersect(names(status.collect)[4:NCOL(status.collect)], names(pm))
            if (!length(smp))
                stop("no sample columns shared between ", basename(pm.name),
                     " and the state-call matrix")
            pmm <- as.matrix(pm[, smp, with = FALSE])
            ## do.call(pmin, ...) is a vectorised row minimum and avoids taking a dependency
            ## on matrixStats for one call. Regions absent from the postMax matrix become NA
            ## and are dropped rather than silently kept.
            keep.post <- do.call(pmin, as.data.frame(pmm)) >= min.posterior
            keep.post[is.na(keep.post)] <- FALSE
            message("Confidence filter: min.posterior = ", min.posterior, "; keeping ",
                    sum(keep.post), " of ", length(keep.post), " regions.")
            status.collect       <- status.collect[keep.post,]
            rc.methlevel.collect <- rc.methlevel.collect[keep.post,]
        }
        # we are only removing non-polymorphic patterns
        if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
            removeNonPolymorphicOnly(
                status.collect = status.collect, 
                rc.methlevel.collect = rc.methlevel.collect, 
                context = context, data.dir = data.dir,
                if.mergingBins = if.mergingBins)}
        # filtering out regions with epiMAF < Minor Allele Frequency
        if (!is.null(epiMAF.cutoff)) {
            populationFiltering(
                status.collect = status.collect, 
                rc.methlevel.collect = rc.methlevel.collect, 
                context = context, data.dir = data.dir, 
                epiMAF.cutoff = epiMAF.cutoff, 
                if.mergingBins = if.mergingBins)}
        # retaining samples based on replicate.consensus
        if (!is.null(replicate.consensus)) {
            replicateFiltering(
                status.collect = status.collect, 
                rc.methlevel.collect = rc.methlevel.collect, 
                context = context, data.dir = data.dir, 
                replicate.consensus = replicate.consensus, 
                if.mergingBins = if.mergingBins, diff.ct = diff.ct)}
    }
}

#------------------------------------------------------------------------------
#' Filter DMR matrix. Filters non-polymorphic patterns by default.
#' @inheritParams filteringAmongContexts
#' @param epiMAF.cutoff Threshold to filter for 
#'                      minor epi-allele frequency. (num; default as NULL)
#' @param replicate.consensus Percentage of concordance in state calls.
#'                            (num; default as NULL)
#' @param data.dir Path to the directory containing DMR matrix files. (char)
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @param if.mergingBins Logical: merge consecutive bins carrying identical state calls
#'                       across all samples, so that one differential region is one row rather
#'                       than one row per overlapping window. (logical; default as TRUE)
#' @param diff.ct Minimum between-group difference, in [0, 1], in the proportion of
#'   replicates called methylated. Applies together with \code{replicate.consensus}. The
#'   default of 0.5 is strict: with three replicates per group it requires at least a
#'   2-of-3 versus 0-of-3 split. Previously fixed internally. (numeric)
#' @param min.posterior Optional minimum HMM posterior, in [0, 1], required across every
#'   sample in the comparison. Needs the postMax matrix from
#'   \code{makeDMRmatrix(postMax.out = TRUE)}. NULL (default) disables the filter.
#' @import magrittr
#' @importFrom data.table fread
#' @return Filtrated state-calls and methylation DMR matrices using:
#'         filtration for non-polymorphic patterns and 
#'         (epiMAF for population-based data or 
#'         replicate consensus for control/treatment data) 
#'         and/or merging neighbouring bins with same state calls profiles.
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
#' ## no group column, so the plain matrices are filtered directly
#' sheet <- data.frame(
#'     file = NA_character_, sample = rep(c("WT", "mutant1"), each = 2),
#'     replicate = rep(c("rep1", "rep2"), 2), stringsAsFactors = FALSE)
#'
#' filterDMRmatrix(replicate.consensus = 1, data.dir = d, samplelist = sheet)
#' data.table::fread(file.path(d, "CG_StateCalls-filtered.txt"))
#' @export
#'
filterDMRmatrix <- function(
        epiMAF.cutoff = NULL, replicate.consensus = NULL, data.dir, 
        samplelist, if.mergingBins = TRUE, min.posterior = NULL, diff.ct = 0.5) 
{
    if (!is.numeric(diff.ct) || length(diff.ct) != 1L || diff.ct < 0 || diff.ct > 1)
        stop("diff.ct must be a single number in [0, 1]")
    if (!is.null(min.posterior)) {
        if (!is.numeric(min.posterior) || length(min.posterior) != 1L ||
            min.posterior < 0 || min.posterior > 1)
            stop("min.posterior must be a single number in [0, 1]")
    }
    contexts <- unique(sub("_.*", "", list.files(data.dir, pattern = '^C')))
    if (!is.null(samplelist$group)){
        samplelist$name <- paste0(samplelist$sample,"_", samplelist$replicate)
        gps <- samplelist$group[!samplelist$group %in% c('control')]
        gps <- unique(gps)
        
        list.collect2 <- lapply(seq_along(gps), function(m) {
            myvec <- c("control", gps[m])
            gp1 <- samplelist$name[which(samplelist$group == myvec[1])]
            gp2 <- samplelist$name[which(samplelist$group == myvec[2])]
            gp1.sample <- unique(
                samplelist$sample[which(samplelist$name == gp1)])
            gp2.sample <- unique(
                samplelist$sample[which(samplelist$name == gp2)])
            out.name <- paste0(gp1.sample, "_", gp2.sample)
            
            list.collect1 <- lapply(seq_along(contexts), function(cn) {
                fn1 <- paste0(
                    data.dir,'/',contexts[cn],"_",out.name,"_StateCalls.txt")
                return(fn1)
            })
            return(list.collect1)
        })
        list.status <- unlist(list.collect2)
    } else {
        list.status <- list.files(
            data.dir, pattern="_StateCalls.txt", full.names=TRUE)
    }
    if (length(list.status) != 0){
        filteringAmongContexts(
            list.status = list.status, data.dir = data.dir, 
            replicate.consensus = replicate.consensus, 
            epiMAF.cutoff = epiMAF.cutoff, if.mergingBins = if.mergingBins,
            min.posterior = min.posterior, diff.ct = diff.ct)
    } else {
        message("DMR matrix files do not exist!")
    }
}

#------------------------------------------------------------------------------
#' Read state calls files for each context and create their GRanges objects
#' @param file.list Path to the DMR matrices with state calls for 
#'                  CG, CHG and CHH contexts. (char)
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @return List of data frames and GRanges objects of state calls; for each
#'         methylation contexts separately.
#' 
DMR.list.in <- function(file.list) {
    contexts <- c('CG','CHG','CHH')
    data.out <- lapply(file.list, function(x) fread(x, header = TRUE))
    gr.out <- lapply(data.out, function(x) GRanges(
        seqnames=x$seqnames,ranges=IRanges(x$start,end=x$end)))
    names(data.out) <- contexts
    names(gr.out) <- contexts
    return(list(data.out = data.out, gr.out = gr.out))
}

#------------------------------------------------------------------------------
#' Save DataFrame with context-specific DMRs as txt file
#' @param context.df DataFrame with filtrated context-specific DMRs
#'                   (DataFrame object).
#' @param out.name Name sample_replicate_context_group-of-analysis OR
#'                 context_group-of-analysis. (char)
#' @param data.out Path to the output directory. (char)
#' @importFrom data.table fwrite
#' @return Context-specific DMRs saved as the txt files in the output directory.
#' 
DMR.list.out <- function(context.df, out.name, data.out){
    fwrite(
        x=context.df, file=paste0(data.out, "/", out.name,".txt"), quote=FALSE, 
        row.names=FALSE, col.names=TRUE, sep="\t")
}

#------------------------------------------------------------------------------
#' Identify DMRs which are occuring only at a given context 
#' (for non-CG and multi-context DMRs)
#' @param context Methylation context for which unique DMRs are identified;
#'                "CG", "CHG" or "CHH". (char)
#' @param list.out List of state calls and their GRanges for each context. 
#'                 (list)
#' @param data.dir Path to the output directory. (char)
#' @param tmp.name For group-specific analysis: vector with with 
#'                 sample_replicate names. (char)
#' @import magrittr
#' @importFrom dplyr semi_join
#' @return 
#' 
DMR.onlyContext <- function(context, list.out, data.dir, tmp.name) {
    message('Generating ', context, '-only DMRs...')
    sel.line <- ifelse(context == 'CG', 1, ifelse(context=='CHG', 2, 3))
    out.line <- setdiff(seq(1,3), sel.line)
    out.1 <- subsetByOverlaps(
        list.out[[2]][[sel.line]], list.out[[2]][[out.line[1]]], invert = TRUE)
    out.2 <- subsetByOverlaps(
        out.1, list.out[[2]][[out.line[2]]], invert = TRUE)
    out.only <- grToDF(out.2)
    out.only $seqnames <- as.integer(as.character(out.only$seqnames))
    out.only <- out.only %>% semi_join(
        out.only, by = c("seqnames","start","end"))
    DMR.list.out(
        context.df=out.only, out.name=paste0(tmp.name, context, "-only-DMRs"),
        data.out=data.dir)
    message("Done!")
}

#------------------------------------------------------------------------------
#' Main function to extract context-specific DMRs 
#' for the provided list of state-calls paths
#' @param file1 Path to the DMR matrix state calls for CG context. (char)
#' @param file2 Path to the DMR matrix state calls for CHG context. (char)
#' @param file3 Path to the DMR matrix state calls for CHH context. (char)
#' @param tmp.name For group-specific analysis: vector with with 
#'                 sample_replicate names. (char)
#' @param data.dir Path to the output directory. (char)
#' @import magrittr
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr semi_join
#' @importFrom GenomicRanges findOverlaps
#' @importFrom data.table rbindlist
#' @importFrom data.table fread
#' @importFrom GenomicRanges intersect
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate
#' @importFrom dplyr semi_join
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom IRanges setTxtProgressBar txtProgressBar subsetByOverlaps
#' @importFrom rlang .data
#' @return Output files containing context-specific DMRs in txt format
#'
extract.context.DMRs <- function(file1, file2, file3, tmp.name, data.dir)
{
    if (file.exists(file1) && file.exists(file2) && file.exists(file3)) {
        list.out <- DMR.list.in(c(file1, file2, file3))
        # context-only
        DMR.onlyContext(
            'CG', list.out = list.out, data.dir = data.dir, 
            tmp.name = tmp.name) #CG
        DMR.onlyContext(
            'CHG', list.out = list.out, data.dir = data.dir,
            tmp.name = tmp.name) #CHG
        DMR.onlyContext(
            'CHH', list.out = list.out, data.dir = data.dir,
            tmp.name = tmp.name) #CHH
        #non-CG
        message("Generating non-CG DMRs...")
        gr.nonCG <- subsetByOverlaps(list.out[[2]][[2]], list.out[[2]][[3]])
        gr.nonCG <- subsetByOverlaps(
            gr.nonCG, list.out[[2]][[1]], invert = TRUE, type = 'any')
        nonCG <- grToDF(reduce(gr.nonCG))
        nonCG$seqnames <- as.integer(as.character(nonCG$seqnames))
        nonCG <- nonCG %>% semi_join(nonCG, by = c("seqnames","start","end"))
        DMR.list.out(
            context.df=nonCG, out.name=paste0(tmp.name, "nonCG-DMRs"),
            data.out=data.dir)
        message("Done!")
        #multi-context
        message("Generating multi-context DMRs...")
        CG.gr <- list.out[[2]][[1]]
        CHG.gr <- list.out[[2]][[2]]
        CHH.gr <- list.out[[2]][[3]]
        multi <- grToDF(intersect(CG.gr, intersect(CHG.gr, CHH.gr)))
        if (nrow(multi) != 0) {
            multi$seqnames <- as.integer(as.character(multi$seqnames))
            multi <- multi %>% semi_join(
                multi, by = c("seqnames","start","end"))
            DMR.list.out(
                context.df=multi,out.name=paste0(
                    tmp.name, "multi-context-DMRs"), data.out=data.dir)
        } else {
            message("No multi-context DMRs found!")
        }
    } else {
        stop("Filtered DMR matrix files for all contexts do not exist!")
    }
}

#------------------------------------------------------------------------------
#' Output context-specific DMRs (only for a given context, nonCG or all context)
#' @param samplelist DataFrame object containing information about
#'                   file, sample, replicate and group. (DataFrame object)
#' @param input.dir Path to the input directory with DMR matrices. (char)
#' @param output.dir Path to the output directory. (char)
#' @param if.filtered Logical to specify if we should use filtered or
#'                    non-filtered matrices. (logical)
#' @importFrom data.table fread
#' @return Output files containing context-specific DMRs for five categories:
#'         CG-only, 
#'         CHG-only, 
#'         CHH-only, 
#'         non-CG (occurs for CHG and/or CHH, but not for CG context), 
#'         multi-context (occurs for CG and CHG and CHH contexts)
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
#' for (cx in c("CG", "CHG", "CHH")) {
#'     mk(paste0(cx, "_StateCalls.txt"), st)
#'     mk(paste0(cx, "_rcMethlvl.txt"), st * 0.9)
#' }
#'
#' sheet <- data.frame(
#'     file = NA_character_, sample = rep(c("WT", "mutant1"), each = 2),
#'     replicate = rep(c("rep1", "rep2"), 2), stringsAsFactors = FALSE)
#'
#' outd <- file.path(tempdir(), "jDMRgrid_context")
#' context.specific.DMRs(samplelist = sheet, input.dir = d,
#'                       output.dir = outd, if.filtered = FALSE)
#' list.files(outd)
#' @export
#' 
context.specific.DMRs <- function(
        samplelist, input.dir, output.dir, if.filtered = FALSE) {
    ## FIX: wrote into a user-supplied directory without creating it, so a fresh run
    ## failed with "No such file or directory".
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(samplelist$group)){
        samplelist$name <- paste0(samplelist$sample,"_", samplelist$replicate)
        gps <- samplelist$group[!samplelist$group %in% c('control')]
        gps <- unique(gps)
        for (m in seq_along(gps)){
            myvec <- c("control", gps[m])
            gp1 <- samplelist$name[which(samplelist$group==myvec[1])]
            gp2 <- samplelist$name[which(samplelist$group==myvec[2])]
            
            gp1.sample <- unique(
                samplelist$sample[which(samplelist$name==gp1)])
            gp2.sample <- unique(
                samplelist$sample[which(samplelist$name==gp2)])
            message("Generating context specific DMRs for ",
                    gp1.sample, "-", gp2.sample,"\n")
            name_string <- ifelse(
                if.filtered==FALSE,"_StateCalls.txt","_StateCalls-filtered.txt")
            CG.f <- paste0(
                input.dir, '/',"CG_", gp1.sample, "_", gp2.sample, name_string)
            CHG.f <- paste0(
                input.dir, '/',"CHG_", gp1.sample, "_", gp2.sample, name_string)
            CHH.f <- paste0(
                input.dir, '/',"CHH_", gp1.sample, "_", gp2.sample, name_string)
            extract.context.DMRs(
                file1=CG.f,file2=CHG.f,file3=CHH.f, tmp.name=paste0(
                    gp1.sample, "_", gp2.sample, "_"), data.dir=output.dir)
        }
    } else {
        message("Generating context specific DMRs. No groups found!")
        name_string <- ifelse(
            if.filtered==FALSE,"_StateCalls.txt","_StateCalls-filtered.txt")
        output<-extract.context.DMRs(
            file1=paste0(input.dir,"/CG", name_string),
            file2=paste0(input.dir,"/CHG", name_string),
            file3=paste0(input.dir,"/CHH", name_string),
            tmp.name="",data.dir=output.dir)
    }
}