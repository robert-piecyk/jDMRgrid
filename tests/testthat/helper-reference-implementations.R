## Reference implementations, preserved verbatim from the code that preceded the vectorised
## versions in 0.3.0. They are the definition of correct behaviour for those rewrites: the
## originals iterated rows and regions and were unusably slow at scale, but their arithmetic
## is what the fast paths must reproduce.
##
## Only the progress bar and Sys.sleep() have been removed from ref_filterEpiMAF; those are
## presentation and only slow the reference down.

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

ref_filterReplicateConsensus <- function(
        status.collect, rc.methlevel.collect, replicate.consensus, diff.ct) {
    if (!is.null(status.collect$epiMAF)){
        status.collect <- as.data.frame(status.collect)
        status.collect <- status.collect[,-nrow(status.collect)]
    }
    # deducing replicates info
    mycol <- names(status.collect)[4:NCOL(status.collect)]
    sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))
    if (length(sampleinfo)==2){
        colnames(sampleinfo) <- c("sample", "replicate")
        dt <- data.frame()
        pb1 <- txtProgressBar(
            min = 1, max = nrow(status.collect), char = "=", 
            style = 3, file = "")
        q <- lapply(
            seq(NROW(status.collect)), function(iter) {
                mypattern <- unlist(
                    status.collect[iter, 4:NCOL(status.collect)])
                df.bind <- cbind(sampleinfo, mypattern)
                out.q <- bindReplicateLists(
                    x = iter, mypattern = mypattern, df.bind = df.bind,
                    sampleinfo = sampleinfo, 
                    replicate.consensus = replicate.consensus, 
                    diff.ct = diff.ct, max.iter = nrow(status.collect))
                setTxtProgressBar(pb1, iter)
                return(out.q)
            })
        close(pb1)
        df.status.collect <- status.collect[unlist(q),]
        if (NROW(df.status.collect) !=0){
            df.rc.methlevel.collect <- rc.methlevel.collect %>%
                semi_join(
                    df.status.collect,
                    by=c("seqnames","start","end"))
            return(list(df.status.collect, df.rc.methlevel.collect))
        } else {
            message("Empty dataframe. Nothing to write!")
            return(NULL)
        }
    } else {
        stop("Column for replicates is missing!!!")
    }
}

ref_merge.bins <- function(rcmethlvl_data, statecalls_data)
{
    states.all <- split(statecalls_data, apply(
        statecalls_data[,4:NCOL(statecalls_data)],1,paste,collapse = "_"))
    rcmethlvl.all <- split(
        rcmethlvl_data, apply(
            statecalls_data[,4:NCOL(statecalls_data)],1,paste,collapse = "_"))
    # merge overlapping bins having the same patterns among samples
    states.list <- lapply(states.all, function(x) {
        data.list <- list()
        # create a GRanges object and merge overlapping bins
        data <- makeGRangesFromDataFrame(x,keep.extra.columns = TRUE)
        # merge overlapping bins
        data.reduced <- reduce(data)
        # find the indices of the overlapping bins
        ovlps <- as.data.frame(findOverlaps(data, data.reduced))
        ## grToDF(), not as.data.frame(): on Bioc devel GRanges extends List, so
        ## as.data.frame() lands in as.data.frame.List and errors. grToDF() is proven
        ## identical to as.data.frame() for GRanges, so the reference stays faithful.
        data.reduced <- jDMRgrid:::grToDF(data.reduced)[,c(1,2,3)]
        data.reduced <- cbind(data.reduced, x[1,4:NCOL(x)])
        data.list[[1]] <- ovlps; data.list[[2]] <- data.reduced
        return(data.list)
    })
    # extract state calls and overlaps from the previous function
    overlaps.data <- lapply(states.list, function(x) x[[1]])
    states.data <- lapply(states.list, function(x) x[[2]])
    # calculate the average methylation for each previously merged bin
    rcmethlvl.data <- lapply(seq_along(rcmethlvl.all), function(x) {
        # extract methylation levels and overlaps for a given run
        rcmethlvl.1 <- rcmethlvl.all[[x]]
        ovlps.1 <- overlaps.data[[x]]
        ovlps.1 <- split(ovlps.1, f = ovlps.1$subjectHits)
        rcmethlvl.out <- lapply(ovlps.1, function(y) {
            myv <- apply(
                rcmethlvl.1[y$queryHits,4:NCOL(rcmethlvl_data)], 2, mean)
            mydf <- data.frame(t(myv), stringsAsFactors = FALSE)
            colnames(mydf) <- names(myv)
            return(mydf)
        })
        # bind mean methylation windows with state calls seqnames, start and end
        rcmethlvl.out <- cbind(
            states.data[[x]][,seq_len(3)],
            rbindlist(rcmethlvl.out))
        return(rcmethlvl.out)
    })
    states.data <- rbindlist(states.data)
    states.data <- states.data[order(seqnames, start, end), ]
    rcmethlvl.data <- rbindlist(rcmethlvl.data)
    rcmethlvl.data <- rcmethlvl.data[order(seqnames, start, end), ]
    return(list(states.data, rcmethlvl.data))
}

ref_filterEpiMAF <- function(mat1, mat2, epiMAF){
    floorDec <- function(valParm ,x){
        y <- function(x, level=1) round(x - 5*10^(-level-1), level)
        res <-y(as.numeric(valParm),x)
        return(res)
    }
    mypatterns <- mat1[, 4:(NCOL(mat1) - 1)]
    epiMAFs <- apply(mypatterns, 1, function(row)
    {
        mycount <- table(row)
        epiMAF.out <- min(mycount) / length(row)
        return(floorDec(as.numeric(as.character(epiMAF.out)), 5))
    })
    mat1$epiMAF <- epiMAFs
    
    
    df.status.collect <- mat1[which(mat1$epiMAF < epiMAF),]
    if (NROW(df.status.collect) !=0){
        df.rc.methlevel.collect <- mat2 %>% semi_join(
            df.status.collect, by=c("seqnames","start","end"))
        return(list(df.status.collect, df.rc.methlevel.collect))
    } else {
        message("Empty dataframe. Nothing to write!")
        return(NULL)
    }
}
