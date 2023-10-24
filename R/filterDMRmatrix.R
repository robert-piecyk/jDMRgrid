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
    Sys.sleep(1/max.iter)
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
    cb <- combn(df.gp$n,2)
    my.diff <- unlist(lapply(cb, function(x) abs(cb[1]-cb[2])))
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
#' @export
#'
filterReplicateConsensus <- function(
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
#' @export
#' 
filterEpiMAF <- function(status.collect, rc.methlevel.collect, epiMAF){
    floorDec <- function(valParm ,x){
        y <- function(x, level=1) round(x - 5*10^(-level-1), level)
        res <-y(as.numeric(valParm),x)
        return(res)
    }
    pb2 <- txtProgressBar(
        min = 1, max = NROW(status.collect), char = "=", style = 3, file = "")
    mypatterns <- status.collect[, 4:(NCOL(status.collect) - 1)]
    epiMAFs <- apply(mypatterns, 1, function(row)
    {
        mycount <- table(row)
        epiMAF.out <- min(mycount) / length(row)
        return(floorDec(as.numeric(as.character(epiMAF.out)), 5))
    })
    status.collect$epiMAF <- epiMAFs
    
    Sys.sleep(1 / NROW(status.collect))
    setTxtProgressBar(pb2, NROW(status.collect))
    close(pb2)
    
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
        data.reduced <- as.data.frame(data.reduced)[,c(1,2,3)]
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
        if.mergingBins) {
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
        # we are only removing non-polymorphic patterns
        if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
            removeNonPolymorphicOnly(
                status.collect = status.collect, 
                rc.methlevel.collect = rc.methlevel.collect, 
                context = context, data.dir = data.dir)}
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
                if.mergingBins = if.mergingBins, diff.ct = 0.5)}
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
#' @param if.mergingBins Logical if bins with same state-calls profile among
#'                       samples should be merged. (logical; default as FALSE)
#' @import magrittr
#' @importFrom data.table fread
#' @return Filtrated state-calls and methylation DMR matrices using:
#'         filtration for non-polymorphic patterns and 
#'         (epiMAF for population-based data or 
#'         replicate consensus for control/treatment data) 
#'         and/or merging neighbouring bins with same state calls profiles.
#' @export
#'
filterDMRmatrix <- function(
        epiMAF.cutoff = NULL, replicate.consensus = NULL, data.dir, 
        samplelist, if.mergingBins = FALSE) 
{
    contexts <- unique(sub("_.*", "", list.files(data.dir, pattern = 'C')))
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
            epiMAF.cutoff = epiMAF.cutoff, if.mergingBins = if.mergingBins)
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
    out.only <- as.data.frame(out.2)
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
        nonCG <- as.data.frame(reduce(gr.nonCG))
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
        multi <- as.data.frame(intersect(CG.gr, intersect(CHG.gr, CHH.gr)))
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
#' @export
#' 
context.specific.DMRs <- function(
        samplelist, input.dir, output.dir, if.filtered = FALSE) {
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