#--------------------------------------------------------------------------
#'
#' @param contexts Methylation context(s) to be analysed. (vector)
#' @param c1 Index of the first context. (num)
#' @param c2 Index of the second context. (num)
#' @param cor.array Correlation array gained from the distanceCorrelation.
#' @param dfs
#' @param maxweights
#' @param param.list
#' @param miny
#' @importFrom minpack.lm nlsLM
#' @importFrom stats na.omit coefficients
#' @return
#'
calculateDecayingConstance <- function(
        contexts, c1, c2, cor.array, dfs, maxweights, params.list, miny, skip) 
{
    context.transition <- paste0(contexts[c1], '-', contexts[c2])
    if (c1 <= c2) {
        df <- data.frame(
            distance = as.numeric(dimnames(cor.array)[[3]]),
            correlation = cor.array[c1,c2,,'correlation'],
            weight = cor.array[c1,c2,,'weight'],
            from = contexts[c1], to = contexts[c2])
        y <- df$correlation[(skip + 1):nrow(df)]
        x <- df$distance[(skip + 1):nrow(df)]
        weight <- df$weight[(skip + 1 ):nrow(df)]
        startvalues <- list(a0 = na.omit(y)[1], D = 50)
        p <- NULL
        if (is.null(p)) {
            startvalues <- list(a0 = na.omit(y)[1])
            p <- tryCatch({
                fit <- nlsLM(
                    y ~ a0 * exp(-x/Inf), start=startvalues, 
                    weights=weight)
                s <- summary(fit)
                c <- coefficients(s)
                params <- c[seq_along(startvalues)]
                names(params) <- names(startvalues)
                params <- as.list(params)
                params$D <- Inf
                params
            }, error = function(e) {
                startvalues$D <- Inf
                startvalues
            })
        }
        if (p$D <= 0) {
            p$D <- Inf }
        params.list[[context.transition]] <- p
        df$correlation.fit <- p$a0 * exp(-df$distance/p$D)
        df$logweight <- log(df$weight+1)
        dfs[[context.transition]] <- df
        maxweights[context.transition] <- max(
            df$logweight, na.rm = TRUE)
    }
    return(list(maxweights = maxweights, dfs = dfs, params.list = params.list))
}

#--------------------------------------------------------------------------
#' Obtain an estimate for the transDist parameter by fitting an exponential 
#' function to the supplied correlations from distanceCorrelation. 
#' Modified function from the methIMPUTE package.
#' @inheritParams calculateDecayingConstance
#' @param distcor The output produced by distanceCorrelation.
#' @param skip Skip the first n cytosines for the fitting.
#' @param plot.parameters Whether to plot fitted params on to the plot or not.
#' @import ggplot2
#' @import magrittr
#' @importFrom minpack.lm nlsLM
#' @importFrom grDevices dev.off pdf
#' @importFrom data.table fread fwrite
#' @importFrom stringr str_replace_all
#' @importFrom dplyr bind_rows filter
#' @importFrom stats na.omit coefficients
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @return
#' 
modified.estimateTransDist <- function(distcor, skip, plot.parameters=TRUE) {
    contexts <- dimnames(distcor$data)[[1]]
    cor.array <- distcor$data
    params.list <- list(); dfs <- list(); maxweights <- numeric()
    miny <- min(cor.array, na.rm = TRUE)
    for (c1 in seq_along(contexts)) {
        for (c2 in seq_along(contexts)) {
            output.list <- calculateDecayingConstance(
                contexts, c1, c2, cor.array, dfs, maxweights, params.list, 
                miny, skip)
        }
    }
    maxweight <- max(output.list$maxweights, na.rm = TRUE)
    df <- do.call(rbind, output.list$dfs) #plot correlations
    params.list <- output.list$params.list
    df$a0 <- vapply(
        params.list[paste0(df$from, '-', df$to)], 
        function(x) round(x$a0, 2), numeric(1))
    df$D <- vapply(
        params.list[paste0(df$from, '-', df$to)], 
        function(x) round(x$D, 0), numeric(1))
    df$params <- paste0("a0 = ", df$a0, ", D = ", df$D)
    ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(
        x='distance', y='correlation', alpha='logweight'))
    ggplt <- ggplt + geom_line(aes_string(
        x='distance', y='correlation.fit'), col='blue')
    if (plot.parameters) {
        if (any(!is.na(df$correlation))) {
            max_corr <- max(df$correlation, na.rm = TRUE)
            max_dist <- max(df$distance, na.rm = TRUE)
            
            ggplt <- ggplt + geom_text(
                aes_string(label='params'), x = max_dist, y = max_corr, 
                vjust = 1, hjust = 1)
        } else {
            message("No non-missing values in df$correlation.")
        }
    }
    ggplt <- ggplt + xlab('distance in [bp]')
    ggplt <- ggplt + facet_grid(from ~ to)
    if (miny < 0) {
        ggplt <- ggplt + geom_hline(aes_string(
            'yintercept'=0), linetype=2, alpha=0.5)
    }
    transDist <- vapply(params.list, function(x) x$D, numeric(1))
    return(list(transDist=transDist, plot=ggplt))
}

#--------------------------------------------------------------------------
#'
#' @param model
#' @param out.dir
#' @param context
#' @param name
#' @import magrittr
#' @importFrom data.table fwrite fread
#' @return 
#' 
modifiedExportMethylome <- function(model, out.dir, context, name) {
    #data <- model$data
    data <- model
    final_dataset <- as(data, 'data.frame')
    final_dataset <- final_dataset[,c(
        'seqnames','start','end','strand','context','counts.methylated',
        'counts.total','posteriorMax','status','rc.meth.lvl')]
    
    # dropping columns
    drops <- c(
        'width','strand','clusterlen','counts.methylated','counts.total',
        'distance', 'transitionContext', 'posteriorMeth','posteriorUnmeth')
    final_dataset <- final_dataset[ , !(names(final_dataset) %in% drops)]
    #------------------------------------------------------------------
    # convert full string into M/U/I
    final_dataset <- statusStringCheck(final_dataset)
    #------------------------------------------------------------------
    # take 4 digit of decimal value posteriorMax column
    final_dataset$posteriorMax <-floorDec(
        as.numeric(as.character(final_dataset$posteriorMax)),5)
    final_dataset$rc.meth.lvl <- floorDec(
        as.numeric(as.character(final_dataset$rc.meth.lvl)),5)
    final_dataset$seqnames <- as.character(final_dataset$seqnames)
    
    saveFile <- paste0(out.dir, "/", basename(name), "_", context, ".txt")
    
    fwrite(
        final_dataset, file = saveFile, quote = FALSE, sep = '\t', 
        row.names = FALSE, col.names = TRUE)
    return (final_dataset)
}

#--------------------------------------------------------------------------
#' How many cytosines are methylated and unmethylated in methimpute file
#' @param overlaps
#' @param overlaps.hits
#' @param data_gr
#' @importFrom stats aggregate
#' @importFrom dplyr bind_rows
#' @return
#' 
findMethylatedOrUnmethylated <- function(overlaps, overlaps.hits, data_gr)
{
    counts <- array(
        NA, dim=c(length(data_gr), 2), 
        dimnames=list(NULL, c("methylated", "total")))
    methylated <- aggregate(
        overlaps.hits$methylated, list(subjectHits(overlaps)), FUN=sum)
    total <- aggregate(
        overlaps.hits$total, list(subjectHits(overlaps)), FUN=sum)
    if (NROW(methylated) != NROW(counts) ){
        missingr <- which(
            !rownames(data.frame(data_gr)) %in% methylated$Group.1)
        methylated <- bind_rows(
            data.frame(Group.1=missingr, x=0), methylated)
        methylated <- methylated[order(methylated$Group.1),]
        total <- bind_rows(
            data.frame(Group.1=missingr,x=0), total)
        total <- total[order(total$Group.1),]
    }
    counts[,"methylated"] <- methylated$x
    counts[,"total"] <- total$x
    data_gr$counts <- counts
    return(data_gr)
}

#--------------------------------------------------------------------------
#' Transform methimpute dataset with single cytosine according to the bins
#' @param df
#' @param context
#' @param refRegion
#' @param mincov
#' @param nCytosines
#' @param if.Bismark
#' @param FASTA.file
#' @import magrittr
#' @importFrom dplyr filter
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges countOverlaps findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats aggregate
#' @importFrom dplyr bind_rows
#' @importFrom methimpute importBismark extractCytosinesFromFASTA 
#' @importFrom mehtimpute inflateMethylome
#' @return
#' 
makeRegionsImpute <- function(
        df, context, refRegion, mincov, nCytosines, if.Bismark, FASTA.file) {
    tmp_reg <- refRegion
    data <- as.data.frame(tmp_reg$reg.obs)
    data <- data %>% filter(data$chr != "M" & data$chr != "C")
    #reference methimpute file
    if (if.Bismark == TRUE) {
        ref_data <- importBismark(df)
        cytosine.positions <- extractCytosinesFromFASTA(
            FASTA.file, contexts = as.character(context))
        ref_data <- inflateMethylome(ref_data, cytosine.positions)
        ref_data <- as.data.frame(ref_data)
        ref_data <- ref_data[,-c(3,4)]
        colnames(ref_data) <- c('V1','V2','V3','V4','V5','V6')
    } else {
        ref_data <- fread(
            df, skip = 1, select = c("V1","V2","V3","V4","V5","V6"))
        ref_data$V4 <- as.character(ref_data$V4)
    }
    #remove Mt and chloroplast coordinates. Following is for Arabidopsis only
    ref_data <- ref_data %>% filter(ref_data$V1 != "M" & ref_data$V1 != "C")
    #filtering by coverage
    if (mincov > 0) {
        message("Filtering for coverage: ", mincov)
        ref_data <- ref_data[which(ref_data$V6 >= mincov),] }
    ref_data <- ref_data[which(ref_data$V4 == context),]
    #create GRanges objects for the methimpute and binned genome files
    data_gr <- GRanges(
        seqnames=data$chr, ranges=IRanges(start=data$start, end=data$end),
        clusterlen=data$cluster.length, context=as.factor(context))
    ref_gr <- GRanges(
        seqnames=ref_data$V1, ranges=IRanges(start=ref_data$V2, width=1),
        context=as.factor(context), methylated=ref_data$V5, total=ref_data$V6)
    #count number of overlapping cytosines within bin
    data_gr$cytosineCount <- countOverlaps(data_gr, ref_gr)
    #find overlaps between binned genome and methimpute dataset
    overlaps <- findOverlaps(ref_gr, data_gr)
    overlaps.hits <- ref_gr[queryHits(overlaps)]
    if (NROW(overlaps.hits) != 0){
        data_gr <- findMethylatedOrUnmethylated(
            overlaps = overlaps, overlaps.hits = overlaps.hits, 
            data_gr = data_gr)
    }
    rm(ref_data, ref_gr, overlaps, overlaps.hits)
    return(data_gr)
}

#--------------------------------------------------------------------------
#' Perform HMM for state calling within binned genome
#' @inheritParams modified.estimateTransDist
#' @inheritParams modifiedExportMethylome
#' @param df DataFrame with binned genome. ()
#' @param context Methylation contexts. (char)
#' @param fit.plot
#' @param fit.name
#' @param refRegion
#' @param include.intermediate
#' @param probability
#' @param out.dir
#' @param name
#' @param mincov
#' @param if.Bismark
#' @param FASTA.file
#' @importFrom methimpute distanceCorrelation callMethylation
#' @return 
#' 
makeMethimpute <- function(
        df, context, fit.plot, fit.name, refRegion, include.intermediate, 
        probability, out.dir, name, mincov, if.Bismark, 
        FASTA.file)
    {
    methylome.data <- makeRegionsImpute(
        df=df, context=context, refRegion=refRegion, mincov=mincov, 
        if.Bismark=if.Bismark, FASTA.file=FASTA.file)
    if (!is.null(methylome.data$counts)) {
        quant.cutoff <- as.numeric(
            quantile(
                methylome.data$counts[,"total"], probs = c(0.96), na.rm=TRUE))
        distcor <- distanceCorrelation(data=methylome.data, distances=0:100)
        fit <- modified.estimateTransDist(distcor=distcor, skip=2)
        if (fit.plot==TRUE){
            message("Generating fit plot for ", name)
            pdf(paste0(out.dir, "/", fit.name, "-fit.pdf", sep = ""))
            message(fit)
            dev.off()
        }
        model <- callMethylation(
            data = methylome.data, transDist = fit$transDist,
            count.cutoff = quant.cutoff, max.time = Inf, max.iter = Inf,
            include.intermediate = include.intermediate, update = probability)
        methFile <- modifiedExportMethylome(
            model = model$data, out.dir = out.dir, context = context, 
            name = name)
        rm(model)
    }
    rm(methylome.data)
}