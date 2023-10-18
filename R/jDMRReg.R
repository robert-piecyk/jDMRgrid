#'
#' @param model
#' @param out.dir
#' @param name
#' @param distcor
#' @param skip
#' @param plot.parameters
#' @param df
#' @param refRegion
#' @param context Cytosine context
#' @param fit.plot
#' @param fit.name
#' @param refRegion
#' @param include.intermediate
#' @param probability
#' @param mincov
#' @param nCytosines
#' @import ggplot2
#' @import magrittr
#' @importFrom minpack.lm nlsLM
#' @importFrom grDevices dev.off pdf
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom stringr str_replace_all
#' @importFrom dplyr bind_rows filter
#' @importFrom stats na.omit coefficients
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methimpute callMethylation
#' @importFrom methimpute distanceCorrelation
#' @importFrom stats coefficients
#' @export

modified.estimateTransDist <- function(distcor, skip=2, plot.parameters=TRUE) {
    
    ## Context correlation fits and plots
    contexts <- dimnames(distcor$data)[[1]]
    cor.array <- distcor$data
    maxweights <- numeric()
    params.list <- list()
    miny <- min(cor.array, na.rm = TRUE)
    dfs <- list()
    for (c1 in seq_along(contexts)) {
        for (c2 in seq_along(contexts)) {
            context.transition <- paste0(contexts[c1], '-', contexts[c2])
            if (distcor$separate.contexts) {
                if (c1 != c2) {
                    next
                }
            }
            if (c1 <= c2) {
                df <- data.frame(
                    distance = as.numeric(dimnames(cor.array)[[3]]),
                    correlation = cor.array[c1,c2,,'correlation'],
                    weight = cor.array[c1,c2,,'weight'],
                    from = contexts[c1], to = contexts[c2])
                
                ## Fit
                y <- df$correlation[(skip+1):nrow(df)]
                x <- df$distance[(skip+1):nrow(df)]
                weight <- df$weight[(skip+1):nrow(df)]
                startvalues <- list(a0 = stats::na.omit(y)[1], D = 50)
                p <- NULL
                
                if (is.null(p)) {
                    startvalues <- list(a0 = stats::na.omit(y)[1])
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
                
                ## Check if we have negative D
                if (p$D <= 0) {
                    p$D <- Inf
                }
                params.list[[context.transition]] <- p
                
                ## Plot
                df$correlation.fit <- p$a0 * exp(-df$distance/p$D)
                df$logweight <- log(df$weight+1)
                dfs[[context.transition]] <- df
                maxweights[context.transition] <- max(
                    df$logweight, na.rm = TRUE)
            }
        }
    }
    maxweight <- max(maxweights, na.rm = TRUE)
    
    ## Plot correlation
    df <- do.call(rbind, dfs)
    #df$a0 <- round(
    #    sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'a0'), 2)
    #df$D <- round(
    #    sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'D'), 0)
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
    
    #transDist <- sapply(params.list, '[[', 'D')
    transDist <- vapply(params.list, function(x) x$D, numeric(1))
    return(list(transDist=transDist, plot=ggplt))
}

#--------------------------------------------------------------------------
#'
#' @import magrittr
#' @importFrom data.table fwrite fread
#' @export
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
#'
#' @import magrittr
#' @importFrom dplyr filter
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges countOverlaps findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats aggregate
#' @importFrom dplyr bind_rows
makeRegionsImpute <- function(df, context, refRegion, mincov, nCytosines) {
    tmp_reg <- refRegion
    data <- as.data.frame(tmp_reg$reg.obs)
    data <- data %>% filter(data$chr != "M" & data$chr != "C")
    
    #reference methimpute file
    ref_data <- fread(df, skip = 1, select = c("V1","V2","V3","V4","V5","V6"))
    
    #remove Mt and chloroplast coordinates. Following is for Arabidopsis only
    ref_data <- ref_data %>% filter(ref_data$V1 != "M" & ref_data$V1 != "C")
    
    #filtering by coverage
    if (mincov>0){
        message("Filtering for coverage: ", mincov)
        ref_data <- ref_data[which(ref_data$V6 >= mincov),]
    }
    
    ref_data <- ref_data[which(ref_data$V4==context),]
    
    data_gr <- GRanges(
        seqnames=data$chr, ranges=IRanges(start=data$start, end=data$end),
        clusterlen=data$cluster.length, context=as.factor(context))
    
    ref_gr <- GRanges(
        seqnames=ref_data$V1, ranges=IRanges(start=ref_data$V2, width=1),
        context=as.factor(context), methylated=ref_data$V5, total=ref_data$V6)
    
    data_gr$cytosineCount <- countOverlaps(data_gr, ref_gr)
    
    counts <- array(
        NA, dim=c(length(data_gr), 2), 
        dimnames=list(NULL, c("methylated", "total")))
    
    overlaps <- findOverlaps(ref_gr, data_gr)
    
    overlaps.hits <- ref_gr[queryHits(overlaps)]
    if (NROW(overlaps.hits) != 0){
        
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
    }
    rm(ref_data, ref_gr, overlaps, overlaps.hits)
    return(data_gr)
}

#--------------------------------------------------------------------------
makeMethimpute <- function(
        df, context, fit.plot, fit.name, refRegion, include.intermediate, 
        probability, out.dir, name, mincov)
    {
    methylome.data <- makeRegionsImpute(df, context, refRegion, mincov)
    if (!is.null(methylome.data$counts)) {
        quant.cutoff <- as.numeric(
            quantile(
                methylome.data$counts[,"total"], probs = c(0.96), na.rm=TRUE))
        distcor <- distanceCorrelation(methylome.data, distances=0:100)
        fit <- modified.estimateTransDist(distcor)
        
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