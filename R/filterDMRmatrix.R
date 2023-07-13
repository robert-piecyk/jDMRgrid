#'
#' @param status.collect
#' @param rc.methlevel.collect
#' @param replicate.consensus
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom dplyr semi_join summarize
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'

filterReplicateConsensus <- function(status.collect, rc.methlevel.collect, replicate.consensus, diff.ct=0.5){

  if (!is.null(status.collect$epiMAF)){
    status.collect <- status.collect[,-c("epiMAF")]
  }
  # deducing replicates info
  mycol <- names(status.collect)[4:NCOL(status.collect)]
  sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))

  if (length(sampleinfo)==2){
    colnames(sampleinfo) <- c("sample", "replicate")
    dt <- data.frame()
    pb1 <- txtProgressBar(min = 1, max = NROW(status.collect), char = "=", style = 3, file = "")

    q <- lapply(1:NROW(status.collect), function(x){
      mypattern <- unlist(status.collect[x, 4:NCOL(status.collect)])
      df.bind <- cbind(sampleinfo, mypattern)
      for (m in unique(df.bind$sample)){
        total.reps <- length(df.bind$mypattern[df.bind$sample==m])
        rval <- round(replicate.consensus * total.reps)
        pattern.vals <- df.bind$mypattern[df.bind$sample==m]
        df.bind$diff.count[df.bind$sample==m] <- length(which(pattern.vals==1))/total.reps
        tt <- table(pattern.vals)
        if (max(tt) >= rval){
          df.bind$count[df.bind$sample==m] <- 0
        } else {
          df.bind_rows <- apply(df.bind, 1, function(row) {
            m <- row[which(colnames(df.bind) == 'sample')]
            total.reps <- length(df.bind$mypattern[df.bind$sample == m])
            rval <- round(replicate.consensus * total.reps)
            pattern.vals <- df.bind$mypattern[df.bind$sample == m]
            row[["diff.count"]] <- length(which(pattern.vals == 1)) / total.reps
            tt <- table(pattern.vals)
            if (max(tt) >= rval) {
              row[["count"]] <- 0
            } else {
              row[["count"]] <- 1
            }
            return(row)
          })
          df.bind_rows <- as.data.frame(t(df.bind_rows))
          df.bind_rows$mypattern <- as.numeric(df.bind_rows$mypattern)
          df.bind_rows$diff.count <- as.numeric(df.bind_rows$diff.count)
          df.bind_rows$count <- as.numeric(df.bind_rows$count)
          df.bind <- df.bind_rows
        }
      }
      Sys.sleep(1/NROW(status.collect))
      setTxtProgressBar(pb1, x)
      close(pb1)

      df.gp <- group_by(df.bind, sample) %>% dplyr::summarize(n = mean(diff.count))
      cb <- combn(df.gp$n,2)
      my.diff <- unlist(lapply(cb, function(x) abs(cb[1]-cb[2])))

      # allowing 50% difference between control and treatment groups
      if ((min(my.diff) >= diff.ct) & (sum(df.bind$count)==0)) {
        dt <- rbind(dt, status.collect[x,])
      }
      #print(df.bind)
      #if (sum(df.bind$count)==0) {
      #  dt <- rbind(dt, status.collect[x,])
      #}
    })
    out <- q[!sapply(q,is.null)]
    #status.collect <- q[-which(sapply(q, is.null))]
    df.status.collect <- rbindlist(out)
    if (NROW(df.status.collect) !=0){
      df.rc.methlevel.collect <- rc.methlevel.collect %>% dplyr::semi_join(df.status.collect,
                                                                           by=c("seqnames","start","end"))
      return(list(df.status.collect, df.rc.methlevel.collect))
    } else {
      message("\nEmpty dataframe. Nothing to write!")
      return(NULL)
    }
  } else {
    stop("Column for replicates is missing!!!")
  }
}


#------------------------------------------------------------------------------------------------

#' @param mat1
#' @param mat2
#' @param epiMAF
#' @import magrittr
#' @importFrom dplyr semi_join
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'

filterEpiMAF <- function(mat1, mat2, epiMAF){
  floorDec <- function(valParm ,x){
    y <- function(x, level=1) round(x - 5*10^(-level-1), level)
    res <-y(as.numeric(valParm),x)
    return(res)
  }
  pb2 <- txtProgressBar(min = 1, max = NROW(mat1), char = "=", style = 3, file = "")
  mypatterns <- mat1[, 4:(NCOL(mat1) - 1)]
  epiMAFs <- apply(mypatterns, 1, function(row)
  {
    mycount <- table(row)
    epiMAF.out <- min(mycount) / length(row)
    return(floorDec(as.numeric(as.character(epiMAF.out)), 5))
  })
  mat1$epiMAF <- epiMAFs

  Sys.sleep(1 / NROW(mat1))
  setTxtProgressBar(pb2, NROW(mat1))
  close(pb2)

  df.status.collect <- mat1[which(mat1$epiMAF < epiMAF),]
  if (NROW(df.status.collect) !=0){
    df.rc.methlevel.collect <- mat2 %>% dplyr::semi_join(df.status.collect, by=c("seqnames","start","end"))
    return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("Empty dataframe. Nothing to write!")
    return(NULL)
  }
}

#------------------------------------------------------------------------------------------------
#'
#' @param rcmethlvl
#' @param statecalls
#' @param gap
#' @import magrittr
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges findOverlaps
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table rbindlist
#' @export
#'
merge.bins <- function(rcmethlvl, statecalls)
{
  # split state call dataset into unique patterns of state calls among samples
  states.all <- split(statecalls, apply(statecalls[, 4:NCOL(statecalls)], 1, paste, collapse = "_"))
  rcmethlvl.all <- split(rcmethlvl, apply(statecalls[, 4:NCOL(statecalls)], 1, paste, collapse = "_"))

  # merge overlapping bins having the same patterns among samples AND save the indices corresponding the merged windows
  states.list <- lapply(states.all, function(x) {
    data.list <- list()
    # create a GRanges object
    data <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
    # merge overlapping bins
    data.reduced <- GenomicRanges::reduce(data)
    # find the indices of the overlapping bins
    ovlps <- as.data.frame(GenomicRanges::findOverlaps(data, data.reduced))
    # prepare dataset for the final step
    data.reduced <- as.data.frame(data.reduced)
    data.reduced <- data.reduced[,c(1,2,3)]
    data.reduced <- cbind(data.reduced, x[1,4:NCOL(x)])
    data.list[[1]] <- ovlps
    data.list[[2]] <- data.reduced
    return(data.list)
  })

  # extract state calls and overlaps from the previous function
  overlaps.data <- lapply(states.list, function(x) x[[1]])
  states.data <- lapply(states.list, function(x) x[[2]])

  # calculate the average methylation for each previously merged bin by using the overlaps dataset from the previous function
  rcmethlvl.data <- lapply(seq_along(rcmethlvl.all), function(x) {
    # extract methylation levels and overlaps for a given run
    rcmethlvl.1 <- rcmethlvl.all[[x]]
    ovlps.1 <- overlaps.data[[x]]
    ovlps.1 <- split(ovlps.1, f = ovlps.1$subjectHits)
    rcmethlvl.out <- lapply(ovlps.1, function(y) {
      myv <- apply(rcmethlvl.1[y$queryHits,4:NCOL(rcmethlvl)], 2, mean)
      mydf <- data.frame(t(myv), stringsAsFactors = FALSE)
      colnames(mydf) <- names(myv)
      return(mydf)
    })
    # bind mean methylation windows with state calls seqnames, start and end
    rcmethlvl.out <- cbind(states.data[[x]][,1:3], data.table::rbindlist(rcmethlvl.out))
    return(rcmethlvl.out)
  })

  # prepare dataset to get it through next step of analysis
  states.data <- data.table::rbindlist(states.data)
  states.data <- states.data[order(seqnames, start, end), ]
  rcmethlvl.data <- data.table::rbindlist(rcmethlvl.data)
  rcmethlvl.data <- rcmethlvl.data[order(seqnames, start, end), ]
  return(list(states.data, rcmethlvl.data))
}
#------------------------------------------------------------------------------------------------

export.out <- function(out.rcmethlvl, out.statecalls, context, out.name1, out.name2, data.out){
  data.table::fwrite(x=out.statecalls,
         file=paste0(data.out, "/", context, "_", out.name1, ".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  if (!is.null(out.rcmethlvl)) {
    data.table::fwrite(x=out.rcmethlvl,
           file=paste0(data.out, "/", context, "_", out.name2, ".txt"),
           quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  } else {
    message("Generate filtered matrix of recalibrated methylation levels set to FALSE!")
  }
}

#------------------------------------------------------------------------------------------------
#' filter DMR matrix
#'
#' Filters non-polymorphic patterns by default.
#' @param epiMAF.cutoff Numeric threshold to filter for minor epi-allele frequency. Applicable for population level data. By default this option is set to NULL. (NULL or numeric value between 0 and 1)
#' @param replicate.consensus Numeric threshold as the percentage of concordance in methylation states among samples with multiple replicates. Applicable for control/treatment data. By default this option is set to NULL. (NULL or numeric value between 0 and 1)
#' @param data.dir Path to the directory containing DMR matrix files. Looks for files with suffix("_StateCalls.txt" and "_rcMethlvl.txt"). (character)
#' @param samplefiles Path to the text file containing path to samples and sample names. For control/treatment data an additional column specifying the replicates is required. (character)
#' @import magrittr
#' @importFrom data.table fread
#' @export
#'
filterDMRmatrix <- function(epiMAF.cutoff=NULL,
                            replicate.consensus=NULL,
                            data.dir,
                            samplefiles) {

  contexts <- unique(sub("_.*", "", list.files(data.dir, pattern = 'C')))
  ft <- data.table::fread(samplefiles)
  if (!is.null(ft$group)){
    ft$name <- paste0(ft$sample,"_", ft$replicate)
    gps <- ft$group[!ft$group %in% c('control')]
    gps <- unique(gps)

    list.collect2 <- lapply(seq_along(gps), function(m) {
      myvec <- c("control", gps[m])
      gp1 <- ft$name[which(ft$group == myvec[1])]
      gp2 <- ft$name[which(ft$group == myvec[2])]
      gp1.sample <- unique(ft$sample[which(ft$name == gp1)])
      gp2.sample <- unique(ft$sample[which(ft$name == gp2)])
      out.name <- paste0(gp1.sample, "_", gp2.sample)

      list.collect1 <- lapply(seq_along(contexts), function(cn) {
        fn1 <- paste0(data.dir, '/', contexts[cn], "_", out.name, "_StateCalls.txt")
        return(fn1)
      })
      return(list.collect1)
    })
    list.status <- unlist(list.collect2)
  } else {
    list.status <- list.files(data.dir, pattern="_StateCalls.txt", full.names=TRUE)
  }
  if (length(list.status) != 0){
    for (i in seq_along(list.status)){
      # extract context for a given state calls file
      context <- gsub("_StateCalls.txt", "", basename(list.status[i]))
      message("\nFiltering DMR matrix for ", context, ' context.')

      # read a given state calls file
      if (file.exists(list.status[i])){
        status.collect  <- data.table::fread(list.status[i], header=T)
      } else {
        stop("Files do not exist or is non-readable!")
      }

      # remove non-polymorphic/unchanged patterns
      message("Removing non-polymorphic patterns...")

      # extract rows for which state calls for each sample are NOT unchanged (either all 0 or 1)
      index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 &
                       rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)

      # update state calls list with the changed patterns
      status.collect <- status.collect[index,]

      # read a corresponding methylation file and update it using indices from the previous step
      rc.methlvl.name <- paste0(data.dir, "/", context, "_rcMethlvl.txt")
      rc.methlevel.collect <- data.table::fread(rc.methlvl.name, header=T)
      rc.methlevel.collect <- rc.methlevel.collect[index,]

      # if epiMAP and replicate cutoffs are set to NULL, we are only removing non-polymorphic patterns
      if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
        message("Both, epiMAF and replicate consensus set to NULL. Merge corresponding bins")
        out1 <- status.collect
        out2 <- rc.methlevel.collect
        out <- merge.bins(statecalls=out1, rcmethlvl=out2)
        export.out(out.statecalls=out[[1]],
                   out.rcmethlvl=out[[2]],
                   context=context,
                   out.name1="StateCalls-filtered",
                   out.name2="rcMethlvl-filtered",
                   data.out=data.dir)
      }

      # filtering out regions with epiMAF < Minor Allele Frequency
      if (!is.null(epiMAF.cutoff)) {
        message("Filtering for epiMAF: ", epiMAF.cutoff, '...')
        mydf <- filterEpiMAF(mat1=status.collect, mat2=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
        if (!is.null(mydf)){
          # for population data remove the epiMAF column
          out1 <- mydf[[1]][,-c("epiMAF")]
          out2 <- mydf[[2]]
          out <- merge.bins(statecalls=out1, rcmethlvl=out2)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        } else {
          message('Filtering for epiMAF returns NULL dataset. Proceeding with the original dataset...')
          out1 <- status.collect
          out2 <- rc.methlevel.collect
          out <- merge.bins(statecalls=out1, rcmethlvl=out2)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      }
      # retaining samples based on replicate.consensus
      if (!is.null(replicate.consensus)) {
        message("Filtering for replicate consensus...")
        mydf <- filterReplicateConsensus(status.collect, rc.methlevel.collect, replicate.consensus)
        if (!is.null(mydf)){
          out1 <- mydf[[1]]
          out2 <- mydf[[2]]
          out <- merge.bins(statecalls=out1, rcmethlvl=out2)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        } else {
          message('Filtering for replicate consensus returns NULL dataset. Proceeding with the original dataset...')
          out1 <- status.collect
          out2 <- rc.methlevel.collect
          out <- merge.bins(statecalls=out1, rcmethlvl=out2)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      }
    }
  } else {
    message("\nDMR matrix files do not exist!")
  }
}

#------------------------------------------------------------------------------------------------

DMR.list.out <- function(context.df, out.name, data.out){
  data.table::fwrite(x=context.df,
         file=paste0(data.out, "/", out.name,".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

#------------------------------------------------------------------------------------------------
#' @param data.dir
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
#' @importFrom IRanges setTxtProgressBar txtProgressBar
#'
#'

extract.context.DMRs <- function(file1, file2, file3, tmp.name, data.dir){
  if (file.exists(file1) && file.exists(file2) && file.exists(file3)){
    CG.out <- data.table::fread(file1, header=TRUE)
    CHG.out <- data.table::fread(file2, header=TRUE)
    CHH.out <- data.table::fread(file3, header=TRUE)

    CG.gr <- GenomicRanges::GRanges(seqnames=CG.out$seqnames, ranges=IRanges(CG.out$start, end=CG.out$end))
    CHG.gr <- GenomicRanges::GRanges(seqnames=CHG.out$seqnames, ranges=IRanges(CHG.out$start, end=CHG.out$end))
    CHH.gr <- GenomicRanges::GRanges(seqnames=CHH.out$seqnames, ranges=IRanges(CHH.out$start, end=CHH.out$end))

    #CG.only
    message("Generating CG-only DMRs...")
    out1 <- IRanges::subsetByOverlaps(CG.gr, CHG.gr, invert = TRUE)
    CG.1 <- IRanges::subsetByOverlaps(out1, CHH.gr, invert = TRUE)
    CG.1 <- as.data.frame(CG.1)
    CG.1$seqnames <- as.integer(as.character(CG.1$seqnames))
    CG.only <- CG.out %>% dplyr::semi_join(CG.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CG.only,
                 out.name=paste0(tmp.name,"CG-only-DMRs"),
                 data.out=data.dir)
    message("Done!\n")

    #CHG.only
    message("Generating CHG-only DMRs...")
    out2 <- IRanges::subsetByOverlaps(CHG.gr, CG.gr, invert = TRUE)
    CHG.1 <- IRanges::subsetByOverlaps(out2, CHH.gr, invert = TRUE)
    CHG.1 <- as.data.frame(CHG.1)
    CHG.1$seqnames <-as.integer(as.character(CHG.1$seqnames))
    CHG.only <- CHG.out %>% dplyr::semi_join(CHG.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CHG.only,
                 out.name=paste0(tmp.name,"CHG-only-DMRs"),
                 data.out=data.dir)
    message("Done!\n")

    #CHH.only
    message("Generating CHH-only DMRs...")
    out3 <- IRanges::subsetByOverlaps(CHH.gr, CG.gr, invert = TRUE)
    CHH.1 <- IRanges::subsetByOverlaps(out3, CHG.gr, invert = TRUE)
    CHH.1 <- as.data.frame(CHH.1)
    CHH.1$seqnames <-as.integer(as.character(CHH.1$seqnames))
    CHH.only <- CHH.out %>% dplyr::semi_join(CHH.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CHH.only,
                 out.name=paste0(tmp.name, "CHH-only-DMRs"),
                 data.out=data.dir)
    message("Done!\n")

    #non-CG
    message("Generating non-CG DMRs...")
    nonCG.collect <- list()
    overlaps.nonCG <- GenomicRanges::findOverlaps(CHG.gr, CHH.gr, ignore.strand=TRUE)
    overlaps.hits.nonCG <- IRanges::subsetByOverlaps(CHG.gr, CHH.gr)
    mcols(overlaps.hits.nonCG)$DMRs.CHH.coord<- IRanges::CharacterList(split(CHH.gr[subjectHits(overlaps.nonCG)], queryHits(overlaps.nonCG)))
    out.nonCG <- IRanges::subsetByOverlaps(overlaps.hits.nonCG, CG.gr, invert = TRUE)
    nonCG <- data.frame(out.nonCG) %>% dplyr::mutate(DMRs.CHH.coord = strsplit(as.character(DMRs.CHH.coord), ",")) %>% tidyr::unnest(c(DMRs.CHH.coord))
    nonCG.clean <- data.frame(lapply(nonCG, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
    nonCG.clean <- nonCG.clean %>% tidyr::separate(DMRs.CHH.coord, c("CHH.seqnames","CHH.start","CHH.stop"), sep = '([-:])')
    nonCG.clean <- nonCG.clean[-c(4,5)]
    colnames(nonCG.clean)[1] <- "CHG.seqnames"
    colnames(nonCG.clean)[2] <- "CHG.start"
    colnames(nonCG.clean)[3] <- "CHG.stop"

    pb4 <- txtProgressBar(min = 1, max = NROW(nonCG.clean), char = "=", style = 3, file = "")

    for (i1 in 1:NROW(nonCG.clean)){
      myrow <- nonCG.clean[i1,]
      a=GenomicRanges::makeGRangesFromDataFrame(myrow[,c("CHG.seqnames","CHG.start","CHG.stop")])
      b=GenomicRanges::makeGRangesFromDataFrame(myrow[,c("CHH.seqnames","CHH.start","CHH.stop")])
      out <- data.frame(GenomicRanges::intersect(a,b))
      myrow$merged.seqnames <- out$seqnames
      myrow$merged.start <- out$start
      myrow$merged.stop <- out$end
      nonCG.collect[[i1]] <- data.frame(myrow)
      Sys.sleep(1/NROW(nonCG.clean))
      setTxtProgressBar(pb4, i1)
    }
    close(pb4)
    DMR.list.out(context.df=data.table::rbindlist(nonCG.collect),
                 out.name=paste0(tmp.name,"nonCG-DMRs"),
                 data.out=data.dir)
    message("Done!\n")

    #multi-context
    message("Generating multi-context DMRs...")
    multi.context.collect <- list()
    overlaps.1 <- GenomicRanges::findOverlaps(CG.gr, CHG.gr, ignore.strand=TRUE)
    overlaps.hits.1 <- IRanges::subsetByOverlaps(CG.gr, CHG.gr)
    mcols(overlaps.hits.1)$DMRs.CHG.coord <- IRanges::CharacterList(split(CHG.gr[subjectHits(overlaps.1)], queryHits(overlaps.1)))
    overlaps.2 <- GenomicRanges::findOverlaps(overlaps.hits.1, CHH.gr, ignore.strand=TRUE)
    overlaps.hits.2 <- IRanges::subsetByOverlaps(overlaps.hits.1, CHH.gr)
    if (NROW(overlaps.hits.2)!=0){
      mcols(overlaps.hits.2)$DMRs.CHH.coord <- IRanges::CharacterList(split(CHH.gr[subjectHits(overlaps.2)], queryHits(overlaps.2)))
      multi.context.1 <- data.frame(overlaps.hits.2) %>% dplyr::mutate(DMRs.CHG.coord = strsplit(as.character(DMRs.CHG.coord), ",")) %>% tidyr::unnest(c(DMRs.CHG.coord))
      multi.context.1.clean <- data.frame(lapply(multi.context.1, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
      multi.context.1.clean <- multi.context.1.clean %>% tidyr::separate(DMRs.CHG.coord, c("CHG.seqnames","CHG.start","CHG.stop"), sep = '([-:])')
      multi.context.2.clean <- data.frame(multi.context.1.clean) %>% dplyr::mutate(DMRs.CHH.coord = strsplit(as.character(DMRs.CHH.coord), ",")) %>% tidyr::unnest(c(DMRs.CHH.coord))
      multi.context.2.clean <- multi.context.2.clean %>% tidyr::separate(DMRs.CHH.coord, c("CHH.seqnames","CHH.start","CHH.stop"), sep = '([-:])')
      multi.context.2.clean <- multi.context.2.clean[-c(4,5)]
      colnames(multi.context.2.clean)[1] <- "CG.seqnames"
      colnames(multi.context.2.clean)[2] <- "CG.start"
      colnames(multi.context.2.clean)[3] <- "CG.stop"

      a <- GenomicRanges::makeGRangesFromDataFrame(multi.context.2.clean[, c("CG.seqnames", "CG.start", "CG.stop")])
      b <- GenomicRanges::makeGRangesFromDataFrame(multi.context.2.clean[, c("CHG.seqnames", "CHG.start", "CHG.stop")])
      c <- GenomicRanges::makeGRangesFromDataFrame(multi.context.2.clean[, c("CHH.seqnames", "CHH.start", "CHH.stop")])
      out2 <- data.frame(GenomicRanges::intersect(GenomicRanges::intersect(a, b), c))
      multi.context.collect <- data.frame(matrix(NA, nrow(out2), 0))
      multi.context.collect$merged.seqnames <- out2$seqnames
      multi.context.collect$merged.start <- out2$start
      multi.context.collect$merged.stop <- out2$end
      DMR.list.out(context.df=multi.context.collect,
                   out.name=paste0(tmp.name, "multi-context-DMRs"),
                   data.out=data.dir)
    } else {
      message("No multi-context DMRs found!")
    }
    message("Done!\n")
  } else {
    stop("Filtered DMR matrix files for all contexts donot exist!")
  }
}

context.specific.DMRs <- function(samplefiles, input.dir, output.dir){
  ft <- fread(samplefiles)
  if (!is.null(ft$group)){
    ft$name <- paste0(ft$sample,"_", ft$replicate)
    gps <- ft$group[!ft$group %in% c('control')]
    gps <- unique(gps)
    for (m in seq_along(gps)){
      myvec <- c("control", gps[m])
      gp1 <- ft$name[which(ft$group==myvec[1])]
      gp2 <- ft$name[which(ft$group==myvec[2])]

      gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
      gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
      message("Generating context specific DMRs for ", gp1.sample, "-", gp2.sample,"\n")
      CG.f <- paste0(input.dir, '/',"CG_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      CHG.f <- paste0(input.dir, '/',"CHG_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      CHH.f <- paste0(input.dir, '/',"CHH_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      extract.context.DMRs(file1=CG.f,
                           file2=CHG.f,
                           file3=CHH.f,
                           tmp.name=paste0(gp1.sample, "_", gp2.sample, "_"),
                           data.dir=output.dir)
    }
  } else {
    message("Generating context specific DMRs. No groups found!\n")
    output <- extract.context.DMRs(file1=paste0(input.dir,"/CG_StateCalls-filtered.txt"),
                                   file2=paste0(input.dir,"/CHG_StateCalls-filtered.txt"),
                                   file3=paste0(input.dir,"/CHH_StateCalls-filtered.txt"),
                                   tmp.name="",
                                   data.dir=output.dir)
  }
}



