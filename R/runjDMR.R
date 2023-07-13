export.bins <- function(mylist, myinfo, out.dir, runName)
{
  lapply(seq_along(myinfo$bin.size), function(z1) {
    window.size.1 <- myinfo$bin.size[z1]
    step.size.1 <- myinfo$step.size[z1]
    mylist.1 <- mylist[[as.character(window.size.1)]]
    context.1 <- myinfo$context[z1]
    out.name.1 <- paste0(out.dir, "/",runName,"_Win", window.size.1,"_Step", step.size.1,"_",context.1,".Rdata",sep="")
    dput(mylist.1, out.name.1)
  })
}
binGenome <- function(methimputefiles,
                      contexts,
                      window,
                      step,
                      min.C,
                      out.dir,
                      runName)
  {
  # message about creating grid
  message('Creating grid...')
  cyt.collect <- list()
  # from one of the methIMPUTE file extract all cytosines positions
  meth.out <- data.table::fread(methimputefiles[1],showProgress=FALSE)
  all.cyt.pos <- meth.out[,c('seqnames','start','strand','context')]
  colnames(all.cyt.pos) <- c('chr','pos','strand','context')
  # create a GRanges object from the cytosines positions
  cyt_gr <- GRanges(seqnames=all.cyt.pos$chr,
                    ranges=IRanges(start=all.cyt.pos$pos, width=1),
                    context=all.cyt.pos$context,
                    strand=all.cyt.pos$strand)
  # get length of chromosomes
  chr.names <- unique(meth.out$seqnames)
  chr.lengths <- unlist(lapply(chr.names, function(x) max(meth.out$start[meth.out$seqnames == x])))
  names(chr.lengths) <- chr.names
  # create a GRanges object from the chromosome (fasta) start-end positions
  gr <- GRanges(seqnames=names(chr.lengths), ranges=IRanges(start=1, end=chr.lengths))
  # main loop function
  if (length(window) != length(step)) {
    message('EXECUTION HALTED! Window and step vectors sizes must have the same length.')
  } else {
    results <- lapply(seq_along(window), function(x1) {
      window.size <- window[x1]
      step.size <- step[x1]
      # Binning genome
      binned.g <- slidingWindows(gr, width = window.size, step = step.size)
      message("Binning genome with windows of: ", window.size, "bp and step-size of: ", step.size, "bp.")
      # Creating a data frame from the binned data
      dd <- data.frame(unlist(binned.g))
      names(dd) <- c("chr", "start", "end", "cluster.length", "strand")
      # Storing the data frame in a list
      new <- list(dd)
      names(new) <- as.numeric(as.character(window.size))
      # Creating a GRange object
      data_gr <- GRanges(seqnames = dd$chr, ranges = IRanges(start = dd$start, end = dd$end), clusterlen = dd$cluster.length)
      # Creating a data frame for bin and step sizes
      mydf <- data.frame(bin.size = window.size, step.size = step.size)
      # Process each context using lapply and combine the results
      results <- lapply(contexts, function(cx) {
        message("Extracting cytosines for ", cx, ".")
        # Filtering cytosines
        ref_gr <- cyt_gr[which(cyt_gr$context == cx),]
        # Counting cytosines in GRanges
        dat.collect <- GenomicRanges::countOverlaps(data_gr, ref_gr)
        # Create a empirical distribution of cytosines within bins and find a threshold based on its min.C percentile
        new.dat.collect <- dat.collect[(which(dat.collect >= as.numeric(quantile(ecdf(dat.collect), min.C/100))))]
        non.empty.bins <- length(new.dat.collect) / length(dat.collect)
        return(non.empty.bins)
      })

      # Add the results to the mydf data frame
      mydf <- cbind(mydf, data.frame(context = contexts,
                                     ratio = unlist(results)))
      return(list(mydf = mydf, collect.bins = new))
    })
    out <- data.table::rbindlist(lapply(results, function(x) x$mydf))
    out <- out[order(out$bin.size, out$step.size),]
    out <- split(out, f = out$context)
    mybins <- out
    data.table::fwrite(data.table::rbindlist(out), file = paste0(out.dir, '/', runName, '_optimal_minC_threshold.csv'))
    #mybins <- lapply(out, function(x) x[which.min(x$ratio),]) #deleted selection of the min ratio per context
    collect.bins <- lapply(results, function(x) x$collect.bins)
    message("Exporting regions...")
    lapply(mybins, function(x) export.bins(mylist=collect.bins[[1]],
                                           myinfo=x,
                                           out.dir=out.dir,
                                           runName=runName))
    return(list.files(out.dir, pattern=paste0(".*", runName, ".*\\.Rdata$"), full.names=TRUE))
  }
  message("Done!")
}
#-----------------------------------------------------------------------------------------
#' Run jDMR on binned genome
#'
#' this function runs a HMM model on a genome binned using a sliding/non-sliding window approach
#' @param out.dir Output directory. (character)
#' @param window Bin size. (numeric vector)
#' @param step Step size. (numeric vector)
#' @param samplefiles Path to the text file containing path to samples and sample names. For control/treatment data an additional column specifying the replicates is required. (character)
#' @param contexts Vector of cytosine contexts selected for DMR calling. By default this option is set for all 3 cytosine contexts CG, CHG and CHH. (character vector)
#' @param min.C Percentile threshold based on empirical distribution of the cytosines across bins. (numeric value between 0 and 100)
#' @param mincov Minimum read coverage over cytosines. By default this option is set as 0. (numeric value between 0 and 1)
#' @param include.intermediate A logical specifying whether or not the intermediate component should be included in the HMM model. By default this option is set as FALSE. (logical)
#' @param runName Name of the operation. By default this option is set to 'GridGenome'. (character)
#' @import magrittr
#' @importFrom data.table fread fwrite
#' @importFrom stringr str_remove_all
#' @importFrom data.table fwrite
#' @importFrom IRanges slidingWindows
#' @export
#'
runjDMRgrid <- function(out.dir,
                        window,
                        step,
                        samplefiles,
                        contexts=c('CG','CHG','CHH'),
                        min.C,
                        mincov=0,
                        include.intermediate=FALSE,
                        runName='GridGenome')
  {
  bin.genome.files <- binGenome(methimputefiles=data.table::fread(samplefiles)$file,
                                contexts=contexts,
                                window=window,
                                step=step,
                                min.C=min.C,
                                out.dir=out.dir,
                                runName=runName)
  bin.select <- lapply(seq_along(contexts), function(x) {
    b <- bin.genome.files[grep(contexts[x], bin.genome.files)]
    if (length(b) != 0) {
      names(b) <- contexts[x]
      return(b)
    }
  })
  bin.select <- Filter(Negate(is.null), bin.select)
  names(bin.select) <- unlist(lapply(bin.select, function(x) names(x)))
  merge.list <- vector(mode="list")
  filelist <- data.table::fread(samplefiles, header=TRUE)

  info_lapply <- lapply(seq_along(bin.select), function(j) {
    refRegion <- dget(bin.select[[j]][[1]])
    refRegion <- list(reg.obs = refRegion)

    info_vapply <- vapply(seq_along(filelist$file), function(i) {
      methfn <- gsub(".*methylome_|\\.txt|_All.txt$", "", filelist$file[i])
      message("Running file: ", methfn, " for context: ", names(bin.select)[j], "\n")
      fileName <- basename(methfn)
      grid.out <- makeMethimpute(
        df = filelist$file[i],
        context = names(bin.select)[j],
        refRegion = refRegion,
        fit.plot = FALSE,
        include.intermediate = include.intermediate,
        probability = "constrained",
        out.dir = out.dir,
        fit.name = paste0(basename(methfn), "_", names(bin.select)[j]),
        name = fileName,
        mincov = mincov
      )
    }, FUN.VALUE=NULL)
  })
}
