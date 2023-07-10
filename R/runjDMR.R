
extractCytclusters <- function(in.fasta, contexts, genome, out.dir){
  all.cyt.pos=CfromFASTAv5(fasta=in.fasta)
  val <- c()
  f <- seqinr::read.fasta(in.fasta)
  for (x in seq_along(f)){
    val[x] <- gsub(">", "", seqinr::getAnnot(f[[x]]))
  }
  mychrs <- val
  for  (chrs in seq_along(mychrs)){
    ref.genome <- all.cyt.pos[which(all.cyt.pos$chr==mychrs[chrs]),]
    makeReg(ref.genome=ref.genome,
            contexts=contexts, # all contexts can be specified as ("CG","CHG","CHH","C")
            makeRegnull=c(TRUE,TRUE,TRUE), # can be set to FALSE if null distribution is not required to be generated
            chr=mychrs[chrs],
            N.boot=10^5,
            N.sim.C="all",
            fp.rate=0.01,
            set.tol=0.01,
            out.dir,
            min.C=5,
            genome)
  }
  out.files <- list.files(out.dir, pattern=paste0(genome,"_regions_"), full.names = TRUE)
  return(out.files)
}

#-----------------------------------------------------------------------------------------
#' Run jDMR for cytosine clusters
#'
#' This function runs a HMM model on identified Cytosine clusters
#' @param out.dir output directory
#' @param fasta.file path to genome fasta file
#' @param samplefiles a text file containing path to samples and sample names. For control/treatment data an additional column specifying the replicates is required.
#' @param genome Genome name as a string .e.g Arabidopsis or Human etc
#' @param contexts cytosine contexts as a vector. By default this option is set for all 3 cytosine contexts CG, CHG and CHH.
#' @param include.intermediate A logical specifying whether or not the intermediate component should be included in the HMM model.By default it is set as FALSE.
#' @param mincov Minimum read coverage over cytosines. Default is set as 0.
#' @param nCytosines Minimum number of cytosines. Default is set as 0.
#' @importFrom data.table fread
#' @import GenomicRanges
#' @importFrom stringr str_remove_all
#' @export
#'
runjDMRregions <- function(out.dir,
                           fasta.file,
                           samplefiles,
                           genome,
                           contexts=c('CG','CHG','CHH'),
                           include.intermediate=FALSE,
                           mincov=0,
                           nCytosines=0) {
  Regionfiles <- extractCytclusters(in.fasta=fasta.file, contexts, genome, out.dir=out.dir)
  df.obs <- list()
  df.sim <- list()
  merge.list <- vector(mode="list")

  # Read the sample file with filenames and file paths
  filelist <- data.table::fread(samplefiles, header=TRUE)

  for (j in seq_along(contexts)){
    Regfiles <- Regionfiles[grep(paste0(contexts[j], ".Rdata"),Regionfiles)]
    cat(paste0("Reading Region files. Merging individual chr data for context ", contexts[j], " ...\n"), sep="")
    cat("\n")
    for (k1 in seq_along(Regfiles)){
      f.file <- dget(Regfiles[k1])
      if (NROW(f.file$reg.obs)==0) {
        cat(paste0("Empty file ", basename(Regfiles[k1]), " ...\n"), sep="")
      } else {
        df.obs[[k1]] <- as.data.frame(f.file$reg.obs)
        df.sim[[k1]] <- as.data.frame(f.file$reg.sim)
      }
    }
    refRegion <- list(reg.obs=do.call(rbind,df.obs),
                      reg.sim=do.call(rbind,df.sim),
                      context=contexts[j])
    #print(refRegion)

    for (k2 in seq_along(filelist$file)){
      if (file.exists(filelist$file[k2])==TRUE){
        methfn <- gsub(".*methylome_|\\.txt|_All.txt$", "", filelist$file[k2])
        cat(paste0("Now running file: ", methfn, " for context ", contexts[j], " ...\n"), sep="")
        regions.out <- makeMethimpute(
          df=filelist$file[k2],
          context=contexts[j],
          refRegion=refRegion,
          fit.plot=FALSE,
          include.intermediate=include.intermediate,
          probability="constrained",
          out.dir=out.dir,
          fit.name=paste0(methfn, "_", contexts[j]),
          name=methfn,
          nCytosines=nCytosines,
          mincov=mincov)
      } else {
        stop("Check filenames! ", filelist$file[k2], " doesnot exist.\n")
      }
    }
  }
}

#-----------------------------------------------------------------------------------------
export.bins <- function(mylist, bin.context, mybin, step, out.dir, genome){
  for (i in 1:length(bin.context))
  {
    out.name <- paste0(out.dir, "/", genome,"_Win", mybin, "_Step", step, "_", bin.context[i], ".Rdata", sep="")
    dput(mylist, out.name)
  }
}

binGenome <- function(fasta.file,
                      contexts,
                      win,
                      min.C,
                      out.dir,
                      genome,
                      selectProperBinSize = FALSE)
  {
  cyt.collect <- list()
  # all cytosine positions
  all.cyt.pos <- CfromFASTAv5(fasta=fasta.file)
  cyt_gr <- GRanges(seqnames=all.cyt.pos$chr,
                    ranges=IRanges(start=all.cyt.pos$pos, width=1),
                    context=all.cyt.pos$context,
                    strand=all.cyt.pos$strand)
  
  # get length of chromosomes
  val <- c()
  f <- seqinr::read.fasta(fasta.file)
  for (x in seq_along(f)){
    val[x] <- seqinr::getLength(f[[x]])
    names(val)[x] <- gsub(">", "", seqinr::getAnnot(f[[x]]))
  }
  gr <- GRanges(seqnames=names(val), ranges=IRanges(start=1, end=val))
  
  mydf.collect <- list()
  collect.bins <- list()
  
  mybins <- c()
  if (nrow(win)>1){
    cat("No bin size specified. Automatically determining bin size!\n")
    for (x1 in 1:nrow(win)){
      cat("----------------------------------------------------------------", "\n")
      window.size <- win$window[x1]
      step.size <- win$step[x1]
      binned.g <- slidingWindows(gr, width = window.size, step = step.size)
      cat(paste0("Binning genome with windows of: ", window.size, 
                 " bp and step-size of: ", step.size, " bp\n"), sep = "")
      dd <- data.frame(unlist(binned.g))
      names(dd)[1] <- "chr"
      names(dd)[2] <- "start"
      names(dd)[3] <- "end"
      names(dd)[4] <- "cluster.length"
      
      new <- list(dd)
      names(new) <- as.numeric(as.character(window.size))
      #names(new) <- "reg.obs"
      collect.bins <- append(collect.bins, new)
      #dput(new, out.name)
      
      # tmp_reg <- new
      # data <- as.data.frame(tmp_reg$reg.obs)
      data_gr <- GRanges(seqnames=dd$chr,
                         ranges=IRanges(start=dd$start, end=dd$end),
                         clusterlen=dd$cluster.length)
      
      mydf <- data.frame(bin.size=window.size,
                         step.size=step.size)
      
      for (cx in seq_along(contexts)){
        cat("Extracting cytosines for ", contexts[cx], "\n")
        ref_gr <- cyt_gr[which(cyt_gr$context==contexts[cx]),]
        data_gr$cytosineCount <- GenomicRanges::countOverlaps(data_gr, ref_gr)
        dat.collect <- data_gr$cytosineCount
        
        new.dat.collect <- length(which(dat.collect>=min.C))
        non.empty.bins <- new.dat.collect/length(dat.collect)
        mydf[,ncol(mydf) + 1] <- non.empty.bins
        colnames(mydf)[ncol(mydf)] <- paste0(contexts[cx])
        rm(dat.collect, new.dat.collect)
      }
      mydf.collect[[x1]] <- mydf
    }
    out <- rbindlist(mydf.collect)
    out <- out[order(out$bin.size, out$step.size),]
    if (selectProperBinSize == TRUE)
    {
      list.min <- list()
      if ('CG' %in% contexts) {
        CG.min <- out[which(out$CG >= 0.9),]
        CG.min <- CG.min[which.min(CG.min$CG),c(1,2)]
        list.min[['CG']] <- CG.min
      } 
      if ('CHG' %in% contexts) {
        CHG.min <- out[which(out$CHG >= 0.9),]
        CHG.min <- mybins[which.min(mybins$CHG),c(1,2)]
        list.min[['CHG']] <- CG.min
      } 
      if ('CHH' %in% contexts) {
        CHH.min <- out[which(out$CHH >= 0.9),]
        CHH.min <- mybins[which.min(mybins$CHH),c(1,2)]
        list.min[['CHH']] <- CHH.min
      }
      mybins <- rbindlist(list.min)
      mybins$context <- names(list.min)
      cat("----------------------------------------------------------------", "\n")
      cat("Exporting regions ....","\n")
      for (xx in 1:nrow(mybins))
      {
        export.df <- collect.bins[which(names(collect.bins)==mybins$bin.size[xx])]
        export.bins(mylist=export.df[[1]], 
                    bin.context=mybins$context[xx], 
                    mybin=mybins$bin.size[xx], 
                    step=mybins$step.size[xx], 
                    out.dir, 
                    genome)
      }
      return(list.files(out.dir, pattern=paste0(".*",genome,".*\\.Rdata$"), full.names=TRUE))
      cat("done!", "\n")      
    } else {
      mybins <- out[rowSums((out[,3:ncol(out)] >= 0.9) == FALSE) == 0,]
      cat("----------------------------------------------------------------", "\n")
      cat("Exporting regions ....","\n")
      for (xx in 1:nrow(mybins))
      {
        export.df <- collect.bins[which(names(collect.bins)==mybins$bin.size[[xx]])]
        export.bins(mylist=export.df[[1]], 
                    bin.context=colnames(mybins)[-c(1,2)], 
                    mybin=mybins$bin.size[xx], 
                    step=mybins$step.size[xx], 
                    out.dir, 
                    genome)
      }
      return(list.files(out.dir, pattern=paste0(".*",genome,".*\\.Rdata$"), full.names=TRUE))
      cat("done!", "\n")
    }
  } else {
    cat("----------------------------------------------------------------", "\n")
    # User input bin size
    window.size <- win$window
    step.size <- win$step     
    binned.g <- slidingWindows(gr, width = window.size, step = step.size)
    cat(paste0("Binning genome with windows of: ", window.size, " bp and step-size of: ", step.size, " bp\n"), sep = "")
    
    dd <- data.frame(unlist(binned.g))
    names(dd)[1] <- "chr"
    names(dd)[2] <- "start"
    names(dd)[3] <- "end"
    names(dd)[4] <- "cluster.length"
    
    new <- list(dd)
    names(new) <- as.numeric(as.character(window.size))
    #names(new) <- "reg.obs"
    collect.bins <- append(collect.bins, new)
    data_gr <- GRanges(seqnames=dd$chr,
                       ranges=IRanges(start=dd$start, end=dd$end),
                       clusterlen=dd$cluster.length)
    mydf <- data.frame(bin.size=window.size,
                       step.size=step.size)
    for (cx in seq_along(contexts)){
      cat("Extracting cytosines for ", contexts[cx], "\n")
      ref_gr <- cyt_gr[which(cyt_gr$context==contexts[cx]),]
      data_gr$cytosineCount <- GenomicRanges::countOverlaps(data_gr, ref_gr)
      dat.collect <- data_gr$cytosineCount
      
      new.dat.collect <- length(which(dat.collect>=min.C))
      non.empty.bins <- new.dat.collect/length(dat.collect)
      mydf[,ncol(mydf) + 1] <- non.empty.bins
      colnames(mydf)[ncol(mydf)] <- paste0(contexts[cx])
      rm(dat.collect, new.dat.collect)
    }
    
    cat("----------------------------------------------------------------", "\n")
    cat("Exporting regions ....","\n")
    export.bins(mylist=collect.bins[[1]], 
                bin.context=contexts, 
                mybin=win$window, 
                step=win$step, 
                out.dir, 
                genome)
    return(list.files(out.dir, pattern=paste0(".*",genome,".*\\.Rdata$"), full.names=TRUE))
    cat("done!", "\n")
    }
}
  # cat("----------------------------------------------------------------", "\n")
  # cat("Exporting regions ....","\n")
  # for (xx in seq_along(mybins)){
  #   export.df <- collect.bins[which(names(collect.bins)==mybins[[xx]])]
  #   export.bins(mylist=export.df, mybin=mybins[[xx]], out.dir, genome)
  # }
  # return(list(CG=min.CG,
  #             CHG=min.CHG,
  #             CHH=min.CHH))
  # cat("done!", "\n")


#-----------------------------------------------------------------------------------------
#' Run jDMR on binned genome
#'
#' this function runs a HMM model on a genome binned using a sliding/non-sliding window approach
#' @param out.dir Output directory
#' @param fasta.file Path to genome fasta file
#' @param samplefiles A text file containing path to samples and sample names. For control/treatment data an additional column specifying the replicates is required.
#' @param win A data frame consisiting of two columns (window and step), specyfing the parameters to bin a genome. If you want to use a sliding window, please use the same step size and bin size.
#' @param min.C Minimum number of cytosines in at least 90 percent of the bins/regions.
#' @param genome Genome name as a string .e.g Arabidopsis or Human etc.
#' @param contexts cytosine contexts as a vector. By default this option is set for all 3 cytosine contexts CG, CHG and CHH.
#' @param mincov Minimum read coverage over cytosines. Default is set as 0.
#' @param include.intermediate A logical specifying whether or not the intermediate component should be included in the HMM model. By default it is set as FALSE.
#' @param selectProperBinSize Logical specifying wheter we should select the proper bin size for each context among all possibilities in win. By default it is set as FALSE.
#' @param nCytosines Minimum number of cytosines. Default is set as 0.
#' @importFrom data.table fread
#' @importFrom stringr str_remove_all
#' @export
#'
runjDMRgrid <- function(out.dir,
                        fasta.file,
                        samplefiles,
                        win=data.frame(window=c(100,200,300,400,500,600,700,800,900,1000),
                                       step=c(100,200,300,400,500,600,700,800,900,1000)/2),
                        genome,
                        contexts=c('CG','CHG','CHH'),
                        min.C=0,
                        mincov=0,
                        include.intermediate=FALSE,
                        selectProperBinSize=FALSE,
                        nCytosines=0){

  #bin.genome.files <- unlist(lapply(contexts, function(g) list.files(out.dir, paste0(g,".Rdata"), full.names=TRUE)))
  bin.genome.files <- list.files(out.dir, pattern=paste0(".*",genome,".*\\.Rdata$"), full.names=TRUE)
  if (length(bin.genome.files) == 0) {
    cat("binned genome files donot exist. Creating for the first time...")
    cat("\n")
    bin.genome.files <- binGenome(fasta.file,
                                  win=win,
                                  contexts=contexts,
                                  out.dir=out.dir,
                                  min.C=min.C,
                                  genome=genome,
                                  selectProperBinSize = selectProperBinSize)
  } else {
    cat(paste0("\nbinned genome files ", bin.genome.files, " exist!\n"))
    cat("\n")
  }
  bin.select <- list()
  for (x in 1:length(contexts)){
    b <- bin.genome.files[grep(contexts[x], bin.genome.files)]
    if ((length(b)) != 0) {
      bin.select[[x]] <- b
      names(bin.select)[[x]] <- contexts[x]
    }
  }
  merge.list <- vector(mode="list")
  filelist <- data.table::fread(samplefiles, header=TRUE)
  for (j in seq_along(bin.select)){
    #cat(paste0("Win ", bin.select[[j]], " Step ", bin.select[[j]], " context ", names(bin.select)[j],  "\n"))
    #cat(bin.select[[j]])
    # refRegion.name <- paste0(out.dir, "/", genome, "_Win", bin.select[[j]],
    #                          "_Step", bin.select[[j]], "_", names(bin.select)[j], ".Rdata", sep="")

    refRegion <- dget(bin.select[[j]][[1]])
    refRegion <- list(reg.obs = refRegion)
    for (i in seq_along(filelist$file)){
      methfn <- gsub(".*methylome_|\\.txt|_All.txt$", "", filelist$file[i])
      cat("----------------------------------------------------------------", "\n")
      cat(paste0("Running file: ", methfn," for context: ", names(bin.select)[j],"\n"), sep = "")
      grid.out <- makeMethimpute(
        df=filelist$file[i],
        context=names(bin.select)[j],
        refRegion=refRegion,
        fit.plot=FALSE,
        include.intermediate=include.intermediate,
        probability="constrained",
        out.dir=out.dir,
        fit.name=paste0(methfn, "_", names(bin.select)[j]),
        name=methfn,
        nCytosines=nCytosines,
        mincov=mincov)
    }
  }
}
