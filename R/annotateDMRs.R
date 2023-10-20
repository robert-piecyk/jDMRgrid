#-----------------------------------------------------------------------------
#' Import gff files
#' @param gff Paths to gff files
#' @importFrom  rtracklayer import.gff3
#' @return merge all supplied gff3 annotations into one
#' 
gff3.in <- function(gff){
    input.gff <- lapply(gff, function(x){
        import.gff3(x, colnames=c("type", "ID"))
    })
    merged.gff <- do.call(c, input.gff)
    return(merged.gff)
}

#-----------------------------------------------------------------------------
#' "chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#' "three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#' "snoRNA","snRNA","rRNA","TE","promoters"
#' Output annotated gff3 files
#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @importFrom GenomicRanges findOverlaps
#' @importFrom  rtracklayer export.gff
#' @importFrom  S4Vectors values values<- elementMetadata<-
#' @return export output files in gff3 format
#' 
gff3.out <- function(annotation, grangesObj, gff, name, out.dir) {
    getgff3 <- lapply(annotation, function(x){
        idx <- which(elementMetadata(gff)[,"type"] == x)
        gff <- gff[idx,]
        hits <- findOverlaps(grangesObj, gff, ignore.strand=FALSE)
        gr.matched <- grangesObj[queryHits(hits)]
        if (NROW(gr.matched) != 0){
            mcols(gr.matched) <- cbind.data.frame(
                mcols(gr.matched),mcols(gff[subjectHits(hits)]))
            values(gr.matched) <- cbind(values(gr.matched), region="DMR")
            names(elementMetadata(gr.matched))[names(
                elementMetadata(gr.matched)) == "type"] <- "annotation"
        }
        return(gr.matched)
    })
    export.gff(do.call(c, getgff3), paste0(
        out.dir,"/", name, "_annotation.gff3"), version="3")
}

#-----------------------------------------------------------------------------
#' Extract annotated regions
#' @inheritParams gff3.in
#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @param getAnno
#' @param mygff
#' @param mygr
#' @param annotation
#' @param gff.files
#' @param gff3.out
#' @param input.dir
#' @import magrittr
#' @importFrom dplyr mutate id
#' @importFrom GenomicRanges findOverlaps subsetByOverlaps seqnames
#' @importFrom S4Vectors elementMetadata mcols  mcols<-
#' @importFrom IRanges CharacterList ranges
#' @importFrom rtracklayer strand
#' @importFrom Biostrings type
#' @importFrom methods as
#' @importFrom rlang .data
#' @importFrom tidyr unnest
#' @return annotated list
#' 
annotate <- function(getAnno, mygff, mygr){
    lapply(getAnno, function(x){
        idx <- which(elementMetadata(mygff)[,"type"] == x)
        mygff_subset <- mygff[idx,]
        hits <- findOverlaps(mygr, mygff_subset, ignore.strand=FALSE)
        myranges <- subsetByOverlaps(mygr, mygff_subset)
        
        mcols(myranges)$id <- CharacterList(
            split(mygff_subset$ID[subjectHits(hits)],queryHits(hits)))
        mcols(myranges)$type <- CharacterList(
            split(mygff_subset$type[subjectHits(hits)],queryHits(hits)))
        mcols(myranges)$chr <- CharacterList(
            split(seqnames(mygff_subset)[subjectHits(hits)],queryHits(hits)))
        mcols(myranges)$coord <- CharacterList(
            split(ranges(mygff_subset)[subjectHits(hits)],queryHits(hits)))
        mcols(myranges)$anno.str <- CharacterList(
            split(strand(mygff_subset)[subjectHits(hits)],queryHits(hits)))
        df <- as(myranges, "data.frame")
        df <- df %>% mutate(id = strsplit(
            as.character(id),","), type = strsplit(
                as.character(type),","), chr = strsplit(
                    as.character(chr),","), coord = strsplit(
                        as.character(coord),","), anno.str = strsplit(
                            as.character(anno.str), ",")) %>%
            unnest(
                c(id,type,chr,coord,anno.str))
        cleandf <- data.frame(
            lapply(df, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
        cleandf$Annotation.coord <- apply(
            cleandf[,c("chr","coord","anno.str")], 1, paste, collapse=":")
        cleandf <- subset(cleandf, select = -c(
            chr, coord, anno.str))
        return(cleandf)
    })
}

#-----------------------------------------------------------------------------
#' ff
#' @param d
#' @param annotation Annotation terms used to annotate DMRs
#' @param tmp.name
#' @param out.dir
#' @import magrittr
#' @importFrom  data.table fwrite uniqueN
#' @importFrom  dplyr group_by summarize
#' @importFrom  stats start end
#' 
writeTxtAnnotations <- function(
        d, gr, annotation, tmp.name, out.dir, anno.list, final.df)
{
    out <- d %>%
        group_by(seqnames, start, end) %>%
        summarize(
            type = paste(type, collapse=","),
            id = paste(id, collapse=","),
        ) %>% as.data.frame()
    if (NROW(out)!=0){
        for (k1 in seq_len(NROW(out))){
            out$unique.anno.type[k1] <- paste0(unlist(lapply(
                strsplit(out$type[k1], ","), unique)), collapse=",")
        }
    }
    # count the DMR overlaps; the output
    out.1 <- out[which(unlist(
        lapply(strsplit(out$type,','), uniqueN))==1),]
    out.2 <- out[which(unlist(
        lapply(strsplit(out$type,','), uniqueN))>=2),]
    #counting unique annotations
    for (k2 in seq_along(annotation)){
        anno.list[[k2]] <- NROW(
            out.1[grep(annotation[k2], out.1$type),])
        names(anno.list)[[k2]] <- annotation[k2]
    }
    #also counting multiple overlaps
    anno.list[length(annotation)+1] <- list(NROW(out.2))
    names(anno.list)[[length(annotation)+1]] <- "multiple.overlaps"
    df.1 <- do.call(cbind, anno.list)
    df <- cbind(sample=tmp.name, total.DMRs=NROW(gr), df.1)
    final.df <- rbind(final.df, data.frame(df))
    out <- out[,-c(4,6)]
    fwrite(
        x=out, file=paste0(out.dir, "/", tmp.name, "_annotation.txt"), 
        quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    fwrite(
        x=final.df, file=paste0(out.dir, "/DMR-counts.txt"), quote=FALSE, 
        row.names=FALSE, col.names=TRUE, sep="\t")
}

#-----------------------------------------------------------------------------
#' Output annotated DMRs in text and gff3 format, and DMR count table.
#' @inheritParams gff3.in
#' @inheritParams gff3.out
#' @inheritParams annotate
#' @param out.dir Output directory
#' @param annotation Annotation terms used to annotate DMRs
#' @param gff.files Multiple gff3 files can be supplied as a vector
#' @param gff3.out A logical specifying whether output annotated files in gff3
#' @param input.dir Input directory containing filtered DMR matrix/matrices
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread rbindlist
#' @export
#' @return output files containing annotated DMRs and DMR counts table
#' 
annotateDMRs <- function(annotation, gff.files, if.gff3, input.dir, out.dir) {
    anno.list <- list(); final.df <- list()
    file.list <- list.files(input.dir, pattern="*.txt", full.names = TRUE)
    for (i in seq_along(file.list)){
        message("Running file ", file.list[i])
        tmp.name <- gsub("\\.txt$", "", basename(file.list[i]))
        file <- fread(file.list[i], select=c(1,2,3))
        if (NROW(file) != 0)
        {
            colnames(file) <- c('V1','V2','V3')
            gr <- GRanges(
                seqnames=file$V1,ranges=IRanges(start=file$V2, end=file$V3))
            gff <- gff3.in(gff.files)
            if (if.gff3==TRUE) {
                gff3.out(
                    annotation = annotation, gff = gff, grangesObj = gr,
                    name = tmp.name, out.dir = out.dir)
            }
            d <- rbindlist(annotate(getAnno=annotation, mygff=gff, mygr=gr))
            if (nrow(d) != 0) {
                writeTxtAnnotations(
                    d = d, gr = gr, annotation = annotation, 
                    tmp.name = tmp.name, out.dir = out.dir, 
                    anno.list = anno.list, final.df = final.df)
            }
        }
    }
}
