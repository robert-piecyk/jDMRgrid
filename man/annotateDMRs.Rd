\name{annotateDMRs}
\alias{annotateDMRs}
\title{Annotate DMRs}
\usage{
annotateDMRs(annotation, gff.files, if.gff3, input.dir, out.dir)
}
\arguments{
\item{annotation}{annotation terms used to annotate DMRs}

\item{gff.files}{multiple gff3 files can be supplied as a vector}

\item{if.gff3}{a logical specifying whether output annotated files in 
gff3 format}

\item{input.dir}{input directory containing filtered DMR matrix/matrices. 
Ideally any file containing 3 columns i.e (chr, start, stop) can be supplied.}

\item{out.dir}{output directory}
}
\value{
output files containing annotated DMRs and DMR counts table.
}
\description{
This function takes gff3 files as input and outputs annotated DMRs in text 
and gff3 format. Additionally, a DMR count table is generated.
}
\examples{
## Get some toy data
data.dir.pop <- "folder_population/annotate_DMRs"
data.dir.rep <- "folder_replicate/annotate_DMRs"
gff.file_promoters <- system.file(
    "extdata/toyData","TAIR10_promoters.gff3", package="jDMRgrid")
gff.file_promoters <- system.file(
    "extdata/toyData", "TAIR10_TE.gff3", package="jDMRgrid")
gff.file_genes <- system.file(
    "extdata/toyData", "TAIR10.gene.chr1.gff3", package="jDMRgrid")
annotateDMRs(gff.files = c(
    gff.file_promoters, gff.file_TE, gff.file_genes), annotation = c(
        "promoters", "TE", "gene"), 
    input.dir = "folder_population/context_DMRs", if.gff3 = FALSE, 
    out.dir = data.dir.pop)
}
