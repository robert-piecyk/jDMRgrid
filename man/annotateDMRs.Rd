\name{annotateDMRs}
\alias{annotateDMRs}
\title{Annotate DMRs}
\usage{
annotateDMRs(annotation, gff.files, if.gff3, input.dir, out.dir)
}
\arguments{
\item{annotation}{Vector with the annotation terms used to annotate DMRs. 
(character)}

\item{gff.files}{Vector with the multiple gff3 files corresponding to 
the annotation terms. (character)}

\item{if.gff3}{Logical specifying whether output annotated files in 
gff3 format. (logical)}

\item{input.dir}{Path to the input directory containing filtered DMR matrices. 
Ideally any file containing 3 columns i.e (chr, start, stop) can be supplied.
(character)}

\item{out.dir}{Path to the output directory. (character)}
}

\value{
output files containing annotated DMRs and DMR counts table.
}
\description{
This function takes gff3 files as input and outputs annotated DMRs in text 
and gff3 format. Additionally, a DMR count table is generated.
}
\examples{
system2(command = "/bin/mkdir", args = "~/annotation.files")
annotateDMRs(
    annotation = c('promoters','TE'), 
    gff.files = c(
        system.file("extdata/toyData","TAIR10_promoters.gff3", 
                    package="jDMRgrid"),
        system.file("extdata/toyData","TAIR10_TE.gff3", 
                    package="jDMRgrid")), 
    if.gff3 = TRUE, 
    input.dir = system.file(
        "extdata", "matrix_CT_replicate", package="jDMRgrid"),
    out.dir = "~/annotation.files")
}
