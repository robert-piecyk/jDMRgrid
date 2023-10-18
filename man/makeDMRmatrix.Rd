\name{makeDMRmatrix}
\alias{makeDMRmatrix}
\title{Builds a DMR matrix for all samples}
\usage{
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = FALSE, samplefiles,
    input.dir, out.dir, include.intermediate = FALSE)
}
\arguments{
\item{samplefiles}{file containing full path of base level methylation calls, 
sample names and replicates(optional)}

\item{input.dir}{input directory containing all region level methylome calls}

\item{out.dir}{output directory}

\item{include.intermediate}{A logical specifying whether or not the 
intermediate component should be included in the HMM model.
By default this option is set as FALSE.}

\item{context}{cytosine context}
}
\value{
output files containing matrices of state calls, methylation levels and 
optionally postMax probabilities.
}
\description{
This function generates a binary matrix, a matrix of recalibrated methylation 
levels and posterior probabilities for all samples.
}
\examples{
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = TRUE, 
    samplefiles = system.file(
        "extdata", "listFiles1.fn", package="jDMRgrid"), 
    input.dir = 'folder_population/grid', out.dir = 'folder_population/matrix',
    include.intermediate = FALSE)
}