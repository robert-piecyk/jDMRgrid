\name{makeDMRmatrix}
\alias{makeDMRmatrix}
\title{Builds a DMR matrix for all samples}
\usage{
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = FALSE, samplefiles,
    input.dir, out.dir, include.intermediate = FALSE)
}
\arguments{
\item{contexts}{Vector of cytosine contexts selected for DMR calling. By 
default this option is set for all 3 cytosine contexts CG, CHG and CHH. 
(character vector; by default c('CG','CHG','CHH'))}

\item{postMax.out}{(A logical specyfing whether or not to output posterior
probability for the identified state call known as postMax probability.
(optional; logical; by default is FALSE))}

\item{samplefiles}{File containing full path of base level methylation calls, 
sample names and replicates(optional; character)}

\item{input.dir}{Path to the input directory containing all region level 
methylome calls (character)}

\item{out.dir}{Path to the output directory (character)}

\item{include.intermediate}{A logical specifying whether or not the 
intermediate component should be included in the HMM model.
(optional; logical; be default is FALSE)}

}
\value{
Output files containing matrices of state calls, methylation levels and 
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