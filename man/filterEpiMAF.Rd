\name{filterEpiMAF}
\alias{filterEpiMAF}
\title{filter DMR matrix according to the minor epi-allele frequency}
\usage{
filterEpiMAF(
    epiMAF.cutoff = NULL, replicate.consensus = NULL, data.dir, samplefiles,
    if.mergingBins = TRUE)
}
\arguments{
\item{epiMAF.cutoff}{Numeric threshold to filter for minor epi-allele 
frequency. Applicable for population level data. By default this option is 
set to NULL. (NULL or numeric value between 0 and 1)}

\item{replicate.consensus}{Numeric threshold as the percentage of concordance 
in methylation states among samples with multiple replicates. 
Applicable for control/treatment data. By default this option is set to NULL. 
(NULL or numeric value between 0 and 1)}

\item{data.dir}{Path to the directory containing DMR matrix files. 
Looks for files with suffix("_StateCalls.txt" and "_rcMethlvl.txt"). 
(character)}

\item{samplefiles}{Path to the text file containing path to samples 
and sample names. For control/treatment data an 
additional column specifying the replicates is required. (character)}

\item{if.mergingBins}{Logical argument if merging consecutive bins 
having the same stateCalls should be performed. 
By default it is set to TRUE. (logical)}
}
\value{
output files containing filtered matrices of state calls, methylation levels 
and optionally postMax probabilities.
}
\description{
Filters non-polymorphic patterns by default.
}
\examples{
TBA
}