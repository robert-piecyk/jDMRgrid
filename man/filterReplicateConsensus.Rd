\name{filterReplicateConsensus}
\alias{filterReplicateConsensus}
\title{filter DMR matrix according to the replicate consensus}
\usage{
filterReplicateConsensus(
    status.collect, rc.methlevel.collect, replicate.consensus, diff.ct)
}

\arguments{

\item{status.collect}{Data frame containing state calls for a given sample.
(DataFrame object)}

\item{rc.methlevel.collect}{Data frame containing methylation levels for 
a given sample. (DataFrame object)}

\item{replicate.consensus}{Numeric threshold as the percentage of concordance 
in methylation states among samples with multiple replicates. 
Applicable for control/treatment data. (numeric value between 0 and 1)}

\item{diff.ct}{Allowed difference between control and treatment groups. 
In filterDMRmatrix function it is set to 0.5 (numeric value between 0 and 1)}

}
\value{
list of data frames for rc-methylation level and state-calls using 
replicate consensus criterion.}
\description{
Filtering of DMR matrices (state-calls and methylation levels) according to the 
percentage of concordance in methylation states among samples with multiple 
replicates. Applicable for control/treatment data.
}
\examples{
list.out <- filterReplicateConsensus(
    status.collect = data.table::fread(
        system.file(
            "extdata", "matrix_CT_replicate/CG_WT_mutant1_StateCalls.txt", 
            package="jDMRgrid")), 
    rc.methlevel.collect = data.table::fread(
        system.file(
            "extdata", "matrix_CT_replicate/CG_WT_mutant1_rcMethlvl.txt", 
            package="jDMRgrid")), 
    replicate.consensus = 0.5,
    diff.ct = 0.5)

}