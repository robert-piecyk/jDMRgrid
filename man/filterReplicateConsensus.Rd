\name{filterReplicateConsensus}
\alias{filterReplicateConsensus}
\title{filter DMR matrix according to the replicate consensus}
\usage{
filterReplicateConsensus(
    status.collect, rc.methlevel.collect, replicate.consensus, diff.ct=0.5)
}

\arguments{

\item{status.collect}{Data frame containing state calls for a given sample.}

\item{rc.methlevel.collect}{Data frame containing methylation levels for 
a given sample.}

\item{replicate.consensus}{Numeric threshold as the percentage of concordance 
in methylation states among samples with multiple replicates. 
Applicable for control/treatment data.}

\item{diff.ct}{I need to discover what it is actually. By default this option 
is set to 0.5.}

}
\value{
List of filtered matrices for state-calls and methylation levels.
}
\description{
Filtering of DMR matrices (state-calls and methylation levels) according to the 
percentage of concordance in methylation states among samples with multiple 
replicates. Applicable for control/treatment data.
}
\examples{
TBA
}