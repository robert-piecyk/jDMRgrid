\name{filterEpiMAF}
\alias{filterEpiMAF}
\title{filter DMR matrix according to the minor epi-allele frequency}
\usage{
filterEpiMAF(status.collect, rc.methlevel.collect, epiMAF)
}
\arguments{
\item{status.collect}{Data frame containing state calls for a given sample. 
(DataFrame object)}

\item{rc.methlevel.collect}{Data frame containing methylation levels for 
a given sample. (DataFrame object)}

\item{epiMAF}{Numeric threshold to filter for minor epi-allele 
frequency. Applicable for population level data. By default this option is 
set to NULL. (NULL or numeric value between 0 and 1)}
}
\value{
list of data frames for rc-methylation level and state-calls using 
minor epi-allele frequency criterion.
}
\description{
Filter for minor epi-allele frequency in population-based data. 
Example: Let's have 4 regions to be filtered having the following frequency of 
the state calls among all sample (M = Meth, I = Intermediate, U = Unmeth)
Region 1: 1M, 2I, 7U
Region 2: 3M, 7U
Region 3: 2M 4I 4U
Region 4: 5M 5U

Let's calculate the frequency of the state call which occurs least often:
Region 1: Meth (10%)
Region 2: Meth (30%)
Region 3: Meth (20%)
Region 4: Meth (50%)

If the epiMAF cutoff is set to 0.33, region 1 and region 2 would be removed
from the final dataset.
}
\examples{
list.out <- filterEpiMAF(
    status.collect = data.table::fread(
        system.file(
            "extdata", "matrix_population/CG_StateCalls.txt", 
            package="jDMRgrid")), 
    rc.methlevel.collect = data.table::fread(
        system.file(
            "extdata", "matrix_population/CG_rcMethlvl.txt", 
            package="jDMRgrid")), 
    epiMAF = 0.3)
}