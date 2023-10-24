\name{filterDMRmatrix}
\alias{filterDMRmatrix}
\title{filter DMR matrix}
\usage{
filterDMRmatrix(
    epiMAF.cutoff = NULL, replicate.consensus = NULL, data.dir, samplelist,
    if.mergingBins = FALSE)
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

\item{samplelist}{DataFrame object containing information about
file, sample and replicate. For control/treatment data an 
additional column specifying the replicates is required. (DataFrame object)}

\item{if.mergingBins}{Logical argument if merging consecutive bins 
having the same stateCalls should be performed. 
By default it is set to FALSE. (logical)}
}
\value{
output files containing filtered matrices of state calls, methylation levels 
and optionally postMax probabilities.
}
\description{
Filters non-polymorphic patterns by default. User can decide whether:
-> filters non-polymorphic patterns only (epiMAF.cutoff = NULL, 
replicate.consensus = NULL, if.mergingBins = FALSE)
-> Set if.mergingBins = TRUE only if the number of samples is relatively small
-> Set epiMAF.cutoff to numeric value if you deal with population-based data
-> Set replicate.consensus to numeric value if you deal with replicate/groups

}
\examples{
#------------------------------------------------------------------------------
# REPLICATE-BASED DATASET; STEP IV
#------------------------------------------------------------------------------
# Step 1.7: Copy required files and create required directories for outputs
load(system.file("data", "listFiles2.RData", package="jDMRgrid"))
listFiles2$file <- system.file("extdata", listFiles2$file,  package="jDMRgrid")
input_file <- system.file("extdata", "matrix_CT_replicate", package="jDMRgrid")
system2(
    command = "/bin/cp", 
    args = paste0("-r ", input_file, " ", "~/"))

# Step 1.8: filterDMRmatrix for replicate-based data
filterDMRmatrix(
    epiMAF.cutoff = NULL, replicate.consensus = 0.5, 
    data.dir = "~/matrix_CT_replicate", 
    samplelist = listFiles2, if.mergingBins = TRUE)

#------------------------------------------------------------------------------
# POPULATION-BASED DATASET; STEP III
#------------------------------------------------------------------------------
# Step 2.5: Copy required files and create required directories for outputs
load(system.file("data", "listFiles1.RData", package="jDMRgrid"))
listFiles1$file <- system.file("extdata", listFiles1$file,  package="jDMRgrid")
input_file <- system.file("extdata", "matrix_population", package="jDMRgrid")
system2(
    command = "/bin/cp", 
    args = paste0("-r ", input_file, " ", "~/"))

# Step 2.6: filterDMRmatrix for population-based data
filterDMRmatrix(
    epiMAF.cutoff = 0.33, replicate.consensus = NULL, 
    data.dir = "~/matrix_population", 
    samplelist = listFiles1, if.mergingBins = TRUE)

}