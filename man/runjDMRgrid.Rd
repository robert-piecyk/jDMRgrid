\name{runjDMRgrid}
\alias{runjDMRgrid}
\title{Run jDMR on binned genome}
\usage{
runjDMRgrid(
    out.dir, window, step, samplelist, contexts = c("CG", "CHG", "CHH"),
    min.C, mincov = 0, include.intermediate = FALSE, runName = "GridGenome",
    numCores = NULL, parallelApply = FALSE, if.Bismark = FALSE, 
    FASTA.file = NULL
)
}
\arguments{
\item{out.dir}{Output directory. (character)}

\item{window}{Bin size. (numeric vector)}

\item{step}{Step size. (numeric vector)}

\item{samplelist}{DataFrame object containing information about
file, sample and replicate. For control/treatment data an 
additional column specifying the replicates is required. (DataFrame object)}

\item{contexts}{Vector of cytosine contexts selected for DMR calling. By 
default this option is set for all 3 cytosine contexts CG, CHG and CHH. 
(character vector)}

\item{min.C}{Percentile threshold based on empirical distribution of the 
cytosines across bins. (numeric value between 0 and 100)}

\item{mincov}{Minimum read coverage over cytosines. 
By default this option is set as 0. (numeric value between 0 and 1)}

\item{include.intermediate}{Logical specifying whether or not the 
intermediate component should be included in the HMM model. 
By default this option is set as FALSE. (logical)}

\item{runName}{Name of the operation. By default this option is 
set to 'GridGenome'. (character)}

\item{numCores}{Number of cores to perform parallel operation using 
foreach loop. Default to NULL. (numeric)}

\item{parallelApply}{Logical specifying if futureapply method should be 
performed to use parallel operation. Default to FALSE. (logical)}

\item{if.Bismark}{Logical if Bismark inputs (CX reports in txt format)
are used. Default as FALSE. (logical)}

\item{FASTA.file}{Path to the FASTA file; required if Bismark outputs are
used. Please use the same FASTA as in Bismark analysis. 
Default as NULL. (char)}
}
\value{
output files containing matrices of state calls, methylation levels 
and optionally postMax probabilities.
}
\description{
this function runs a HMM model on a genome binned using a 
sliding/non-sliding window approach
}
\examples{
#------------------------------------------------------------------------------
# REPLICATE-BASED DATASET; STEP I
#------------------------------------------------------------------------------
# Step 1.1: Copy required files and create required directories for outputs
load(system.file("data", "listFiles2.RData", package="jDMRgrid"))
listFiles2$file <- system.file("extdata", listFiles2$file,  package="jDMRgrid")

system2(command = "/bin/mkdir", args = "~/folder_replicate")
system2(command = "/bin/mkdir", args = "~/folder_replicate/grid-single")
system2(command = "/bin/mkdir", args = "~/folder_replicate/grid-foreach")
system2(command = "/bin/mkdir", args = "~/folder_replicate/grid-future")

# Step 1.2.1: runjDMRgrid in single mode
runjDMRgrid(out.dir = '~/folder_replicate/grid-single',
            window = 200, step = 50, samplelist = listFiles2,
            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
            include.intermediate = TRUE, runName = "Arabidopsis")
## Step 1.2.2: runjDMRgrid in for-each mode
#runjDMRgrid(out.dir = '~/folder_replicate/grid-foreach',
#            window = 200, step = 50, samplelist = listFiles2,,
#            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
#            include.intermediate = TRUE, runName = "Arabidopsis", numCores = 4)
# Step 1.2.3: runjDMRgrid in future-apply mode
runjDMRgrid(out.dir = '~/folder_replicate/grid-future',
            window = 200, step = 50, samplelist = listFiles2,,
            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
            include.intermediate = TRUE, runName = "Arabidopsis",
            parallelApply = TRUE)

#------------------------------------------------------------------------------
# POPULATION-BASED DATASET; STEP I
#------------------------------------------------------------------------------
# Step 2.1: Copy required files and create required directories for outputs
load(system.file("data", "listFiles1.RData", package="jDMRgrid"))
listFiles1$file <- system.file("extdata", listFiles1$file,  package="jDMRgrid")

system2(command = "/bin/mkdir", args = "~/folder_population")
system2(command = "/bin/mkdir", args = "~/folder_population/grid-single")
system2(command = "/bin/mkdir", args = "~/folder_population/grid-foreach")
system2(command = "/bin/mkdir", args = "~/folder_population/grid-future")

# Step 2.2.1: runjDMRgrid in single mode
runjDMRgrid(out.dir = '~/folder_population/grid-single',
            window = 200, step = 50, samplelist = listFiles1,
            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
            include.intermediate = TRUE, runName = "Arabidopsis")
## Step 2.2.2: runjDMRgrid in for-each mode
#runjDMRgrid(out.dir = '~/folder_population/grid-foreach',
#            window = 200, step = 50, samplelist = listFiles1,
#            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
#            include.intermediate = TRUE, runName = "Arabidopsis", numCores = 4)
# Step 2.2.3: runjDMRgrid in future-apply mode
runjDMRgrid(out.dir = '~/folder_population/grid-future',
            window = 200, step = 50, samplelist = listFiles1,
            contexts = c("CG", "CHG", "CHH"), min.C = 10, mincov = 0,
            include.intermediate = TRUE, runName = "Arabidopsis",
            parallelApply = TRUE)
}