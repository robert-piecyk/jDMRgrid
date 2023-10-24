\name{makeDMRmatrix}
\alias{makeDMRmatrix}
\title{Builds a DMR matrix for all samples}
\usage{
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = FALSE, samplelist,
    input.dir, out.dir, include.intermediate = FALSE)
}
\arguments{
\item{contexts}{Vector of cytosine contexts selected for DMR calling. By 
default this option is set for all 3 cytosine contexts CG, CHG and CHH. 
(character vector; by default c('CG','CHG','CHH'))}

\item{postMax.out}{(A logical specyfing whether or not to output posterior
probability for the identified state call known as postMax probability.
(optional; logical; by default is FALSE))}

\item{samplelist}{DataFrame object containing information about
file, sample and replicate. For control/treatment data an 
additional column specifying the replicates is required. (DataFrame object)}

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
#------------------------------------------------------------------------------
# REPLICATE-BASED DATASET; STEP II
#------------------------------------------------------------------------------
# Step 1.3: Copy required files and create required directories for outputs
load(system.file("data", "listFiles2.RData", package="jDMRgrid"))
listFiles2$file <- system.file("extdata", listFiles2$file,  package="jDMRgrid")
input_file <- system.file("extdata", "grid_replicate", package="jDMRgrid")
system2(
    command = "/bin/cp", 
    args = paste0("-r ", input_file, " ", "~/"))
system2(command = "/bin/mkdir", args = "~/matrix_replicate")

# Step 1.4: makeDMRmatrix for replicate-based data
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = TRUE, 
    samplelist = listFiles2, input.dir = '~/grid_replicate', 
    out.dir = '~/matrix_replicate', include.intermediate = FALSE)

#------------------------------------------------------------------------------
# POPULATION-BASED DATASET; STEP II
#------------------------------------------------------------------------------
# Step 1.3: Copy required files and create required directories for outputs
load(system.file("data", "listFiles1.RData", package="jDMRgrid"))
listFiles1$file <- system.file("extdata", listFiles1$file,  package="jDMRgrid")
input_file <- system.file("extdata", "grid_population", package="jDMRgrid")
system2(
    command = "/bin/cp", 
    args = paste0("-r ", input_file, " ", "~/"))
system2(command = "/bin/mkdir", args = "~/matrix_population")

# Step 1.4: makeDMRmatrix for population-based data
makeDMRmatrix(
    contexts = c("CG", "CHG", "CHH"), postMax.out = TRUE, 
    samplelist = listFiles1, input.dir = '~/grid_population', 
    out.dir = '~/matrix_population', include.intermediate = FALSE)
}