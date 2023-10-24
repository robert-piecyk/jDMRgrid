\name{splitGroups}
\alias{splitGroups}
\title{Split DMR matrix into pairwise groups control/treatment groups}
\usage{
splitGroups(
    samplelist, postMax.out = FALSE, contexts = c("CG", "CHG", "CHH"), 
    input.dir, out.dir)
}
\arguments{
\item{samplelist}{
DataFrame object containing information about
file, sample, and replicate. For control/treatment data, an 
additional column specifying the replicates is required. (DataFrame object)
}

\item{postMax.out}{
Logical specifying whether or not the 
postMax probabilities should be divided into control/treatment groups. 
By default, this option is set to FALSE. (logical)
}

\item{contexts}{
Vector of cytosine contexts selected for DMR calling. By 
default, this option is set for all 3 cytosine contexts CG, CHG, and CHH. 
(character vector)
}

\item{input.dir}{
Path to the input directory with DMR matrices. (character)
}

\item{out.dir}{
Path to the output directory with DMR matrices split by
control/treatment groups. (character)
}
}
\value{
DMR matrices split by sample/treatment groups, saved in the 
output directory with DMR matrices.
}
\description{
This function splits DMR matrices for state-calls, rc.meth.lvl, and optionally
postMax probabilities by control/treatment groups.
}
\examples{
#------------------------------------------------------------------------------
# REPLICATE-BASED DATASET; STEP III
#------------------------------------------------------------------------------
# Step 1.5: Copy required files and create required directories for outputs
load(system.file("data", "listFiles2.RData", package = "jDMRgrid"))
listFiles2$file <- system.file("extdata", listFiles2$file, package = "jDMRgrid")
system2(command = "/bin/mkdir", args = "~/matrix_CT_replicate")
input_file <- system.file("extdata", "matrix_replicate", package = "jDMRgrid")
system2(
    command = "/bin/cp", 
    args = paste0("-r ", input_file, " ", "~/"))

# Step 1.6: Split control/treatment for replicate-based data
splitGroups(
    samplelist = listFiles2,
    postMax.out = TRUE, contexts = c("CG", "CHG", "CHH"), 
    input.dir = '~/matrix_replicate', out.dir = '~/matrix_CT_replicate')
}