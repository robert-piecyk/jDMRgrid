\name{context.specific.DMRs}
\alias{context.specific.DMRs}
\title{subset context specific DMRs}
\usage{
context.specific.DMRs(samplelist, input.dir, output.dir, if.filtered = FALSE)
}
\arguments{
\item{samplelist}{DataFrame object containing information about
file, sample and replicate. For control/treatment data an 
additional column specifying the replicates is required. (DataFrame object)}

\item{input.dir}{Path to the input directory containing DMR matrix files. 
Looks for files with suffix("_StateCalls.txt" and "_rcMethlvl.txt"). 
(character)}

\item{output.dir}{Path to the output directory, which will contain set of DMRs
containing context-specific DMRs, non-CG DMRs and multi-context DMRs. 
(character)}

\item{if.filtered}{Logical to specify if we should use filtered or
non-filtered matrices. By default it is set to FALSE. (logical)}
}

\value{
output txt files containing context-specific DMRs for five categories:
-> CG-only, 
-> CHG-only,
-> CHH-only, 
-> non-CG (occurs for CHG and/or CHH, but not for CG context),
-> multi-context (occurs for CG and CHG and CHH contexts)
}
\description{
Output context-specific DMRs (only for a given context, nonCG or all context)
}
\examples{
load(system.file("data", "listFiles2.RData", package="jDMRgrid"))
listFiles2$file <- system.file("extdata", listFiles2$file,  package="jDMRgrid")
system2(command = "/bin/mkdir", args = "~/context.specific.DMRs")
context.specific.DMRs(
    samplelist = listFiles2, 
    input.dir = system.file(
        "extdata", "matrix_CT_replicate/", package="jDMRgrid"), 
    output.dir = '~/context.specific.DMRs', 
    if.filtered = FALSE)
}