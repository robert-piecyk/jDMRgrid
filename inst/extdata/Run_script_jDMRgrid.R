#-------------------------------------------------------------------------------
# Install jDMRgrid
#-------------------------------------------------------------------------------
rm(list=ls())
library(devtools)
devtools::install_github("robert-piecyk/jDMRgrid")
library(jDMRgrid)
home_dir <- getwd()
#-------------------------------------------------------------------------------
# Step 1: Run jDMR using grid approach
#-------------------------------------------------------------------------------
sample.files.1 <- system.file("extdata", "listFiles1.fn", package="jDMRgrid")
sample.files.2 <- system.file("extdata", "listFiles2.fn", package="jDMRgrid")
runjDMRgrid(out.dir = 'jDMRgrid_test/folder_population/grid',
            window = 200,
            step = 50,
            samplefiles = sample.files.1,
            contexts = c("CG", "CHG", "CHH"),
            min.C = 10,
            mincov = 0,
            include.intermediate = FALSE,
            runName = "Arabidopsis")

runjDMRgrid(out.dir = 'jDMRgrid_test/folder_replicate/grid',
            window = 200,
            step = 50,
            samplefiles = sample.files.2,
            contexts = c("CG", "CHG", "CHH"),
            min.C = 10,
            mincov = 0,
            include.intermediate = TRUE,
            runName = "Arabidopsis")

#-------------------------------------------------------------------------------
# Step 2: Generate DMR matrix
#-------------------------------------------------------------------------------
makeDMRmatrix(contexts = c("CG", "CHG", "CHH"),
              postMax.out = TRUE,
              samplefiles = sample.files.1,
              input.dir = 'jDMRgrid_test/folder_population/grid',
              out.dir = 'jDMRgrid_test/folder_population/matrix',
              include.intermediate = FALSE)

makeDMRmatrix(contexts = c("CG", "CHG", "CHH"),
              postMax.out = TRUE,
              samplefiles = sample.files.2,
              input.dir = 'jDMRgrid_test/folder_replicate/grid',
              out.dir = 'jDMRgrid_test/folder_replicate/matrix',
              include.intermediate = FALSE)

#-------------------------------------------------------------------------------
# Run this step only if you have multiple treatment groups and
# you want to perform pairwise comparisons with
# each of the treatment groups with control.
# This function will split the DMR matrix into
# pairwise control-treatment groups.
# Refer to section 1.2 in the manual.
# The metadata file requires an additional column "group"
#-------------------------------------------------------------------------------
split.groups(samplefiles = sample.files.2,
             input.dir = "jDMRgrid_test/folder_replicate/matrix",
             out.dir = "jDMRgrid_test/folder_replicate/matrix")

#-------------------------------------------------------------------------------
# Step 3: Filter the DMR matrix
#-------------------------------------------------------------------------------
filterDMRmatrix(epiMAF.cutoff = 0.33,
                replicate.consensus = NULL,
                data.dir = "jDMRgrid_test/folder_population/matrix",
                samplefiles = sample.files.1)

filterDMRmatrix(epiMAF.cutoff = NULL,
                replicate.consensus = 0.8,
                data.dir = "jDMRgrid_test/folder_replicate/matrix",
                samplefiles = sample.files.2)

#-------------------------------------------------------------------------------
# Step 4: Generate context-specific DMRs
#-------------------------------------------------------------------------------
output.dir.1 <- "jDMRgrid_test/folder_population/context_DMRs"
output.dir.2 <- "jDMRgrid_test/folder_replicate/context_DMRs"
context.specific.DMRs(samplefiles = sample.files.1,
                      output.dir = output.dir.1,
                      input.dir = "jDMRgrid_test/folder_population/matrix")

context.specific.DMRs(samplefiles = sample.files.2,
                      output.dir = output.dir.2,
                      input.dir = "jDMRgrid_test/folder_replicate/matrix")

#-------------------------------------------------------------------------------
# Step 5: Annotate DMRs. Please create a new folder and move all files to
# be annotated into the new folder
#-------------------------------------------------------------------------------
data.dir.pop <- "jDMRgrid_test/folder_population/annotate_DMRs"
data.dir.rep <- "jDMRgrid_test/folder_replicate/annotate_DMRs"
gff.file_promoters <- system.file("extdata/toyData",
                                  "TAIR10_promoters.gff3",
                                  package="jDMRgrid")
gff.file_TE <- system.file("extdata/toyData",
                           "TAIR10_TE.gff3",
                           package="jDMRgrid")
gff.file_genes <- system.file("extdata/toyData",
                              "TAIR10.gene.chr1.gff3",
                              package="jDMRgrid")

annotateDMRs(gff.files = c(gff.file_promoters, gff.file_TE, gff.file_genes),
             annotation = c("promoters", "TE", "gene"),
             input.dir = "jDMRgrid_test/folder_population/context_DMRs",
             gff3.out = FALSE,
             out.dir = data.dir.pop)
annotateDMRs(gff.files = c(gff.file_promoters, gff.file_TE, gff.file_genes),
             annotation = c("promoters", "TE", "gene"),
             input.dir = "jDMRgrid_test/folder_replicate/context_DMRs",
             gff3.out = FALSE,
             out.dir = data.dir.rep)
