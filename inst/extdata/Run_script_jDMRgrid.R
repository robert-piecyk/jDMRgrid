#-------------------------------------------------------------------------------
# Install jDMRgrid
#-------------------------------------------------------------------------------
rm(list=ls())
library(devtools)
devtools::install_github("robert-piecyk/jDMRgrid")
library(jDMR)

#-------------------------------------------------------------------------------
# Step 1: Run jDMR using grid approach
#-------------------------------------------------------------------------------
runjDMRgrid(out.dir = 'folder_population/grid',
            window = 200,
            step = 50,
            samplefiles = system.file("extdata", "listFiles1.fn", package="jDMRgrid"),
            contexts = c("CG", "CHG", "CHH"),
            min.C = 10,
            mincov = 0,
            include.intermediate = FALSE,
            runName = "Arabidopsis")

runjDMRgrid(out.dir = 'folder_replicate/grid',
            window = 200,
            step = 50,
            samplefiles = system.file("extdata", "listFiles2.fn", package="jDMRgrid"),
            contexts = c("CG", "CHG", "CHH"),
            min.C = 10,
            mincov = 0,
            include.intermediate = FALSE,
            runName = "Arabidopsis")

#-------------------------------------------------------------------------------
# Step 2: Generate DMR matrix
#-------------------------------------------------------------------------------
makeDMRmatrix(contexts = c("CG", "CHG", "CHH"),
              postMax.out = TRUE,
              samplefiles = system.file("extdata", "listFiles1.fn", package="jDMRgrid"),
              input.dir = 'folder_population/grid',
              out.dir = 'folder_population/matrix',
              include.intermediate = FALSE)

makeDMRmatrix(contexts = c("CG", "CHG", "CHH"),
              postMax.out = TRUE,
              samplefiles = system.file("extdata", "listFiles2.fn", package="jDMRgrid"),
              input.dir = 'folder_replicate/grid',
              out.dir = 'folder_replicate/matrix',
              include.intermediate = FALSE)

#-------------------------------------------------------------------------------
# Run this step only if you have multiple treatment groups and you want to
# perform pairwise comparisons with each of the treatment groups with control.
# This function will split the DMR matrix into pairwise control-treatment groups.
# Refer to section 1.2 in the manual. The metadata file requires an additional column
# "group"
#-------------------------------------------------------------------------------
split.groups(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMRgrid"),
             input.dir="folder_replicate/matrix",
             out.dir="folder_replicate/matrix")

#-------------------------------------------------------------------------------
# Step 3: Filter the DMR matrix
#-------------------------------------------------------------------------------
filterDMRmatrix(epiMAF.cutoff = 0.33,
                replicate.consensus = NULL,
                data.dir = "folder_population/matrix",
                samplefiles = system.file("extdata", "listFiles1.fn", package="jDMRgrid"))

filterDMRmatrix(epiMAF.cutoff = NULL,
                replicate.consensus = 0.8,
                data.dir = "folder_replicate/matrix",
                samplefiles = system.file("extdata", "listFiles2.fn", package="jDMRgrid"))

#-------------------------------------------------------------------------------
# Step 4: Generate context-specific DMRs
#-------------------------------------------------------------------------------
context.specific.DMRs(samplefiles=system.file("extdata", "listFiles1.fn", package="jDMR"),
                      data.dir="folder_population/matrix")

context.specific.DMRs(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
                      data.dir="folder_replicate/matrix")

#-------------------------------------------------------------------------------
# Step 5: Annotate DMRs. Please create a new folder and move all files to be annotated
# into the new folder
#-------------------------------------------------------------------------------
data.dir <- "/folder_population/annotate_DMRs"
gff.AT <- "/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3"
gff.TE <- "/Annotations/TAIR10_TE.gff3"
gff.pr <- "/Annotations/TAIR10_promoters.gff3"

annotateDMRs(gff.files=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","TE","promoters"), #string containing annotation types
             input.dir=data.dir,
             gff3.out=FALSE,
             out.dir=data.dir)
