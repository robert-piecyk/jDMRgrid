# jDMRgrid

## Description

jDMRgrid calls differentially methylated regions (DMRs) in population-level or
control/treatment WGBS data. The genome is divided into sliding windows, a Hidden Markov Model
assigns a methylation state to each window for every sample and sequence context, and those
calls are assembled into a DMR matrix. A region is reported where the state differs between
groups.

Input can be METHimpute output or Bismark cytosine (CX) reports.

## Installation

Requires R >= 4.0.

```r
install.packages("devtools")
devtools::install_github("robert-piecyk/jDMRgrid", ref = "v0.3.0")
```

## Quick start

```r
library(jDMRgrid)

runjDMRgrid(out.dir = "grid", window = 200, step = 50, samplelist = samples,
            contexts = c("CG", "CHG", "CHH"), min.C = 10,
            min.C.type = "percentile", include.intermediate = FALSE)

makeDMRmatrix(contexts = c("CG", "CHG", "CHH"), samplelist = samples,
              input.dir = "grid", out.dir = "matrix")

splitGroups(samplelist = samples, input.dir = "matrix", out.dir = "matrix")

filterDMRmatrix(replicate.consensus = 0.8, data.dir = "matrix",
                samplelist = samples, min.posterior = 0.9)

rankDMRs(data.dir = "matrix", samplefiles = samples, contexts = "CG",
         method = "betabinom", out.dir = "ranked", if.filtered = TRUE)
```

`min.C.type` has no default and must be given. It decides whether `min.C` is an absolute count
of cytosines per window or a percentile of the observed distribution, and the two differ by
orders of magnitude.

## Documentation

The vignette covers both designs end to end, the output files, region definition modes, scale,
and troubleshooting.

```r
browseVignettes("jDMRgrid")
```

The source is `vignettes/manual.Rmd`. `R CMD build` renders the HTML; the PDF is produced with

```r
rmarkdown::render("vignettes/manual.Rmd",
                  output_format = BiocStyle::pdf_document(toc = TRUE, number_sections = TRUE))
```

Do not pass `output_file` to that call. BiocStyle inserts the author block by post-processing,
and a custom output name skips it.

## What is new in 0.3.0

* `rankDMRs()` orders surviving regions by effect size, a beta-binomial Wald statistic, or a
  permutation null.
* `filterDMRmatrix(min.posterior=)` uses the HMM posterior, which earlier versions wrote and
  never read.
* `region.mode` on `runjDMRgrid()` offers `grid`, `cluster` and `cluster-grid`, with
  `max.gap="adaptive"`.
* Defaults changed: `if.mergingBins` and `postMax.out` are now `TRUE`, and `min.C.type` is
  required. See `NEWS.md`, these alter the output of an unmodified script.
* Bismark CX reports work end to end, and gzipped input is accepted.
* Test suite added.

## License

GPL-3. See the `License` field in `DESCRIPTION`.

## Authors

Robert Stefan Piecyk (author and maintainer).

Rashmi Rekha Hazarika, Yadollah Shahryary Dizaji, Frank Johannes and Ming Zhou
contributed to earlier versions.
