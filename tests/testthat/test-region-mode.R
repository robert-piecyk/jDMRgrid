## Regression tests for region.mode.
##
## These exist because of a bug that shipped in the region.mode code: clustering ran on the
## cytosine GRanges, which carries each cytosine's strand, so GenomicRanges::reduce() grouped
## within strand and returned a "+" and a "-" copy of every interval.
##
## Two properties of that bug shape these tests:
##
##   1. It only surfaced from the SECOND sample onward -- the first sample's output was
##      written and the run then died inside methimpute with "each range must have an end
##      that is greater or equal to its start minus one". A single-sample smoke test passes
##      against the bug, so every case here runs two samples.
##   2. region.mode = "cluster" never errored at all. It silently returned twice the regions
##      it should have. Only a check on the region table itself catches that, so the tests
##      assert on the duplicate coordinates rather than merely on completion.

skip_if_not_installed("methimpute")

two_samples <- function() {
    e <- new.env()
    load(system.file("data", "listFiles2.RData", package = "jDMRgrid"), envir = e)
    sl <- e$listFiles2
    sl$file <- system.file("extdata", sl$file, package = "jDMRgrid")
    ## one replicate from each of the two groups
    sl <- sl[!duplicated(sl$sample), , drop = FALSE]
    stopifnot(nrow(sl) == 2L, all(file.exists(sl$file)))
    sl
}

run_mode <- function(mode, out.dir) {
    runjDMRgrid(
        out.dir = out.dir, window = 200, step = 50,
        samplelist = two_samples(), contexts = "CG",
        min.C = 5, mincov = 0, include.intermediate = FALSE,
        runName = gsub("-", "_", mode), min.C.type = "count",
        region.mode = mode,
        max.gap = if (identical(mode, "grid")) 150L else "adaptive",
        min.sites = 5)
    invisible(out.dir)
}

grid_table <- function(out.dir) {
    f <- list.files(out.dir, pattern = "\\.Rdata$", full.names = TRUE)
    expect_length(f, 1L)
    dget(f)
}

for (mode in c("grid", "cluster", "cluster-grid")) {

    test_that(sprintf("region.mode = '%s' completes for every sample", mode), {
        out <- file.path(tempfile("rm_"), gsub("-", "_", mode))
        dir.create(out, recursive = TRUE)
        expect_no_error(suppressWarnings(suppressMessages(run_mode(mode, out))))

        ## The regression: the run used to die after the first sample, leaving one file.
        calls <- list.files(out, pattern = "_CG\\.txt$")
        expect_length(calls, 2L)
        for (f in calls)
            expect_gt(nrow(data.table::fread(file.path(out, f), showProgress = FALSE)), 0L)
    })

    test_that(sprintf("region.mode = '%s' yields unstranded, unique regions", mode), {
        out <- file.path(tempfile("rm_"), gsub("-", "_", mode))
        dir.create(out, recursive = TRUE)
        suppressWarnings(suppressMessages(run_mode(mode, out)))
        g <- grid_table(out)

        expect_true(all(c("chr", "start", "end") %in% names(g)))
        expect_gt(nrow(g), 0L)
        ## degenerate ranges are what methimpute choked on
        expect_true(all(g$end >= g$start))
        ## clustering per strand produced a "+" and a "-" copy of each interval
        expect_equal(sort(unique(as.character(g$strand))), "*")
        expect_false(any(duplicated(g[, c("chr", "start", "end")])))
    })
}

test_that("cluster modes restrict regions to cytosine-bearing sequence", {
    outs <- vapply(c("grid", "cluster"), function(m) {
        o <- file.path(tempfile("rm_"), m); dir.create(o, recursive = TRUE)
        suppressWarnings(suppressMessages(run_mode(m, o)))
        nrow(grid_table(o))
    }, numeric(1))
    ## the grid tiles the whole slice; clusters cover only where cytosines occur
    expect_lt(outs[["cluster"]], outs[["grid"]])
})

test_that("an impossible min.C fails loudly rather than emptily", {
    out <- file.path(tempfile("rm_"), "empty"); dir.create(out, recursive = TRUE)
    expect_error(
        suppressWarnings(suppressMessages(runjDMRgrid(
            out.dir = out, window = 200, step = 50, samplelist = two_samples(),
            contexts = "CG", min.C = 10000, mincov = 0, include.intermediate = FALSE,
            runName = "empty", min.C.type = "count", region.mode = "grid"))),
        "no regions survived")
})
