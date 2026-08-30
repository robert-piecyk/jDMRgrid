## Regression tests for rankDMRs' test statistic.
##
## The betabinom statistic took Var(level) = rho * p * (1 - p) at each group's OWN
## proportion. That is exactly zero when a group sits at 0% or 100%, which is the case of
## perfect separation -- the strongest evidence a region can carry. The zero became an NA se,
## then a score of 0, so a region called 1,1 against 0,0 ranked LAST. The variance is now
## taken at the pooled proportion.

write_matrix <- function(dir, cx, coords, values, samples) {
    d <- data.frame(coords, values, check.names = FALSE)
    names(d) <- c("seqnames", "start", "end", samples)
    data.table::fwrite(d, file.path(dir, cx), sep = "\t")
}

fixture <- function() {
    d <- tempfile("rank_"); dir.create(d)
    coords  <- data.frame(seqnames = c(1L, 1L, 1L),
                          start = c(100L, 300L, 500L),
                          end   = c(299L, 499L, 699L))
    samples <- c("WT_rep1", "WT_rep2", "mut_rep1", "mut_rep2")
    ## row 1: perfect separation (1,1 vs 0,0)
    ## row 2: a strong but non-degenerate difference
    ## row 3: a weak difference
    lv <- rbind(c(1.0, 1.0, 0.0, 0.0),
                c(0.9, 0.8, 0.2, 0.1),
                c(0.55, 0.50, 0.45, 0.40))
    st <- rbind(c(1L, 1L, 0L, 0L), c(1L, 1L, 0L, 0L), c(1L, 1L, 0L, 0L))
    pm <- matrix(0.99, 3, 4)
    write_matrix(d, "CG_rcMethlvl-filtered.txt", coords, lv, samples)
    write_matrix(d, "CG_StateCalls-filtered.txt", coords, st, samples)
    write_matrix(d, "CG_postMax.txt", coords, pm, samples)
    list(dir = d, sheet = data.frame(
        file = NA_character_, sample = c("WT", "WT", "mut", "mut"),
        replicate = c("rep1", "rep2", "rep1", "rep2"),
        group = c("control", "control", "treatment", "treatment"),
        stringsAsFactors = FALSE))
}

test_that("perfect separation gets a finite statistic and ranks first", {
    f <- fixture()
    out <- tempfile("ranked_"); dir.create(out)
    suppressMessages(suppressWarnings(
        rankDMRs(data.dir = f$dir, samplefiles = f$sheet, contexts = "CG",
                 method = "betabinom", out.dir = out, if.filtered = TRUE,
                 level.from = "states")))
    r <- data.table::fread(file.path(out, "CG_rankedDMRs.txt"))

    perfect <- r[r$start == 100L, ]
    expect_equal(nrow(perfect), 1L)
    ## the defect: se was NA here, so stat was NA and score 0
    expect_true(is.finite(perfect$se))
    expect_true(is.finite(perfect$stat))
    expect_gt(perfect$se, 0)
    ## and it must come out on top, not last
    expect_equal(perfect$rank, 1L)
})

test_that("the statistic decreases as the group difference shrinks", {
    f <- fixture()
    out <- tempfile("ranked_"); dir.create(out)
    suppressMessages(suppressWarnings(
        rankDMRs(data.dir = f$dir, samplefiles = f$sheet, contexts = "CG",
                 method = "betabinom", out.dir = out, if.filtered = TRUE,
                 level.from = "states")))
    r <- data.table::fread(file.path(out, "CG_rankedDMRs.txt"))
    data.table::setkey(r, start)
    expect_true(all(is.finite(r$stat)))
    ## rows were built as perfect > strong > weak
    expect_gt(r[start == 100L]$stat, r[start == 300L]$stat)
    expect_gt(r[start == 300L]$stat, r[start == 500L]$stat)
})

test_that("no ranked row is emitted for regions absent from the matrix", {
    f <- fixture()
    out <- tempfile("ranked_"); dir.create(out)
    suppressMessages(suppressWarnings(
        rankDMRs(data.dir = f$dir, samplefiles = f$sheet, contexts = "CG",
                 method = "betabinom", out.dir = out, if.filtered = TRUE,
                 level.from = "states")))
    r <- data.table::fread(file.path(out, "CG_rankedDMRs.txt"))
    expect_equal(nrow(r), 3L)
    expect_false(any(is.na(r$seqnames)))
})
