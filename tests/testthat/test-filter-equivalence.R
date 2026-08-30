## Equivalence tests for the vectorised filters introduced in 0.3.0.
##
## filterReplicateConsensus(), merge.bins() and filterEpiMAF() were rewritten from per-row and
## per-region loops into vectorised form. The rewrites were verified once against the original
## implementations, but that check lived outside the package and so did not run again. The
## originals are kept in helper-reference-implementations.R and the fast paths are checked
## against them here on every build.
##
## Case counts are smaller than the original one-off sweeps (which used 162 and 108
## combinations) because the reference implementations are slow by construction. The
## parameters still span the cases that previously exposed defects: three state levels, more
## than two sample groups, and bins that overlap, abut or are disjoint.

skip_if_not_installed("GenomicRanges")

norm <- function(x) {
    if (is.null(x)) return(NULL)
    lapply(x, function(d) {
        d <- as.data.frame(d)
        d$seqnames <- as.character(d$seqnames)
        d <- d[order(d$seqnames, d$start, d$end), , drop = FALSE]
        rownames(d) <- NULL
        for (k in names(d)) if (is.numeric(d[[k]])) d[[k]] <- round(d[[k]], 9)
        d
    })
}

same <- function(a, b) {
    if (is.null(a) && is.null(b)) return(TRUE)
    if (is.null(a) != is.null(b)) return(FALSE)
    isTRUE(all.equal(norm(a), norm(b), check.attributes = FALSE))
}

test_that("filterReplicateConsensus matches the per-region implementation", {
    set.seed(7)
    mismatches <- 0L; n <- 0L
    for (nstate in c(2L, 3L))
      for (ngrp in c(2L, 3L))
        for (nrep in c(2L, 3L))
          for (cons in c(0.5, 0.8, 1.0))
            for (dct in c(0, 0.5)) {
                cols <- unlist(lapply(seq_len(ngrp),
                    function(g) paste0("S", g, "_rep", seq_len(nrep))))
                nreg <- 25L
                M <- matrix(sample(0:(nstate - 1L), nreg * length(cols), TRUE),
                            nreg, length(cols), dimnames = list(NULL, cols))
                co <- data.frame(seqnames = "1", start = seq_len(nreg) * 100L,
                                 end = seq_len(nreg) * 100L + 99L)
                sc <- data.frame(co, M, check.names = FALSE)
                rc <- data.frame(co, matrix(round(runif(nreg * length(cols)), 4),
                                 nreg, dimnames = list(NULL, cols)), check.names = FALSE)
                n <- n + 1L
                a <- ref_filterReplicateConsensus(sc, rc, cons, dct)
                b <- jDMRgrid::filterReplicateConsensus(sc, rc, cons, dct)
                if (!same(a, b)) mismatches <- mismatches + 1L
            }
    expect_gt(n, 40L)          # the sweep actually ran
    expect_equal(mismatches, 0L)
})

test_that("merge.bins matches the split-and-reduce implementation", {
    set.seed(11)
    merge_bins <- get("merge.bins", envir = asNamespace("jDMRgrid"))
    mismatches <- 0L; n <- 0L
    for (nsamp in c(2L, 3L))
      for (nstate in c(2L, 3L))
        for (nchr in c(1L, 2L))
          for (step in c(50L, 200L, 300L)) {          # overlapping, abutting, disjoint
              nb <- 40L
              starts <- unlist(lapply(seq_len(nchr),
                  function(k) seq(1L, by = step, length.out = nb)))
              chrs <- rep(as.character(seq_len(nchr)), each = nb)
              cols <- paste0("s", seq_len(nsamp), "_rep1")
              M <- matrix(sample(0:(nstate - 1L), length(starts) * nsamp, TRUE),
                          ncol = nsamp, dimnames = list(NULL, cols))
              co <- data.frame(seqnames = chrs, start = starts, end = starts + 199L)
              sc <- data.frame(co, M, check.names = FALSE)
              rc <- data.frame(co, matrix(round(runif(length(starts) * nsamp), 4),
                               ncol = nsamp, dimnames = list(NULL, cols)),
                               check.names = FALSE)
              n <- n + 1L
              if (!same(ref_merge.bins(rc, sc), merge_bins(rc, sc)))
                  mismatches <- mismatches + 1L
          }
    expect_gt(n, 20L)
    expect_equal(mismatches, 0L)
})

test_that("filterEpiMAF reproduces the reference computation over the same samples", {
    ## The reference reads columns 4:(NCOL - 1). The epiMAF column is appended AFTER that
    ## line, so in real use there is no trailing column and the reference silently excluded
    ## the last SAMPLE. Handing it one spare column makes both implementations see the same
    ## sample set, and the arithmetic can then be compared directly.
    set.seed(13)
    mismatches <- 0L; n <- 0L
    for (nstate in c(2L, 3L))
      for (nsamp in c(4L, 6L))
        for (maf in c(0.1, 0.25, 0.34)) {
            nreg <- 30L
            cols <- paste0("s", seq_len(nsamp))
            M <- matrix(sample(0:(nstate - 1L), nreg * nsamp, TRUE), nreg, nsamp,
                        dimnames = list(NULL, cols))
            co <- data.frame(seqnames = "1", start = seq_len(nreg) * 100L,
                             end = seq_len(nreg) * 100L + 99L)
            sc <- data.frame(co, M, check.names = FALSE)
            rc <- data.frame(co, matrix(round(runif(nreg * nsamp), 4), nreg,
                             dimnames = list(NULL, cols)), check.names = FALSE)
            n <- n + 1L
            a <- ref_filterEpiMAF(cbind(sc, spare = NA_real_), rc, maf)
            b <- jDMRgrid::filterEpiMAF(sc, rc, maf)
            ok <- if (is.null(a) && is.null(b)) TRUE
                  else if (is.null(a) != is.null(b)) FALSE
                  else isTRUE(all.equal(a[[1]]$epiMAF, b[[1]]$epiMAF)) &&
                       identical(a[[1]]$start, b[[1]]$start) &&
                       isTRUE(all.equal(norm(list(a[[2]]))[[1]],
                                        norm(list(b[[2]]))[[1]], check.attributes = FALSE))
            if (!ok) mismatches <- mismatches + 1L
        }
    expect_gt(n, 10L)
    expect_equal(mismatches, 0L)
})

test_that("filterEpiMAF uses every sample, unlike the reference", {
    ## Regression for the off-by-one column range. With four samples the reference could only
    ## reach 1/3 = 0.333, so an epiMAF.cutoff of 0.33 matched nothing however the data looked.
    set.seed(17)
    nreg <- 20L; nsamp <- 4L
    cols <- paste0("s", seq_len(nsamp))
    M <- matrix(sample(0:1L, nreg * nsamp, TRUE), nreg, nsamp, dimnames = list(NULL, cols))
    M[1, ] <- c(1L, 0L, 0L, 0L)                     # exactly one minor call in four samples
    co <- data.frame(seqnames = "1", start = seq_len(nreg) * 100L,
                     end = seq_len(nreg) * 100L + 99L)
    sc <- data.frame(co, M, check.names = FALSE)
    rc <- data.frame(co, matrix(round(runif(nreg * nsamp), 4), nreg,
                     dimnames = list(NULL, cols)), check.names = FALSE)

    b <- jDMRgrid::filterEpiMAF(sc, rc, 0.30)
    expect_false(is.null(b))
    expect_true(any(b[[1]]$epiMAF == 0.25))         # 1 of 4, unreachable for the reference

    a <- ref_filterEpiMAF(sc, rc, 0.30)             # same call, unfixed range
    expect_true(is.null(a) || all(a[[1]]$epiMAF > 0.30))
})
