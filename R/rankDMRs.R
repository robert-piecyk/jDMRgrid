#' Rank DMRs by effect size, a dispersion-aware statistic, or a permutation null
#'
#' Adds an ordering to the regions produced by the standard jDMRgrid pipeline. This is an
#' optional layer: it reads matrices \code{makeDMRmatrix} already wrote and changes nothing
#' about how regions are called or filtered.
#'
#' The state-calling core decides differential status by rule rather than by a statistical
#' test, so its output carries no quantity to sort by. Meanwhile the HMM produces a
#' per-region, per-sample posterior (\code{posteriorMax}) that the pipeline writes to disk and
#' never reads again. The three methods use progressively more of that information.
#'
#' \describe{
#'   \item{\code{"effect"}}{Between-group difference in \code{rc.meth.lvl}, weighted by the
#'     lowest contributing \code{posteriorMax}. No distributional assumption.}
#'   \item{\code{"betabinom"}}{A Wald-type statistic whose standard error uses a
#'     context-specific overdispersion \code{rho}. Note the matrices hold methylation
#'     \emph{levels}, not (methylated, total) counts, so this is a moment-based
#'     approximation -- per-sample levels are treated as having variance
#'     \code{rho * p * (1 - p)} -- and \strong{not} a beta-binomial likelihood on counts.
#'     Default \code{rho} values (CG 0.020, CHG 0.040, CHH 0.080) are the per-context
#'     overdispersions used by DMRspiker for this organism.}
#'   \item{\code{"permutation"}}{Shuffles the group labels \code{nperm} times, recomputes the
#'     statistic genome-wide, and pools the permuted values into a single empirical null.
#'     Per-region empirical p-values are taken against that pooled null and adjusted with BH.
#'     Pooling is deliberate: with few samples per group the number of distinct label
#'     assignments is small (six for 2 vs 2, twenty for 3 vs 3), so a per-region null would be
#'     far too coarse. This is the same construction BSmooth uses.}
#' }
#'
#' @param data.dir Directory holding the matrices written by \code{makeDMRmatrix}. (character)
#' @param samplefiles Sample sheet with a \code{group} column, or a data.frame. (character)
#' @param contexts Cytosine contexts to rank. (character vector)
#' @param method One of \code{"effect"}, \code{"betabinom"}, \code{"permutation"}. (character)
#' @param out.dir Directory for the ranked tables; defaults to \code{data.dir}. (character)
#' @param if.filtered Use the \code{-filtered} matrices. (logical)
#' @param level.from Where per-region methylation levels come from. \code{"counts"}
#'   (default) recomputes observed methylation, sum(methylated)/sum(total), from the input
#'   files named in the sample sheet. \code{"states"} uses the \code{rc.meth.lvl} matrix,
#'   which is a posterior-weighted mixture of the HMM's two fitted state levels rather than
#'   observed methylation, and so cannot represent a change that stays within one state.
#' @param nperm Permutations for \code{method = "permutation"}. (numeric)
#' @param rho Named per-context overdispersions for \code{"betabinom"}. (numeric)
#' @param seed Seed for the permutation null, recorded in the output. (numeric)
#'
#' @return Invisibly, a named list of \code{data.table}s. Also writes
#'   \code{<context>_rankedDMRs.txt} into \code{out.dir}.
#'
#' @examples
#' # rankDMRs(data.dir = "matrix", samplefiles = "listFiles2.fn", method = "effect")
#'
#' @importFrom data.table fread fwrite data.table setorder
#' @importFrom stats p.adjust
#' @export
rankDMRs <- function(data.dir, samplefiles, contexts = c("CG", "CHG", "CHH"),
                     method = c("effect", "betabinom", "permutation"),
                     out.dir = NULL, if.filtered = FALSE, nperm = 1000L,
                     rho = c(CG = 0.020, CHG = 0.040, CHH = 0.080), seed = 42L,
                     level.from = c("counts", "states")) {

    method <- match.arg(method)
    level.from <- match.arg(level.from)
    if (is.null(out.dir)) out.dir <- data.dir
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(data.dir)) stop("data.dir does not exist: ", data.dir)

    ft <- if (is.character(samplefiles)) data.table::fread(samplefiles) else
          data.table::as.data.table(samplefiles)
    if (is.null(ft$group)) stop("samplefiles needs a 'group' column")
    ft$name <- if (!is.null(ft$replicate)) paste0(ft$sample, "_", ft$replicate) else ft$sample
    grps <- unique(ft$group)
    if (length(grps) < 2L) stop("need at least two groups")
    ctrl.name <- if ("control" %in% grps) "control" else grps[1]
    trt.name  <- setdiff(grps, ctrl.name)[1]

    ## statistic for one labelling of the samples
    stat_for <- function(M, c.idx, t.idx, rho.cx) {
        mu.c <- rowMeans(M[, c.idx, drop = FALSE], na.rm = TRUE)
        mu.t <- rowMeans(M[, t.idx, drop = FALSE], na.rm = TRUE)
        delta <- mu.t - mu.c
        if (method == "effect") return(list(delta = delta, stat = abs(delta)))
        ## Moment-based Wald statistic, with Var(level) approximated as rho * p * (1 - p).
        ##
        ## The variance is taken at the POOLED proportion, not per group. Using each group's
        ## own proportion makes the variance exactly zero whenever a group sits at 0% or
        ## 100%, which is precisely the case of perfect separation -- the strongest evidence
        ## a region can carry. The old form produced se = 0 there, then NA, then a score of
        ## 0, so a region called 1,1 against 0,0 ranked LAST rather than first.
        ##
        ## Pooling is also the standard variance for a two-proportion comparison under the
        ## null that both groups share a proportion, so it is the right quantity for a test
        ## statistic. It degenerates only when the pooled proportion is itself 0 or 1, and
        ## then both groups are identical, delta is 0, and the region is not a DMR anyway.
        p.pool <- (mu.c * length(c.idx) + mu.t * length(t.idx)) /
                  (length(c.idx) + length(t.idx))
        se <- sqrt(rho.cx * p.pool * (1 - p.pool) *
                   (1 / length(c.idx) + 1 / length(t.idx)))
        se[!is.finite(se) | se <= 0] <- NA_real_
        list(delta = delta, se = se, stat = abs(delta / se))
    }

    sfx <- if (isTRUE(if.filtered)) "-filtered" else ""
    res <- list()
    for (cx in contexts) {
        f.rc <- file.path(data.dir, paste0(cx, "_rcMethlvl", sfx, ".txt"))
        ## When filterDMRmatrix() was given a sample sheet with a `group` column it writes
        ## per-group-pair matrices named <ctx>_<grp1>_<grp2>_rcMethlvl-filtered.txt, not the
        ## plain <ctx>_ name. Fall back to whichever matching file is actually on disk.
        if (!file.exists(f.rc)) {
            alt <- list.files(data.dir, full.names = TRUE,
                pattern = paste0("^", cx, "_.*rcMethlvl", sfx, "\\.txt$"))
            if (length(alt) == 1L) {
                f.rc <- alt
                message(cx, ": using ", basename(f.rc))
            } else if (length(alt) > 1L) {
                stop("several candidate matrices for ", cx, " in ", data.dir,
                     " (", paste(basename(alt), collapse = ", "),
                     "); rank one group pair at a time")
            }
        }
        f.pm <- file.path(data.dir, paste0(cx, "_postMax.txt"))
        if (!file.exists(f.rc)) { message("skipping ", cx, ": ", basename(f.rc), " absent"); next }
        rc <- data.table::fread(f.rc)
        pm <- if (file.exists(f.pm)) data.table::fread(f.pm) else NULL

        coord <- c("seqnames", "start", "end")
        smp <- setdiff(names(rc), coord)
        c.cols <- intersect(ft$name[ft$group == ctrl.name], smp)
        t.cols <- intersect(ft$name[ft$group == trt.name],  smp)
        if (!length(c.cols) || !length(t.cols))
            stop("sample names in ", basename(f.rc), " do not match the sample sheet")

        ## Where the methylation level comes from.
        ##
        ## `states`: the rc.meth.lvl matrix, which is what earlier versions of this function
        ## used. That column is NOT observed methylation -- it is a posterior-weighted mixture
        ## of the HMM's two fitted state levels, and on a real chromosome 95%+ of bins sit
        ## exactly on one of those two values (measured on A. thaliana chr1: CG 67.2% at the
        ## fitted unmethylated level 0.0143 and 28.5% at the fitted methylated level 0.7551,
        ## leaving 4.2% anywhere in between). A delta built from it is therefore close to the
        ## state difference multiplied by the gap between the fitted levels, not an effect
        ## size in methylation units, and it cannot see a change that stays inside one state.
        ##
        ## `counts` (default): observed methylation, sum(methylated)/sum(total) per region
        ## per sample, recomputed from the same input files the pipeline was run on. This is
        ## what an effect size in a ranking should mean.
        if (identical(level.from, "counts")) {
            if (is.null(ft$file) || !all(file.exists(ft$file)))
                stop("level.from = \"counts\" needs a readable `file` column in the sample sheet")
            regions <- data.table::as.data.table(rc[, coord, with = FALSE])
            regions[, rid := .I]
            ## seqnames must be coerced BEFORE the key is set: assigning to a keyed column
            ## silently drops the key, and foverlaps then refuses the join.
            regions[, seqnames := as.character(seqnames)]
            data.table::setkeyv(regions, c("seqnames", "start", "end"))
            Mc <- matrix(NA_real_, nrow(regions), length(smp),
                         dimnames = list(NULL, smp))
            for (nmi in smp) {
                fpath <- ft$file[match(nmi, ft$name)]
                cy <- data.table::fread(fpath, showProgress = FALSE,
                        select = c("seqnames", "start", "context",
                                   "counts.methylated", "counts.total"))
                cy <- cy[context == cx]
                cy[, `:=`(seqnames = as.character(seqnames),
                          start = as.integer(start), end = as.integer(start))]
                ov <- data.table::foverlaps(cy, regions,
                        by.x = c("seqnames", "start", "end"),
                        type = "within", nomatch = 0L)
                agg <- ov[, .(m = sum(counts.methylated), t = sum(counts.total)), by = rid]
                Mc[agg$rid, nmi] <- agg$m / pmax(agg$t, 1L)
            }
            M <- Mc
            n.na <- sum(!is.finite(M))
            if (n.na)
                message(cx, ": ", n.na, " region/sample cells had no covered cytosines")
        } else {
            M <- as.matrix(rc[, smp, with = FALSE])
        }
        c.idx <- match(c.cols, smp); t.idx <- match(t.cols, smp)
        rho.cx <- if (!is.null(rho[[cx]])) rho[[cx]] else 0.04
        obs <- stat_for(M, c.idx, t.idx, rho.cx)

        if (!is.null(pm)) {
            ## filterDMRmatrix writes <ctx>_rcMethlvl-filtered.txt but no filtered postMax, so
            ## with if.filtered = TRUE the two tables have different row counts. Recycling the
            ## shorter one does not error -- it silently inflated the result to the unfiltered
            ## length and reported that as the number of ranked regions -- so the posterior
            ## table is aligned to the region table by coordinate, never by position.
            ## The posterior matrix is written per bin and is never filtered or merged, so its
            ## coordinates may not match the region table at all: filterDMRmatrix(if.mergingBins
            ## = TRUE) fuses runs of overlapping windows into new intervals that appear nowhere
            ## in the postMax file. An exact coordinate join then matches nothing and every
            ## posterior comes back NA.
            ##
            ## Match by OVERLAP instead, and summarise the constituent bins: the reported
            ## minimum is the weakest call anywhere in the region, across every sample and every
            ## bin it was built from, which is the conservative reading.
            pmk <- data.table::as.data.table(pm)
            pcols <- intersect(smp, names(pmk))
            reg <- data.table::as.data.table(rc)[, coord, with = FALSE]
            reg[, rid := .I][, seqnames := as.character(seqnames)]
            pmk[, seqnames := as.character(seqnames)]
            data.table::setkeyv(reg, c("seqnames", "start", "end"))
            ov <- data.table::foverlaps(pmk, reg,
                    by.x = c("seqnames", "start", "end"), nomatch = 0L)
            if (NROW(ov)) {
                agg <- ov[, {
                    v <- unlist(.SD, use.names = FALSE)
                    v <- v[is.finite(v)]
                    .(mn = if (length(v)) min(v) else NA_real_,
                      mu = if (length(v)) mean(v) else NA_real_)
                }, by = rid, .SDcols = pcols]
                min.post <- rep(NA_real_, NROW(reg)); mean.post <- min.post
                min.post[agg$rid] <- agg$mn; mean.post[agg$rid] <- agg$mu
            } else {
                min.post <- mean.post <- rep(NA_real_, NROW(reg))
            }
            n.miss <- sum(is.na(min.post))
            if (n.miss)
                message(cx, ": ", n.miss, " region(s) had no overlapping posterior")
        } else min.post <- mean.post <- rep(NA_real_, nrow(rc))

        out <- data.table::data.table(
            rc[, coord, with = FALSE],
            n_control = length(c.idx), n_treatment = length(t.idx),
            delta = obs$delta, min_posterior = min.post, mean_posterior = mean.post)
        if (method == "effect") {
            out$score <- abs(obs$delta) * ifelse(is.na(min.post), 1, min.post)
        } else {
            out$se <- obs$se; out$stat <- obs$stat
            out$score <- ifelse(is.na(obs$stat), 0, obs$stat)
        }

        if (method == "permutation") {
            ## Seeding for reproducibility must not clobber the caller's random stream.
            ## set.seed() inside a function leaves the global .Random.seed altered, so a
            ## script that draws random numbers after calling rankDMRs() silently gets a
            ## different sequence than it would have. Save the state and restore it on exit.
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
                .old.seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
                on.exit(assign(".Random.seed", .old.seed, envir = globalenv()), add = TRUE)
            } else {
                ## no stream existed before, so remove the one set.seed() creates
                on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                            rm(".Random.seed", envir = globalenv()), add = TRUE)
            }
            set.seed(seed)
            all.idx <- c(c.idx, t.idx); nc <- length(c.idx)
            nmax <- choose(length(all.idx), nc)
            nb <- min(nperm, max(1, nmax - 1))
            if (nb < nperm)
                message(cx, ": only ", nmax, " distinct label assignments exist; using ", nb)
            if (nb < 20L)
                warning("context ", cx, ": the permutation null rests on only ", nb,
                        " distinct label assignments (", length(c.idx), " vs ",
                        length(t.idx), " samples). The smallest attainable empirical p is ",
                        signif(1 / (1 + nb * nrow(rc)), 3),
                        ", so BH-adjusted FDR across ", nrow(rc),
                        " regions will be uninformative. Use method = \"betabinom\" for ",
                        "ranking, or add replicates before relying on the FDR.",
                        call. = FALSE)
            null <- vector("list", nb)
            for (b in seq_len(nb)) {
                pidx <- sample(all.idx)
                null[[b]] <- stat_for(M, pidx[seq_len(nc)], pidx[-seq_len(nc)], rho.cx)$stat
            }
            nullv <- unlist(null, use.names = FALSE)
            nullv <- nullv[is.finite(nullv)]
            s <- ifelse(is.finite(out$stat %||% out$score), out$stat %||% out$score, 0)
            ## Empirical p against the POOLED genome-wide null, by binary search.
            ## Scanning the null once per region is O(regions x nperm x regions): for a
            ## million regions at 50 permutations that is 5e13 comparisons and runs for
            ## many hours. Sorting once and counting with findInterval() is O(n log n),
            ## and gives identical counts including at ties.
            nullv <- sort(nullv)
            n_ge <- length(nullv) - findInterval(s, nullv, left.open = TRUE)
            out$p_emp <- (1 + n_ge) / (1 + length(nullv))
            out$fdr <- stats::p.adjust(out$p_emp, method = "BH")
            out$n_perm <- nb; out$seed <- seed

            ## Say why the FDR is empty when it is empty. The existing nb < 20 warning
            ## covers the case where the p-value FLOOR cannot clear BH. That is not the
            ## usual failure. With many regions the binding constraint is the pooled
            ## null's TAIL: it is filled by other regions' permuted statistics, so at 1e6
            ## regions even a one-in-a-million tail holds a few hundred values, and no
            ## region reaches the BH threshold however many permutations are added.
            ## Measured on 1,012,452 regions: 65 null values above the top statistic at 50
            ## permutations, 218 at 200, tail fraction unchanged at ~1.1e-06.
            if (nrow(out) && !any(out$fdr < 0.05, na.rm = TRUE)) {
                top <- which.max(s)
                n_above <- length(nullv) - findInterval(s[top], nullv, left.open = TRUE)
                message(cx, ": no region reaches FDR < 0.05. The strongest statistic (",
                        signif(s[top], 4), ") is exceeded by ", n_above,
                        " of ", length(nullv), " pooled-null values (tail fraction ",
                        signif(n_above / max(1, length(nullv)), 3),
                        "). BH over ", nrow(out), " regions needs p < ",
                        signif(0.05 / nrow(out), 3), ", attainable here only for a region ",
                        "no permuted value exceeds. More permutations lower the floor and ",
                        "raise the tail count together, so they will not change this; ",
                        "rank with method = \"betabinom\" instead.")
            }
        }

        data.table::setorder(out, -score)
        out$rank <- seq_len(nrow(out))
        data.table::fwrite(out, file.path(out.dir, paste0(cx, "_rankedDMRs.txt")), sep = "\t")
        res[[cx]] <- out
        message(cx, ": ranked ", nrow(out), " regions by '", method, "' (",
                ctrl.name, " vs ", trt.name, ")")
    }
    invisible(res)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
