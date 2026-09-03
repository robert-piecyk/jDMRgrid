# jDMRgrid 0.3.1

## Bug fixes

* `rankDMRs(method = "permutation")` computed the pooled-null p-value by scanning the
  whole null once per region, which is quadratic in the region count. On a real epiRIL
  matrix of 1,012,452 CG regions the call ran 22 hours at 200 permutations without
  finishing, and 12 hours at 50. Sorting the null once and counting with `findInterval()`
  is O(n log n); the same matrix now takes 7.6 minutes at 200 permutations and 4.4 at 50.
  The counts are identical, ties included, so p-values and FDR are unchanged.

* The permutation branch had no test coverage; the suite exercised only `betabinom`.
  Adds the counting identity against a direct scan, an end-to-end check that a larger
  statistic never earns a larger p against the same pooled null, and the RNG-restoration
  guarantee that only this branch calls.

## Diagnostics

* When no region reaches FDR < 0.05, `rankDMRs()` now reports why. The existing warning
  fires only when fewer than 20 label assignments exist, which is the condition for the
  p-value floor to be attainable; that is rarely the binding constraint. With many regions
  the pooled null's tail is filled by other regions' permuted statistics, so at 1e6 regions
  even a one-in-a-million tail holds a few hundred values and no region clears BH however
  many permutations are added. The message gives the strongest statistic, how many null
  values exceed it, the tail fraction, and the BH threshold the region count implies.

  Measured on that matrix: tail fraction 1.28e-06 at 50 permutations and 1.08e-06 at 200,
  with 0 of 1,012,452 regions at FDR < 0.05 in both. Rank genome-wide with
  `method = "betabinom"`.

## Packaging

* `DESCRIPTION` declares the `LICENSE` file, which `R CMD check` flagged as an undeclared
  top-level file.

# jDMRgrid 0.3.0

Changes in version 0.3.0
* New function rankDMRs(): orders surviving regions by effect size, a
  beta-binomial Wald statistic, or a permutation null
* New argument filterDMRmatrix(min.posterior=): filters on the HMM posterior,
  which earlier versions wrote and never read
* New argument runjDMRgrid(region.mode=): "grid", "cluster" or "cluster-grid",
  with max.gap="adaptive" deriving the clustering gap per context
* filterDMRmatrix(diff.ct=) is exposed, keeping its former fixed value of 0.5
* Default changed: filterDMRmatrix(if.mergingBins=) is now TRUE, so one
  differential region is one row instead of one row per overlapping window
* Default changed: makeDMRmatrix(postMax.out=) and splitGroups(postMax.out=) are
  now TRUE, as min.posterior and rankDMRs both require that matrix
* runjDMRgrid(min.C.type=) has no default and must be given; the meaning of
  min.C changed between versions and the two readings differ by orders of
  magnitude
* filterEpiMAF used columns 4:(NCOL-1) while the epiMAF column is appended
  afterwards, so the last sample was excluded; with four samples it used three
  and an epiMAF.cutoff of 0.33 could never match
* Fixed: rankDMRs ranked perfectly separated regions last with NA statistics;
  the Wald variance is now taken at the pooled proportion
* Fixed: rankDMRs combined the filtered region table with the unfiltered
  posterior matrix; posteriors are matched by overlap
* Fixed: rankDMRs no longer leaves the caller's random stream altered
* Fixed: region.mode clustering ran per strand, returning a "+" and a "-" copy
  of every interval; cluster-grid failed after the first sample
* Fixed: filterReplicateConsensus errored on an empty matrix
* Fixed: Bismark CX report input works end to end
* deriveSampleName() accepts gzipped input
* filterReplicateConsensus, merge.bins and filterEpiMAF rewritten in vectorised
  form, verified against the previous implementations
* Toy data replaced with real coverage: 91,971 cytosines over 250 kb of
  A. thaliana chromosome 1 at 11.5x, run through METHimpute
* inst/extdata reduced from 11 MB to 3.7 MB
* Test suite added, run with R CMD check
* Vignette rewritten, in HTML and PDF

# jDMRgrid 0.2.4

Changes in version 0.2.4 (2023-10-30)
* jDMRgrid accepts Bismark CX reports (together with FASTA file; used in Bismark analysis)
* fixed errors for R v4.0 and up
