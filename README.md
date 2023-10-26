# jDMRgrid v0.2.3

## Description
jDMRgrid is a component of the jDMR toolkit, an essential tool for conducting efficient and heuristic DMR (Differentially Methylated Region) calling in large-scale epigenomic studies involving population-level analyses and control/treatment experiments. Its functionality relies on a grid-based methodology that involves dividing the genome into sliding windows or bins of customizable sizes. It determines methylation state calls for each individual sliding window based on a Hidden Markov Model (HMM) approach. Subsequently, these calls are consolidated into a matrix format by merging data from all samples. Additionally, user can perform a filtration of DMR matrix to remove non-polymorphic patterns from the merged matrix.

## Installing from the Github

To install from GitHub (development version), follow the steps given below.

##### Step 1 — Install a last version of R (>=3.6)

##### Step 2 — In R, please install all dependencies and execute the following commands:
```sh
install.packages("devtools")
library(devtools)
devtools::install_github("robert-piecyk/jDMRgrid")
```
## How to use?
Please open the [vignette](https://github.com/robert-piecyk/jDMRgrid/blob/master/vignettes/manual.pdf) file.

## Contributors
##### Robert Piecyk - robert.piecyk@tum.de
##### Rashmi Hazarika - rashmi.hazarika@tum.de
##### Yadi Shahryary - y.shahryary@tum.de
##### Frank Johannes - frank@johanneslab.org

