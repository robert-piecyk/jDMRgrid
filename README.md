# jDMRgrid

## Description
jDMRgrid is a tool for conducting efficient and heuristic DMR (Differentially Methylated Region) calling in population-level or control/treatment WGBS data. Its functionality relies on a grid-based methodology that involves dividing the genome into sliding windows or bins of customizable sizes. For each sample and sliding window methylation state calls are obtained via a Hidden Markov Model (HMM). These calls are then consolidated into a DMR matrix format. 

## Installing from the Github

To install from GitHub (development version), follow the steps given below.

##### Step 1 — Install a last version of R (>=3.6)

##### Step 2 — In R, please install all dependencies and execute the following commands:
```sh
install.packages("devtools")
library(devtools)
devtools::install_github("jlab-code/jDMRgrid")
```
## How to use?
Please open the [vignette](https://github.com/robert-piecyk/jDMRgrid/blob/master/vignettes/manual.pdf) file.

## Contributors
##### Robert Piecyk - robert.piecyk@tum.de
##### Rashmi Hazarika - rashmi.hazarika@tum.de
##### Yadi Shahryary - y.shahryary@tum.de
##### Frank Johannes - frank@johanneslab.org
##### Ming Zhou - ming.zhou@tum.de 

