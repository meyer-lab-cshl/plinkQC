
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/plinkQC)](https://cran.r-project.org/package=plinkQC)
[![Build
Status](https://travis-ci.org/meyer-lab-cshl/plinkQC.svg?branch=master)](https://travis-ci.org/meyer-lab-cshl/plinkQC)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/plinkQC?color=blue)](https://CRAN.R-project.org/package=plinkQC)

## <i class="fa fa-map" aria-hidden="true"></i> plinkQC

**plinkQC** is a R/CRAN package for **genotype quality control** in
genetic association studies. It makes PLINK basic statistics
(e.g.missing genotyping rates per individual, allele frequencies per
genetic marker) and relationship functions easily accessible from within
R and allows for automatic evaluation of the results.

Full documentation is available at

<http://meyer-lab-cshl.github.io/plinkQC/>.


**plinkQC** generates a per-individual and per-marker quality control
report. A step-by-step guide on how to run these analyses can be found
[here](https://meyer-lab-cshl.github.io/plinkQC/articles/plinkQC.html).

Individuals and markers that fail the quality control can subsequently
be removed with **plinkQC** to generate a new, clean dataset.

**plinkQC** facilitates an ancestry check for study individuals based on
comparison to reference datasets. The processing of the reference
datasets is documented in detail
[here](https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html).

Removal of individuals based on relationship status via **plinkQC** is
optimised to retain as many individuals as possible in the study.

## <i class="fa fa-rocket" aria-hidden="true"></i> Installation

The current github version of **plinkQC** is: 0.3.0 and can be installed
via

``` r
library(devtools)
install_github("meyer-lab-cshl/plinkQC")
```

The current CRAN version of **plinkQC** is: 0.2.3 and can be installed
via

``` r
install.packages("plinkQC")
```

A log of version changes can be found

[here](https://github.com/meyer-lab-cshl/plinkQC/blob/master/NEWS.md).

## <i class="fa fa-pencil" aria-hidden="true"></i> Citation
Meyer HV (2018) plinkQC: Genotype quality control in genetic association
studies. doi:10.5281/zenodo.3373798

