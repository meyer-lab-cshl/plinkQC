
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/plinkQC)](https://cran.r-project.org/package=plinkQC)
[![Build
Status](https://travis-ci.org/meyer-lab-cshl/plinkQC.svg?branch=master)](https://travis-ci.org/meyer-lab-cshl/plinkQC)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/plinkQC?color=blue)](https://cran.r-project.org/package=plinkQC)

## <i class="fa fa-map" aria-hidden="true"></i> plinkQC

**plinkQC** is a R/CRAN package for **genotype quality control** including a pre-trained ancestry classifier and relatedness filter optimized to retain the maximally unrelated sample set with highest quality.

**plinkQC** generates a per-individual and per-marker quality control report.
A step-by-step guide on how to run these analyses can be found [here](https://meyer-lab-cshl.github.io/plinkQC/articles/plinkQC.html). Individuals and markers that fail the quality control can subsequently be removed with **plinkQC** to generate a new, clean dataset.

## <i class="fa fa-rocket" aria-hidden="true"></i> Installation

The current github version of **plinkQC** is: 1.0.0 and can be installed
via

``` r
library(devtools)
install_github("meyer-lab-cshl/plinkQC")
```

The current CRAN version of **plinkQC** is: 1.0.0 and can be installed via

``` r
install.packages("plinkQC")
```

A log of version changes can be found [here](https://meyer-lab-cshl.github.io/plinkQC/news/index.html).

## <i class="fa fa-pencil" aria-hidden="true"></i> Citation

Syed M, Walter C and Meyer HV (2025) plinkQC: An Integrated Tool for Ancestry Inference, Sample Selection, and Quality Control in Population Genetics. [[bioRxiv]](https://www.biorxiv.org/content/10.1101/2025.11.25.690541v1)
