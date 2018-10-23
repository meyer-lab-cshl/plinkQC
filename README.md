
<!-- README.md is generated from README.Rmd. Please edit that file -->
<i class="fa fa-map" aria-hidden="true"></i> plinkQC
----------------------------------------------------

**plinkQC** is a R/CRAN package for **genotype quality control** in genetic association studies. It makes PLINK basic statistics (e.g.missing genotyping rates per individual, allele frequencies per genetic marker) and relationship functions easily accessible from within R and allows for automatic evaluation of the results.

**plinkQC** facilitates an ancestry check for study individuals based on comparison to reference datasets. The processing of the reference datasets is documented in detail [here](articles/AncestryCheck.html)

**plinkQC** generates a per-individual and per-marker quality control report.

Individuals and markers that fail the quality control can subsequently be removed with **plinkQC** to generate a new, clean dataset.

Removal of individuals based on relationship status via **plinkQC** is optimised to retain as many individuals as possible in the study.

<i class="fa fa-cog" aria-hidden="true"></i> Getting Started
------------------------------------------------------------

The current github version of **plinkQC** is: 0.1.0 and can be installed via

``` r
library(devtools)
install_github("HannahVMeyer/PhenotypeSimulator")
```

**plinkQC** will soon be available on CRAN!

Full documentation is available at <http://HannahVMeyer.github.io/plinkQC/>

A log of version changes can be found [here](https://github.com/HannahVMeyer/plinkQC/blob/master/NEWS.md).
