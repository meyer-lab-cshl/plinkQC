## ----setup, include = FALSE----------------------------------------------
devtools::load_all()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set parameters------------------------------------------------------
package.dir <- find.package('plinkQC')
qcdir <- "~/RPackagesDevel/plinkQC/inst/extdata"
#qcdir <- file.path(package.dir, 'extdata')
alg <- 'data'
path2plink <- "~/bin"

## ----run qc--------------------------------------------------------------


## ----sample QC,  fig.height=12, fig.width=9------------------------------
fail_samples <- perSampleQC(qcdir=qcdir, alg=alg,
                            refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",
                                                 sep=""), 
                            refColorsFile=paste(qcdir, "/HapMap_PopColors.txt",
                                                 sep=""),
                            prefixMergedDataset="data.HapMapIII",
                            interactive=TRUE, verbose=TRUE)


