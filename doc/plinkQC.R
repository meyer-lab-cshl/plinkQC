## ----setup, include = FALSE---------------------------------------------------
library(plinkQC)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set parameters-----------------------------------------------------------
package.dir <- find.package('plinkQC')
indir <- file.path(package.dir, 'extdata')
qcdir <- tempdir()
name <- 'data'
path2plink <- "/Users/hannah/bin/plink"

## ----individual QC,  eval=FALSE, fig.height=12, fig.width=9-------------------
# fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#                             path2plink=path2plink,
#                             interactive=TRUE, verbose=TRUE)
# 

## ----load individualQC, echo=FALSE,out.width = "500px", fig.align='center'----
par(mfrow=c(2,1), las=1)
knitr::include_graphics("individualQC.pdf")

## ----overview individual QC,fig.width=7, fig.height=7, eval=FALSE-------------
# overview_individuals <- overviewPerIndividualQC(fail_individuals,
#                                                 interactive=TRUE)

## ----load overviewIndividualQC, out.width = "500px", echo=FALSE---------------
par(mfrow=c(2,1), las=1)
knitr::include_graphics("overviewQC.png")

## ----marker QC, eval=FALSE----------------------------------------------------
# fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
#                             path2plink=path2plink,
#                             verbose=TRUE, interactive=TRUE,
#                             showPlinkOutput=FALSE)

## ----load markerQC, echo=FALSE, out.width = "500px", fig.align='center'-------
par(mfrow=c(2,1), las=1)
knitr::include_graphics("markerQC.png")

## ----overview marker QC, eval=FALSE-------------------------------------------
# overview_marker <- overviewPerMarkerQC(fail_markers, interactive=TRUE)

## ----load overviewMarkerQC, out.width = "500px", echo=FALSE-------------------
par(mfrow=c(2,1), las=1)
knitr::include_graphics("overviewMarkerQC.png")

## ----clean data, eval=FALSE---------------------------------------------------
# Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
#                             verbose=TRUE, showPlinkOutput=FALSE)

## ----eval = FALSE-------------------------------------------------------------
# name = "data.hg38"
# path2load_mat = "path/to/load_mat/merged_chrs.postQC.train.pca"
# convert_to_plink2(indir=indir, qcdir=qcdir, name=name,
#                   path2plink2 = path2plink2)
# rename_variant_identifiers(indir=qcdir, qcdir=qcdir, name=name,
#                            path2plink2 = path2plink2)
# name <- paste0(name, ".renamed")

## ----ancestry identification, eval = FALSE------------------------------------
# ancestries <- superpop_classification(indir=qcdir, qcdir=qcdir, name=name,
#                                       path2plink2 = path2plink2,
#                                       path2load_mat = path2load_mat)

## ----load ancestry, out.width = "500px", echo=FALSE---------------------------
par(mfrow=c(2,1), las=1)
knitr::include_graphics("ancestral_prediction.pdf")

## ----check sex, eval=FALSE, out.width = "500px", fig.align='center'-----------
# fail_sex <- check_sex(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
#                       verbose=TRUE, path2plink=path2plink)

## ----load checkSex, out.width = "500px", echo=FALSE, fig.align='center'-------
knitr::include_graphics("checkSex.png")

## ----check het miss, eval=FALSE, fig.height=3, fig.width=5, fig.align='center'----
# fail_het_imiss <- check_het_and_miss(indir=indir, qcdir=qcdir, name=name,
#                                      interactive=TRUE, path2plink=path2plink)

## ----load checkHetImiss, out.width = "500px", echo=FALSE, fig.align='center'----
knitr::include_graphics("checkHetImiss.png")

## ----check related, eval=FALSE, fig.height=3, fig.width=5, fig.align='center'----
# exclude_relatedness <- check_relatedness(indir=indir, qcdir=qcdir, name=name,
#                                          interactive=TRUE,
#                                          path2plink=path2plink)

## ----load checkRelatedness, out.width = "500px", echo=FALSE, fig.align='center'----
knitr::include_graphics("checkRelatedness.png")

## ----check snp missing, eval=FALSE--------------------------------------------
# fail_snpmissing <- check_snp_missingness(indir=indir, qcdir=qcdir, name=name,
#                                          interactive=TRUE,
#                                          path2plink=path2plink,
#                                          showPlinkOutput=FALSE)

## ----load snp missing, out.width = "500px", echo=FALSE, fig.align='center'----
knitr::include_graphics("snpmissingness.png")

## ----check hwe, eval=FALSE----------------------------------------------------
# fail_hwe <- check_hwe(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
#                       path2plink=path2plink, showPlinkOutput=FALSE)

## ----load hwe, out.width = "500px", echo=FALSE, fig.align='center'------------
knitr::include_graphics("hwe.png")

## ----check maf, eval=FALSE----------------------------------------------------
# fail_maf <- check_maf(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
#                       path2plink=path2plink, showPlinkOutput=FALSE)

## ----load  maf, out.width = "500px", echo=FALSE, fig.align='center'-----------
knitr::include_graphics("maf.png")

