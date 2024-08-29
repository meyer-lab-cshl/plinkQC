## ----setup knitr, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----check ancestry, eval=FALSE, fig.height=3, fig.width=5, fig.align='center'----
#  library(plinkQC)
#  indir <- system.file("extdata", package="plinkQC")
#  qcdir <- "~/qcdir"
#  name <- 'data'
#  refname <- 'all_hg38'
#  prefixMergedDataset <- paste(name, ".", refname, sep="")
#  
#  exclude_ancestry <-
#      evaluate_check_ancestry(indir=qcdir, name=name,
#                              prefixMergedDataset=prefixMergedDataset,
#                              refSamplesFile=paste(indir,
#                                                   "/Genomes1000_ID2Pop.txt",
#                                                   sep=""),
#                              refColorsFile=paste(indir,
#                                                  "/Genomes1000_PopColors.txt",
#                                                  sep=""),
#                              interactive=TRUE)

## ----load ancestry, out.width = "500px", echo=FALSE, fig.align='center'-------
knitr::include_graphics("checkAncestry.png")

