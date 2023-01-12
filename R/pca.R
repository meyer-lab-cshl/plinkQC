library(plinkQC)
indir <- system.file("extdata", package="plinkQC")
qcdir='~/qcdir'
name <- 'data.no_ac_gt_snps.hg38'
refname <- 'all_hg38'
prefixMergedDataset <- paste(name, ".", refname, sep="")

exclude_ancestry <-
  evaluate_check_ancestry(indir=qcdir, name=name,
                          prefixMergedDataset=prefixMergedDataset,
                          refSamplesFile=paste(indir, "/Genomes1000_ID2Pop.txt",
                                               sep=""), 
                          refColorsFile=paste(indir, "/Genomes1000_PopColors.txt",
                                              sep=""),
                          interactive=TRUE)
