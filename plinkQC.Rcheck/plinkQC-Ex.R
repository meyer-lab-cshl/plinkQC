pkgname <- "plinkQC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('plinkQC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("check_ancestry")
### * check_ancestry

flush(stderr()); flush(stdout())

### Name: check_ancestry
### Title: Identification of individuals of divergent ancestry
### Aliases: check_ancestry

### ** Examples

## Not run: 
##D indir <- system.file("extdata", package="plinkQC")
##D name <- "data"
##D fail_ancestry <- check_ancestry(indir=indir, name=name,
##D refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",sep=""),
##D refColorsFile=paste(indir, "/HapMap_PopColors.txt", sep=""),
##D prefixMergedDataset="data.HapMapIII", interactive=FALSE,
##D run.check_ancestry=FALSE)
## End(Not run)



cleanEx()
nameEx("check_het_and_miss")
### * check_het_and_miss

flush(stderr()); flush(stdout())

### Name: check_het_and_miss
### Title: Identification of individuals with outlying missing genotype or
###   heterozygosity rates
### Aliases: check_het_and_miss

### ** Examples

 ## Not run: 
##D indir <- system.file("extdata", package="plinkQC")
##D name <- "data"
##D fail_het_miss <- check_het_and_miss(indir=indir, name=name,
##D run.check_het_and_miss=FALSE, interactive=FALSE)
## End(Not run)



cleanEx()
nameEx("check_hwe")
### * check_hwe

flush(stderr()); flush(stdout())

### Name: check_hwe
### Title: Identification of SNPs showing a significant deviation from
###   Hardy-Weinberg- equilibrium (HWE)
### Aliases: check_hwe

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D fail_hwe <- check_hwe(indir=indir, qcdir=qcdir, name=name, interactive=FALSE,
##D verbose=TRUE, path2plink=path2plink)
## End(Not run)



cleanEx()
nameEx("check_maf")
### * check_maf

flush(stderr()); flush(stdout())

### Name: check_maf
### Title: Identification of SNPs with low minor allele frequency
### Aliases: check_maf

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D fail_maf <- check_maf(indir=indir, qcdir=qcdir, name=name, macTh=15,
##D interactive=FALSE, verbose=TRUE, path2plink=path2plink)
## End(Not run)



cleanEx()
nameEx("check_relatedness")
### * check_relatedness

flush(stderr()); flush(stdout())

### Name: check_relatedness
### Title: Identification of related individuals
### Aliases: check_relatedness

### ** Examples

## Not run: 
##D indir <- system.file("extdata", package="plinkQC")
##D name <- 'data'
##D relatednessQC <- check_relatedness(indir=indir, name=name, interactive=FALSE,
##D run.check_relatedness=FALSE)
## End(Not run)



cleanEx()
nameEx("check_sex")
### * check_sex

flush(stderr()); flush(stdout())

### Name: check_sex
### Title: Identification of individuals with discordant sex information
### Aliases: check_sex

### ** Examples

 ## Not run: 
##D indir <- system.file("extdata", package="plinkQC")
##D name <- "data"
##D fail_sex <- check_sex(indir=indir, name=name, run.check_sex=FALSE,
##D interactive=FALSE, verbose=FALSE)
## End(Not run)



cleanEx()
nameEx("check_snp_missingness")
### * check_snp_missingness

flush(stderr()); flush(stdout())

### Name: check_snp_missingness
### Title: Identification of SNPs with high missingness rate
### Aliases: check_snp_missingness

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D fail_snp_missingness <- check_snp_missingness(qcdir=qcdir, indir=indir,
##D name=name, interactive=FALSE, verbose=TRUE, path2plink=path2plink)
## End(Not run)



cleanEx()
nameEx("cleanData")
### * cleanData

flush(stderr()); flush(stdout())

### Name: cleanData
### Title: Create plink dataset with individuals and markers passing
###   quality control
### Aliases: cleanData

### ** Examples

package.dir <- find.package('plinkQC')
indir <- file.path(package.dir, 'extdata')
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D # Run individual QC checks
##D fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
##D refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
##D refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
##D prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE)
##D 
##D # Run marker QC checks
##D fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name)
##D 
##D # Create new dataset of indiviudals and markers passing QC
##D ids_all <- cleanData(indir=indir, qcdir=qcdir, name=name, macTh=15,
##D verbose=TRUE, path2plink=path2plink, filterAncestry=FALSE,
##D filterRelated=TRUE)
## End(Not run)



cleanEx()
nameEx("evaluate_check_ancestry")
### * evaluate_check_ancestry

flush(stderr()); flush(stdout())

### Name: evaluate_check_ancestry
### Title: Evaluate results from PLINK PCA on combined study and reference
###   data
### Aliases: evaluate_check_ancestry

### ** Examples

## Not run: 
##D qcdir <- system.file("extdata", package="plinkQC")
##D name <- "data"
##D fail_ancestry <- evaluate_check_ancestry(indir=qcdir, name=name,
##D refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
##D refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
##D prefixMergedDataset="data.HapMapIII", interactive=FALSE)
## End(Not run)



cleanEx()
nameEx("evaluate_check_het_and_miss")
### * evaluate_check_het_and_miss

flush(stderr()); flush(stdout())

### Name: evaluate_check_het_and_miss
### Title: Evaluate results from PLINK missing genotype and heterozygosity
###   rate check.
### Aliases: evaluate_check_het_and_miss

### ** Examples

qcdir <- system.file("extdata", package="plinkQC")
name <- "data"
## Not run: 
##D fail_het_miss <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,
##D interactive=FALSE)
## End(Not run)



cleanEx()
nameEx("evaluate_check_relatedness")
### * evaluate_check_relatedness

flush(stderr()); flush(stdout())

### Name: evaluate_check_relatedness
### Title: Evaluate results from PLINK IBD estimation.
### Aliases: evaluate_check_relatedness

### ** Examples

qcdir <- system.file("extdata", package="plinkQC")
name <- 'data'
## Not run: 
##D relatednessQC <- evaluate_check_relatedness(qcdir=qcdir, name=name,
##D interactive=FALSE)
## End(Not run)



cleanEx()
nameEx("evaluate_check_sex")
### * evaluate_check_sex

flush(stderr()); flush(stdout())

### Name: evaluate_check_sex
### Title: Evaluate results from PLINK sex check.
### Aliases: evaluate_check_sex

### ** Examples

qcdir <- system.file("extdata", package="plinkQC")
name <- "data"
## Not run: 
##D fail_sex <- evaluate_check_sex(qcdir=qcdir, name=name, interactive=FALSE,
##D verbose=FALSE)
## End(Not run)



cleanEx()
nameEx("overviewPerIndividualQC")
### * overviewPerIndividualQC

flush(stderr()); flush(stdout())

### Name: overviewPerIndividualQC
### Title: Overview of per sample QC
### Aliases: overviewPerIndividualQC

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
## Not run: 
##D fail_individuals <- perIndividualQC(qcdir=qcdir, indir=indir, name=name,
##D refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
##D refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
##D prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
##D do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
##D do.run_check_sex=FALSE, do.run_check_ancestry=FALSE)
##D 
##D overview <- overviewPerIndividualQC(fail_individuals)
## End(Not run)



cleanEx()
nameEx("overviewPerMarkerQC")
### * overviewPerMarkerQC

flush(stderr()); flush(stdout())

### Name: overviewPerMarkerQC
### Title: Overview of per marker QC
### Aliases: overviewPerMarkerQC

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
# All quality control checks
## Not run: 
##D fail_markers <- perMarkerQC(qcdir=qcdir, indir=indir, name=name,
##D interactive=FALSE, verbose=TRUE, path2plink=path2plink)
##D overview <- overviewPerMarkerQC(fail_markers)
## End(Not run)



cleanEx()
nameEx("perIndividualQC")
### * perIndividualQC

flush(stderr()); flush(stdout())

### Name: perIndividualQC
### Title: Quality control for all individuals in plink-dataset
### Aliases: perIndividualQC

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
# All quality control checks
## Not run: 
##D fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
##D refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
##D refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
##D prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
##D do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
##D do.run_check_sex=FALSE, do.run_check_ancestry=FALSE)
##D 
##D # Only check sex and missingness/heterozygosity
##D fail_sex_het_miss <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
##D dont.check_ancestry=TRUE, dont.check_relatedness=TRUE,
##D interactive=FALSE, verbose=FALSE)
## End(Not run)



cleanEx()
nameEx("perMarkerQC")
### * perMarkerQC

flush(stderr()); flush(stdout())

### Name: perMarkerQC
### Title: Quality control for all markers in plink-dataset
### Aliases: perMarkerQC

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
name <- "data"
path2plink <- '/path/to/plink'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
# All quality control checks
## Not run: 
##D fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
##D interactive=FALSE, verbose=TRUE, path2plink=path2plink)
## End(Not run)



cleanEx()
nameEx("run_check_ancestry")
### * run_check_ancestry

flush(stderr()); flush(stdout())

### Name: run_check_ancestry
### Title: Run PLINK principal component analysis
### Aliases: run_check_ancestry

### ** Examples

indir <- system.file("extdata", package="plinkQC")
qcdir <- tempdir()
prefixMergedDataset <- 'data.HapMapIII'
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D run <- run_check_ancestry(indir=indir, qcdir=qcdir, prefixMergedDataset)
## End(Not run)



cleanEx()
nameEx("run_check_heterozygosity")
### * run_check_heterozygosity

flush(stderr()); flush(stdout())

### Name: run_check_heterozygosity
### Title: Run PLINK heterozygosity rate calculation
### Aliases: run_check_heterozygosity

### ** Examples

indir <- system.file("extdata", package="plinkQC")
name <- 'data'
qcdir <- tempdir()
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D run <- run_check_heterozygosity(indir=indir, qcdir=qcdir, name=name)
## End(Not run)



cleanEx()
nameEx("run_check_missingness")
### * run_check_missingness

flush(stderr()); flush(stdout())

### Name: run_check_missingness
### Title: Run PLINK missingness rate calculation
### Aliases: run_check_missingness

### ** Examples

indir <- system.file("extdata", package="plinkQC")
name <- 'data'
qcdir <- tempdir()
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D run <- run_check_missingnessness(indir=indir, qcdir=qcdir, name=name)
## End(Not run)



cleanEx()
nameEx("run_check_relatedness")
### * run_check_relatedness

flush(stderr()); flush(stdout())

### Name: run_check_relatedness
### Title: Run PLINK IBD estimation
### Aliases: run_check_relatedness

### ** Examples

indir <- system.file("extdata", package="plinkQC")
name <- 'data'
qcdir <- tempdir()
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name)
## End(Not run)



cleanEx()
nameEx("run_check_sex")
### * run_check_sex

flush(stderr()); flush(stdout())

### Name: run_check_sex
### Title: Run PLINK sexcheck
### Aliases: run_check_sex

### ** Examples

indir <- system.file("extdata", package="plinkQC")
name <- 'data'
qcdir <- tempdir()
# the following code is not run on package build, as the path2plink on the
# user system is not known.
## Not run: 
##D run <- run_check_sex(indir=indir, qcdir=qcdir, name=name)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
