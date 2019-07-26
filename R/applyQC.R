#' Create plink dataset with individuals and markers passing quality control
#'
#' Individuals that fail per-individual QC and markers that fail
#' per-marker QC are removed from indir/name.bim/.bed/.fam and a new, dataset
#' with the remaining individuals and markers is created as
#' qcdir/name.clean.bim/.bed/.fam.
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where results will be written to.
#' If \code{\link{perIndividualQC}} was conducted, this directory should be the
#' same as qcdir specified in \code{\link{perIndividualQC}}, i.e. it contains
#' name.fail.IDs with IIDs of individuals that failed QC. User needs writing
#' permission to qcdir. Per default, qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param filterRelated [logical] Set to exclude samples that failed relatedness
#' check (via \code{\link{check_relatedness}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/name.fail-IBD.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.evaluate_check_relatedness set to TRUE).
#' @param filterAncestry [logical] Set to exclude samples that failed ancestry
#' check (via \code{\link{check_ancestry}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/name.fail-ancestry.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.check_ancestry set to TRUE).
#' @param filterHeterozygosity [logical] Set to exclude samples that failed
#' check for outlying heterozygosity rates (via
#' \code{\link{check_het_and_miss}} or
#' \code{\link{perIndividualQC}}). Requires file qcdir/name.fail-het.IDs
#' (automatically created by \code{\link{perIndividualQC}} if
#' do.evaluate_check_het_and_miss set to TRUE).
#' @param filterSampleMissingness [logical] Set to exclude samples that failed
#' check for excessive missing genotype rates (via
#' \code{\link{check_het_and_miss}} or
#' \code{\link{perIndividualQC}}). Requires file qcdir/name.fail-imiss.IDs
#' (automatically created by \code{\link{perIndividualQC}} if
#' do.evaluate_check_het_and_miss set to TRUE).
#' @param filterSex [logical] Set to exclude samples that failed the sex
#' check (via \code{\link{check_sex}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/name.fail-sexcheck.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.evaluate_check_sex set to TRUE).
#' @param filterHWE [logical] Set to exclude markers that fail HWE exact test
#' (via \code{\link{check_hwe}} or \code{\link{perMarkerQC}}). Requires hweTh to
#' be set.
#' @param filterMAF [logical] Set to exclude markers that fail minor allele
#' frequency or minor allele count threshold (via \code{\link{check_maf}} or
#' \code{\link{perMarkerQC}}). Requires mafTh or macTh to be set.
#' @param filterSNPMissingness [logical] Set to exclude markers that have
#' excessive missing rates across samples (via
#' \code{\link{check_snp_missingness}} or \code{\link{perMarkerQC}}). Requires
#' lmissTh to be set.
#' @inheritParams checkPlink
#' @inheritParams check_maf
#' @inheritParams check_hwe
#' @inheritParams check_snp_missingness
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @return names [list] with i) passIDs, containing a [data.frame] with family
#' [FID] and individual [IID] IDs of samples that pass the QC, ii) failIDs,
#' containing a [data.frame] with family [FID] and individual [IID] IDs of
#' samples that fail the QC.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir, 'extdata')
#' qcdir <- tempdir()
#' name <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # Run individual QC checks
#' fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE)
#'
#' # Run marker QC checks
#' fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name)
#'
#' # Create new dataset of indiviudals and markers passing QC
#' ids_all <- cleanData(indir=indir, qcdir=qcdir, name=name, macTh=15,
#' verbose=TRUE, path2plink=path2plink, filterAncestry=FALSE,
#' filterRelated=TRUE)
#' }

cleanData <- function(indir, name, qcdir=indir,
                      filterSex=TRUE, filterHeterozygosity=TRUE,
                      filterSampleMissingness=TRUE,
                      filterAncestry=TRUE, filterRelated=TRUE,
                      filterSNPMissingness=TRUE, lmissTh=0.01,
                      filterHWE=TRUE, hweTh=1e-5,
                      filterMAF=TRUE, macTh=20, mafTh=NULL,
                      path2plink=NULL, verbose=FALSE,
                      showPlinkOutput=TRUE) {
    sampleFilter <- c(filterAncestry, filterRelated, filterSex,
                      filterHeterozygosity, filterSampleMissingness)
    markerFilter <- c(filterHWE, filterMAF, filterSNPMissingness)

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    if (!any(c(sampleFilter, markerFilter))) {
        stop("No per-sample and per-marker filters chosen")
    }
    if (!any(sampleFilter)) {
        message("No per-sample filter chosen, carry on with removing markers ",
                "that fail per-marker QCs")
    }
    if (!any(markerFilter)) {
        message("No per-marker filter chosen, carry on with removing samples ",
                "that fail per-samples QCs")
    }
    checkFormat(prefix)
    path2plink <- checkPlink(path2plink)

    # Remove remove.IDs file if already existing
    removeIDs <- NULL
    IDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                   data.table=FALSE, stringsAsFactors=FALSE,
                                   header=FALSE)
    if (any(sampleFilter)) {
        if (filterRelated) {
            if (!file.exists(paste(out, ".fail-IBD.IDs", sep=""))){
                stop("filterRelated is TRUE but file ", out,
                     ".fail-IBD.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed relatedness check")
                }
                removeIDs <- rbind(removeIDs,
                                   data.table::fread(paste(out, ".fail-IBD.IDs",
                                                     sep=""),
                                                     data.table=FALSE,
                                                     stringsAsFactors=FALSE,
                                                     header=FALSE))
            }
        }
        if (filterAncestry) {
            if (!file.exists(paste(out, ".fail-ancestry.IDs", sep=""))){
                stop("filterAncestry is TRUE but file ", out,
                     ".fail-ancestry.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed ancestry check")
                }
                removeIDs <- rbind(removeIDs,
                                   data.table::fread(paste(out,
                                                           ".fail-ancestry.IDs",
                                                           sep=""),
                                                     data.table=FALSE,
                                                     stringsAsFactors=FALSE,
                                                     header=FALSE))
            }
        }
        if (filterHeterozygosity) {
            if (!file.exists(paste(out, ".fail-het.IDs", sep=""))){
                stop("filterHeterozygosity is TRUE but file ", out,
                     ".fail-het.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed heterozygosity ",
                            "check")
                }
                removeIDs <- rbind(removeIDs,
                                   data.table::fread(paste(out, ".fail-het.IDs",
                                                           sep=""),
                                                     data.table=FALSE,
                                                     stringsAsFactors=FALSE,
                                                     header=FALSE))
            }
        }
        if (filterSampleMissingness) {
            if (!file.exists(paste(out, ".fail-imiss.IDs", sep=""))){
                stop("filterSampleMissingness is TRUE but file ", out,
                     ".fail-imiss.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed missingness check")
                }
                removeIDs <- rbind(removeIDs,
                                   data.table::fread(paste(out,
                                                           ".fail-imiss.IDs",
                                                           sep=""),
                                                     data.table=FALSE,
                                                     stringsAsFactors=FALSE,
                                                     header=FALSE))
            }
        }
        if (filterSex) {
            if (!file.exists(paste(out, ".fail-sexcheck.IDs", sep=""))){
                stop("filterSex is TRUE but file ", out,
                     ".fail-sexcheck.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed sex check")
                }
                removeIDs <- rbind(removeIDs,
                                   data.table::fread(paste(out,
                                                           ".fail-sexcheck.IDs",
                                                           sep=""),
                                                     data.table=FALSE,
                                                     stringsAsFactors=FALSE,
                                                     header=FALSE))
            }
        }
        # ensure unique IDs in remove.IDs
        removeIDs <- removeIDs[!duplicated(removeIDs),]
        if (nrow(removeIDs) == nrow(IDs)) {
                stop("All samples are flagged as .fail.IDs ",
                     "no samples remaining to generate the QCed dataset.")
        }
        if (verbose) message("Write file with remove IDs")
        write.table(removeIDs, paste(out, ".remove.IDs", sep=""),
                    col.names=FALSE, row.names=FALSE, quote=FALSE)
        keepIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                     data.table=FALSE, stringsAsFactors=FALSE,
                                     header=FALSE)
        keepIDs <- keepIDs[!keepIDs[,2] %in% removeIDs[,2],1:2]
        remove <- c("--remove", paste(out, ".remove.IDs", sep=""))
        fail_samples <- nrow(removeIDs)
    } else {
        keepIDs <- IDs
        removeIDs <- NULL
        remove <- NULL
        fail_samples <- 0
    }

    hwe <- NULL
    maf <- NULL
    missing <- NULL

    if (filterHWE) {
        if (is.null(hweTh)) {
            stop("filterHWE is TRUE but hweTh not specified")
        } else {
           hwe <- c("--hwe", hweTh)
        }
    }
    if (filterMAF) {
        if (is.null(mafTh) && is.null(macTh)) {
            stop("filterMAF is TRUE but neither mafTh or macTh are provided")
        }

        all_samples <-  R.utils::countLines(paste(prefix, ".fam", sep=""))
        keep_samples <- as.numeric(all_samples) - fail_samples

        if (!is.null(mafTh) && !is.null(macTh)) {
            if (verbose) {
                message("Both mafTh and macTh provided, macTh=", macTh,
                        " is used (corresponds to mafTh=", round(mafTh, 6), ")")
            }
        } else if (!is.null(mafTh)) {
            if(is.null(macTh)) macTh <- mafTh*(2*keep_samples)
            if (verbose) {
                message("The mafTh is ", mafTh, " which corresponds to a mcfTh=",
                        macTh)
            }
        } else {
            if(is.null(mafTh)) mafTh <- macTh/(2*keep_samples)
            if (verbose) {
                message("The macTh is ", macTh," which corresponds to a mafTh=",
                        round(mafTh, 6))
            }
        }
        maf <- c("--maf", mafTh)
    }
    if (filterSNPMissingness) {
        if (is.null(lmissTh)) {
            stop("filterSNPMissingness is TRUE but lmissTh not specified")
        } else {
            missing <- c("--geno", lmissTh)
        }
    }
    if (verbose) message("Remove individual IDs and markers IDs that failed QC")
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, remove, maf, hwe, missing,
                          "--make-bed", "--out", paste(out, ".clean", sep="")),
           std_out=showPlinkOutput, std_err=showPlinkOutput)
    return(list(passIDs=keepIDs, failIDs=removeIDs))
}
