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
#' @inheritParams checkFiltering
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
#' # Run qc on all samples and markers in the dataset
#' ## Run individual QC checks
#' fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' path2plink=path2plink)
#'
#' ## Run marker QC checks
#' fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
#' path2plink=path2plink)
#'
#' ## Create new dataset of individuals and markers passing QC
#' ids_all <- cleanData(indir=indir, qcdir=qcdir, name=name, macTh=15,
#' verbose=TRUE, path2plink=path2plink, filterAncestry=FALSE,
#' filterRelated=TRUE)
#'
#' # Run qc on subset of samples and markers in the dataset
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#'
#' fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' dont.check_ancestry = TRUE, interactive=FALSE, verbose=FALSE,
#' highlight_samples = highlight_samples[,2], highlight_type = "label",
#' remove_individuals = remove_individuals_file, path2plink=path2plink)
#'
#' ## Run marker QC checks
#' fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
#' path2plink=path2plink)
#'
#' ## Create new dataset of individuals and markers passing QC
#' ids_all <- cleanData(indir=indir, qcdir=qcdir, name=name, macTh=15,
#' verbose=TRUE, path2plink=path2plink, filterAncestry=FALSE,
#' remove_individuals = remove_individuals_file)
#' }

cleanData <- function(indir, name, qcdir=indir,
                      filterSex=TRUE, filterHeterozygosity=TRUE,
                      filterSampleMissingness=TRUE,
                      filterAncestry=TRUE, filterRelated=TRUE,
                      filterSNPMissingness=TRUE, lmissTh=0.01,
                      filterHWE=TRUE, hweTh=1e-5,
                      filterMAF=TRUE, macTh=20, mafTh=NULL,
                      path2plink=NULL, verbose=FALSE,
                      keep_individuals=NULL,
                      remove_individuals=NULL,
                      exclude_markers=NULL,
                      extract_markers=NULL,
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
    args_filter <- checkFiltering(extract_markers=extract_markers,
                                  exclude_markers=exclude_markers)
    removeIDs <- checkRemoveIDs(prefix=prefix,
                                remove_individuals=remove_individuals,
                                keep_individuals=keep_individuals)
    allIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                data.table=FALSE, stringsAsFactors=FALSE,
                                header=FALSE)
    allIDs <- allIDs[,1:2]

    if (any(sampleFilter)) {
        if (filterRelated) {
            fail_ibd_ids <- paste0(out, ".fail-IBD.IDs")
            if (!file.exists(fail_ibd_ids)){
                stop("filterRelated is TRUE but file ", out,
                     ".fail-IBD.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed relatedness check")
                }
                if (file.size(fail_ibd_ids) == 0) {
                    if (verbose) {
                        message("No individuals failed relatedness check")
                    }
                } else {
                    removeIDs <- rbind(removeIDs,
                                       data.table::fread(fail_ibd_ids,
                                                         data.table=FALSE,
                                                         stringsAsFactors=FALSE,
                                                         header=FALSE))
                }
            }
        }
        if (filterAncestry) {
            fail_ancestry_ids <- paste0(out, ".fail-ancestry.IDs")
            if (!file.exists(fail_ancestry_ids)){
                stop("filterAncestry is TRUE but file ", out,
                     ".fail-ancestry.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed ancestry check")
                }
                if (file.size(fail_ancestry_ids) == 0) {
                    if (verbose) {
                        message("No individuals failed ancestry check")
                    }
                } else {
                    removeIDs <- rbind(removeIDs,
                                       data.table::fread(fail_ancestry_ids,
                                                         data.table=FALSE,
                                                         stringsAsFactors=FALSE,
                                                         header=FALSE))
                }
            }
        }
        if (filterHeterozygosity) {
            fail_het_ids <- paste0(out, ".fail-het.IDs")
            if (!file.exists(fail_het_ids)){
                stop("filterHeterozygosity is TRUE but file ", out,
                     ".fail-het.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed heterozygosity ",
                            "check")
                }
                if (file.size(fail_het_ids) == 0) {
                    if (verbose) {
                        message("No individuals failed heterozygosity check")
                    }
                } else {
                    removeIDs <- rbind(removeIDs,
                                       data.table::fread(fail_het_ids,
                                                         data.table=FALSE,
                                                         stringsAsFactors=FALSE,
                                                         header=FALSE))
                }
            }
        }
        if (filterSampleMissingness) {
            fail_imiss_ids <- paste0(out, ".fail-imiss.IDs")
            if (!file.exists(fail_imiss_ids)){
                stop("filterSampleMissingness is TRUE but file ", out,
                     ".fail-imiss.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed missingness check")
                }
                if (file.size(fail_imiss_ids) == 0) {
                    if (verbose) {
                        message("No individuals failed missingness check")
                    }
                } else {
                    removeIDs <- rbind(removeIDs,
                                       data.table::fread(fail_imiss_ids,
                                                         data.table=FALSE,
                                                         stringsAsFactors=FALSE,
                                                         header=FALSE))
                }
            }
        }
        if (filterSex) {
            fail_sexcheck_ids <- paste0(out, ".fail-sexcheck.IDs")
            if (!file.exists(fail_sexcheck_ids)){
                stop("filterSex is TRUE but file ", out,
                     ".fail-sexcheck.IDs does not exist")
            } else {
                if (verbose) {
                    message("Read individual IDs that failed sex check")
                }
                if (file.size(fail_sexcheck_ids) == 0) {
                    if (verbose) {
                        message("No individuals failed sex check")
                    }
                } else {
                    removeIDs <- rbind(removeIDs,
                                       data.table::fread(fail_sexcheck_ids,
                                                         data.table=FALSE,
                                                         stringsAsFactors=FALSE,
                                                         header=FALSE))
                }
            }
        }
        # ensure unique IDs in remove.IDs
        if (!is.null(removeIDs)) {
            removeIDs <- removeIDs[!duplicated(removeIDs),]
            if (nrow(removeIDs) == nrow(allIDs)) {
                    stop("All samples are flagged as .fail.IDs ",
                         "no samples remaining to generate the QCed dataset.")
            }
            if (verbose) message("Write file with remove IDs")
            write.table(removeIDs, paste(out, ".remove.IDs", sep=""),
                        col.names=FALSE, row.names=FALSE, quote=FALSE)
            remove <- c("--remove", paste(out, ".remove.IDs", sep=""))
            fail_samples <- nrow(removeIDs)
        } else {
            remove <- NULL
            fail_samples <- 0
        }
    } else {
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
           hwe <- c("--hwe midp", hweTh)
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
                          "--make-bed", "--out", paste(out, ".clean", sep=""),
                          args_filter),
           std_out=showPlinkOutput, std_err=showPlinkOutput)
    keepIDs <- data.table::fread(paste0(out, ".clean.fam"), data.table=FALSE,
                                     stringsAsFactors=FALSE,
                                     header=FALSE)
    keepIDs <- keepIDs[, 1:2]
    colnames(keepIDs) <- c("FID", "IID")
    colnames(removeIDs) <- c("FID", "IID")
    return(list(passIDs=keepIDs, failIDs=removeIDs))
}
