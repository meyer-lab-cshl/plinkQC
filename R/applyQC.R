#' Create plink dataset with samples and markers passing quality control
#'
#' Remove samples that fail per-sample QC and remove markers that fail
#' per-marker QC.
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param filterRelated [logical] Set to exclude samples that failed relatedness
#' check (via \code{\link{check_relatedness}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/alg.fail-IBD.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.evaluate_check_relatedness set to TRUE).
#' @param filterAncestry [logical] Set to exclude samples that failed ancestry
#' check (via \code{\link{check_ancestry}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/alg.fail-ancestry.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.check_ancestry set to TRUE).
#' @param filterHeterozygosity [logical] Set to exclude samples that failed
#' check for outlying heterozygosity rates (via
#' \code{\link{check_het_and_miss}} or
#' \code{\link{perIndividualQC}}). Requires file qcdir/alg.fail-het.IDs
#' (automatically created by \code{\link{perIndividualQC}} if
#' do.evaluate_check_het_and_miss set to TRUE).
#' @param filterSampleMissingness [logical] Set to exclude samples that failed
#' check for excessive missing genotype rates (via
#' \code{\link{check_het_and_miss}} or
#' \code{\link{perIndividualQC}}). Requires file qcdir/alg.fail-imiss.IDs
#' (automatically created by \code{\link{perIndividualQC}} if
#' do.evaluate_check_het_and_miss set to TRUE).
#' @param filterSex [logical] Set to exclude samples that failed the sex
#' check (via \code{\link{check_sex}} or \code{\link{perIndividualQC}}).
#' Requires file qcdir/alg.fail-sexcheck.IDs (automatically created by
#' \code{\link{perIndividualQC}} if do.evaluate_check_sex set to TRUE).
#' @param filterHWE [logical] Set to exclude markers that fail HWE exact test
#' (via \code{\link{check_hwe}} or \code{\link{perMarkerQC}}). Requires hweTh to
#' be set.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param filterMAF [logical] Set to exclude markers that fail minor allele
#' freuqency or minor allele count threshold (via \code{\link{check_maf}} or
#' \code{\link{perMarkerQC}}). Requires mafTh or macTh to be set.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specified, macTh is used (macTh = mafTh\*2\*NrSamples)
#' @param filterSNPMissingness [logical] Set to exclude markers that have
#' excessive missing rates across samples (via
#' \code{\link{check_snp_missingness}} or \code{\link{perMarkerQC}}). Requires
#' lmissTh to be set.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
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
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/2/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # All checks
#' ids_all <- cleanData(qcdir=qcdir, alg=alg, macTh=15, interactive=FALSE,
#' verbose=TRUE, path2plink=path2plink, filterAncestry=TRUE, filterRelated=TRUE)
#'
#' # Only filter samples failing check sex and heterozygosity/missingness
#' ids_sex_miss_het <- cleanData(qcdir=qcdir, alg=alg, macTh=15,
#' interactive=FALSE, verbose=TRUE, path2plink=path2plink, filterAncestry=FALSE,
#' filterRelated=FALSE)
#' }

cleanData <- function(qcdir, alg, lmissTh=0.01, hweTh=1e-5, macTh=20,
                       mafTh=0.01, filterAncestry=TRUE, filterRelated=TRUE,
                       filterSex=TRUE, filterHeterozygosity=TRUE,
                       filterSampleMissingness=TRUE,
                       filterSNPMissingness=TRUE, filterHWE=TRUE,
                       filterMAF=TRUE, path2plink=NULL, verbose=FALSE,
                       showPlinkOutput=TRUE) {
    sampleFilter <- c(filterAncestry, filterRelated, filterSex,
                      filterHeterozygosity, filterSampleMissingness)
    markerFilter <- c(filterHWE, filterMAF, filterSNPMissingness)
    prefix <- paste(qcdir,"/", alg, sep="")

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
    if (!file.exists(paste(prefix, ".fam",sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim",sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed",sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    # Remove remove.IDs file if already existing
    system(paste("rm -f ", prefix, ".remove.IDs", sep="" ))
    if (any(sampleFilter)) {
        if (filterRelated) {
            if (!file.exists(paste(prefix, ".fail-IBD.IDs", sep=""))){
                stop("filterRelated is TRUE but file ", prefix,
                     ".fail-IBD.IDs does not exist")
            } else {
                system(paste("cat ", prefix, ".fail-IBD.IDs >> ",
                             prefix, ".remove.IDs", sep=""))
            }
        }
        if (filterAncestry) {
            if (!file.exists(paste(prefix, ".fail-ancestry.IDs", sep=""))){
                stop("filterAncestry is TRUE but file ", prefix,
                     ".fail-ancestry.IDs does not exist")
            } else {
                system(paste("cat ", prefix, ".fail-ancestry.IDs >> ",
                             prefix, ".remove.IDs", sep=""))
            }
        }
        if (filterHeterozygosity) {
            if (!file.exists(paste(prefix, ".fail-het.IDs", sep=""))){
                stop("filterHeterozygosity is TRUE but file ", prefix,
                     ".fail-het.IDs does not exist")
            } else {
                system(paste("cat ", prefix, ".fail-het.IDs >> ",
                             prefix, ".remove.IDs", sep=""))
            }
        }
        if (filterSampleMissingness) {
            if (!file.exists(paste(prefix, ".fail-imiss.IDs", sep=""))){
                stop("filterSampleMissingness is TRUE but file ", prefix,
                     ".fail-imiss.IDs does not exist")
            } else {
                system(paste("cat ", prefix, ".fail-imiss.IDs >> ",
                             prefix, ".remove.IDs", sep=""))
            }
        }
        if (filterSex) {
            if (!file.exists(paste(prefix, ".fail-sexcheck.IDs", sep=""))){
                stop("filterSex is TRUE but file ", prefix,
                     ".fail-sexcheck.IDs does not exist")
            } else {
                system(paste("cat ", prefix, ".fail-sexcheck.IDs >> ",
                             prefix, ".remove.IDs", sep=""))
            }
        }
        # ensure unique IDs in remove.IDs
        removeIDs <- data.table::fread(paste(prefix, ".remove.IDs",
                                             sep=""),
                                       header=FALSE,
                                       data.table=FALSE, stringsAsFactors=FALSE)
        removeIDs <- removeIDs[!duplicated(removeIDs),]
        write.table(removeIDs, paste(prefix, ".remove.IDs", sep=""),
                    col.names=FALSE, row.names=FALSE, quote=FALSE)
        keepIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                     data.table=FALSE, stringsAsFactors=FALSE,
                                     header=FALSE)
        keepIDs <- keepIDs[keepIDs[,2] %in% removeIDs[,2],1:2]
        removeString <- paste(" --remove ", prefix, ".remove.IDs", sep="")
    } else {
        keepIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                     data.table=FALSE, stringsAsFactors=FALSE,
                                     header=FALSE)
        removeIDs <- NULL
        removeString <- ""
    }

    hwe <- NULL
    maf <- NULL
    missing <- NULL

    if (filterHWE) {
        if (is.null(hweTh)) {
            stop("filterHWE is TRUE but hweTh not specified")
        } else {
           hwe <- paste(" --hwe", hweTh)
        }
    }
    if (filterMAF) {
        if (is.null(mafTh) && is.null(macTh)) {
            stop("filterMAF is TRUE but neither mafTh or macTh are provided")
        }
        if(!is.null(macTh)) {
            all_samples <-  R.utils::countLines(paste(prefix, ".fam",
                                                      sep=""))
            keep_samples <- as.numeric(all_samples) - nrow(removeIDs)
            mafTh <- macTh/(2*keep_samples)
        }
        if (verbose) {
            if (!is.null(mafTh) && !is.null(macTh)) {
                message("Both mafTh and macTh provided, macTh=", macTh,
                        " is used (which corresponds to a mafTh=",
                        round(mafTh, 6), ")")
            } else if (!is.null(mafTh)) {
                message("The mafTh is ", mafTh)
            } else {
                message("The macTh is ", macTh," which corresponds to a mafTh=",
                        round(mafTh, 6), ")")
            }
        }
        maf <- paste(" --maf", mafTh)
    }
    if (filterSNPMissingness) {
        if (is.null(lmissTh)) {
            stop("filterSNPMissingness is TRUE but lmissTh not specified")
        } else {
            missing <- paste(" --geno", lmissTh)
        }
    }
    checkPlink(path2plink)
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }

    system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                 removeString,
                 maf, hwe, missing,
                 " --make-bed --out ", prefix, ".clean",
                 sep=""),
           ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
    return(list(passIDs=keepIDs, failIDs=removeIDs))
}
