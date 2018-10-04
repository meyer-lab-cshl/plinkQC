#' Quality control for all markers in plink-dataset
#'
#' perMarkerQC checks the markers in the plink dataset for their missingness
#' rates across samples, their deviation from Hardy-Weinberg-Equilibrium (HWE)
#' and their minor allele frequencies (MAF). Per default, it assumes that IDs of
#' individuals that have failed \link{\code{perIndividualQC}}have been written to
#' qcdir/alg.failIDs and removes these individuals when computing missingness
#' rates, HWE p-values and MAF. If the qcdir/alg.failIDs file does not exist, a
#' message is written to stdout but the analyses will continue for all samples
#' in the alg.fam/al.bed/alg.bim dataset.
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specified, macTh is used (macTh = mafTh*2*NrSamples)
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) perMarkerQC(arguments) dev.off().
#' If TRUE, depicts i) SNP missingness rates (stratified by minor allele
#' frequency) as histograms, ii) p-values of HWE exact test (stratified by all
#' and low p-values) as histograms and iii) the minor allele frequency
#' distribution as a histogram.
#' @return named [list] with i) SNP_missingness containing a [vector] with
#' SNP IDs failing the missingness threshold lmissTh, ii) hew containing a
#' [vector] with SNP IDs failing the HWE exact test threshold hweTh and iii)
#' maf containing a [vector] with SNPs failing the MAF threshold mafTh.
#' @details perMarkerQC wraps around the marker QC functions
#' \link{\code{check_snp_missingness}}, \link{\code{check_hwe}} and
#' \link{\code{check_maf}}. For details on the parameters and outputs, check
#' these function documentations.
perMarkerQC <- function(qcdir, alg, mafTh=0.01, macTh=20, hweTh=1e-5,
                        lmissTh=0.01, plot=FALSE, verbose=TRUE) {
    if (verbose) message("Identification of SNPs with high missingness rate")
    fail_snp_missingness <- check_snp_missingness(qcdir=qcdir, alg=alg,
                                                  lmissTh=lmissTh, plot=plot)
    if (verbose) message("Identification of SNPs with deviation from HWE")
    fail_hwe <- check_hwe(qcdir=qcdir, alg=alg, hweTh=hweTh, plot=plot)
    if (verbose) message("Remove markers with a low minor allele frequency")
    fail_maf <- check_maf(qcdir=qcdir, alg=alg, mafTh=mafTh, macTh=macTh,
                          plot=plot)
    fail_list <- list(SNP_missingness=fail_snp_missingness$SNP,
                      hwe=fail_hwe$SNP, maf=fail_maf$SNP)
    return(fail_list)
}

#' Identification of SNPs with high missingness rates
#'
#' Uses plink --remove alg.fail.IDs --missing --freq to calculate the rates of
#' missing genotype calls and frequency for all variants in the individuals that
#' passed the perIndividualQC.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_snp_missingness(arguments) dev.off().
#' If TRUE, depicts SNP missingness rates (stratified by minor allele frequency)
#' as histograms.
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier),
#' CLST (Cluster identifier. Only present with --within/--family), N_MISS
#' (Number of missing genotype call(s), not counting obligatory missings),
#' N_CLST  (Cluster size; does not include nonmales on Ychr; Only present with
#' --within/--family), N_GENO (Number of potentially valid call(s)), F_MISS
#' (Missing call rate) for all SNPs failing the lmissTh.
check_snp_missingness <- function(qcdir, alg, lmissTh=0.01, plot=FALSE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_SNP_missingness for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --bfile ", qcdir, "/", alg,
                     "--missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }

    lmiss <- read.table(paste(qcdir, "/", alg, suffix, ".lmiss",sep=""),
                       as.is=TRUE, header=TRUE)
    if (plot) {
        frq <- read.table(paste(qcdir, "/", alg, suffix, ".frq", sep=""),
                      header=TRUE, as.is=TRUE)
        lmiss_frq <- merge(lmiss, frq)
        dataLowMAF <- lmiss_frq[lmiss_frq$MAF < 0.05,]
        hist(log10(dataLowMAF$F_MISS), axes=F, xlim=c(-4,0), col="burlywood",
             ylab="Number of SNPs", xlab="Proportion of missing data",
             main="Missingness for SNPs with MAF<0.05")
        axis(side=2, labels=TRUE)
        axis(side=1, labels=FALSE)
        mtext(c("0.0001","0.001", lmissTh, "0.01", "0.1", "1"), side=1,
              at=c(-4, -3, log10(lmissTh), -2, -1, 0), line=1)
        abline(v=log10(lmissTh), lty=2, col="red")
        dataHighMAF <- lmiss_frq[lmiss_frq$MAF >=0.05,]
        hist(log10(dataHighMAF$F_MISS), axes=FALSE, xlim=c(-4,0),
             col="burlywood",
             ylab="Number of SNPs", xlab="Proportion of missing data",
             main="Missingness for SNPs with MAF>=0.05")
        axis(side=2, labels=TRUE)
        axis(side=1, labels=FALSE)
        mtext(c("0.0001","0.001", lmissTh, "0.01", "0.1", "1"), side=1,
              at=c(-4, -3, log10(lmissTh), -2, -1, 0), line=1)
        abline(v=log10(lmissTh), lty=2, col="red")
    }
    fail_missingness <- lmiss[lmiss$F_MISS > lmissTh,]
    return(fail_missingness)
}

#' Identification of SNPs showing a significant deviation from Hardy-Weinberg-
#' equilibrium (HWE)
#'
#' Uses plink --remove alg.fail.IDs --hwe to calculate the observed and expected
#' heterozygote frequencies per SNP in the individuals that
#' passed the perIndividualQC and compute the deviation of the frequencies from
#' Hardy-Weinberg equilibrium (HWE) by HWE exact test.
#' If plot is TRUE, the p-values of the HWE exact test are displayed at
#' histograms, where the hweTh is used to depict the QC cut-off for SNPs.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_hwe(arguments) dev.off().
#' If TRUE, depicts p-values of HWE exact test (stratified by all and low
#' p-values) as histograms.
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier),
#' TEST (Type of test: one of {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}), A1
#' (Allele 1; usually minor), A2 (Allele 2; usually major), GENO ('/'-separated
#' genotype counts: A1 hom, het, A2 hom), O(HET) (Observed heterozygote frequency
#' E(HET) (Expected heterozygote frequency), P (Hardy-Weinberg equilibrium exact
#' test p-value) for all SNPs that failed the hweTh.
check_hwe <- function(qcdir, alg, hweTh=1e-5, plot=FALSE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_HWE for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --bfile ", qcdir, "/", alg,
                     "--hardy ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir, "/", alg, ".fail.IDs --hardy ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }
    hwe <- read.table(paste(qcdir, "/", alg, suffix, ".hwe", sep=""),
                       header=TRUE, as.is=TRUE)
    hwe <- hwe[grepl("ALL", hwe$TEST),]
    if (plot) {
        hist(-log10(hwe$P), ylab="Number of SNPs",col="palegreen",
             xlab="-log(P_hwe)", main="HWE_Pvalue for All SNPs")
        hwe <- hwe[hwe$P <0.01,]
        hist(-log10(hwe$P), ylab="Number of SNPs", col="palegreen",
             xlab="-log(P_hwe)",
             main="HWE_Pvalue for SNPs with HWE_Pvalue <0.01")
        abline(v=-log10(hweTh), lty=2, col="red")
    }
    fail_hwe <- hwe[hwe$P < hweTh,]
    return(fail_hwe)
}

#' Identification of SNPs with a very low minor allele frequency
#'
#' Uses plink --remove alg.fail.IDs --freq to calculate the minor allele
#' frequencies for all variants in the individuals that passed the perSampleQC.
#' If plot is TRUE, the minor allele frequency distributions is plotted as a
#' histogram.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specied, macTh is used (macTh = mafTh*2*NrSamples)
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_maf(arguments) dev.off().
#' If TRUE, the minor allele frequency distribution is plotted in a histogram.
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier), A1
#' (Allele 1; usually minor), A2 (Allele 2; usually major), MAF (Allele 1
#' frequency), NCHROBS (Number of allele observations) for all SNPs that failed
#' the mafTh/macTh.
check_maf <- function(qcdir, alg, mafTh=0.01, macTh=20, verbose=FALSE,
                      plot=FALSE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_maf for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --noweb --bfile ", qcdir, "/", alg,
                     " --freq --make-bed ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --noweb --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --freq ",
                     "--make-bed ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }
    maf <- read.table(paste(qcdir,"/",alg, suffix, ".frq",sep=""),
                       header=TRUE, as.is=TRUE)
    samples <- read.table(paste(qcdir,"/",alg, suffix, ".fam",sep=""),
                          header=TRUE, as.is=TRUE)
    if (is.null(mafTh) && is.null(macTh)) {
        stop("Either mafTh or macTh need to be provided")
    }
    if(!is.null(macTh)) {
        mafTh <- macTh/(2*nrow(samples))
    }
    if (verbose) {
        if (!is.null(mafTh) && !is.null(macTh)) {
            message("Both mafTh and macTh provided, macTh=", macTh,
                    " is used (corresponds to mafTh=", mafTh, ")")
        } else if (!is.null(mafTh)) {
            message("The mafTh is ", mafTh)
        } else {
            message("The macTh is ", macTh," which corresponds to a mafTh=",
                    mafTh, ")")
        }
    }
    if (plot) {
        hist(maf$MAF, col="lightblue1", ylab="Number of SNPs", breaks=50,
             xlab="Minor Allele Frequency",
             main="Minor Allele Frequency (MAF) for All SNPs")
        abline(v=mafTh, lty=2, col="red")
    }
    failed_maf <- maf[maf$MAF < mafTh,]
    return(failed_maf)
}
