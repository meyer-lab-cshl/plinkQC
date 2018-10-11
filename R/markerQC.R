#' Quality control for all markers in plink-dataset
#'
#' perMarkerQC checks the markers in the plink dataset for their missingness
#' rates across samples, their deviation from Hardy-Weinberg-Equilibrium (HWE)
#' and their minor allele frequencies (MAF). Per default, it assumes that IDs of
#' individuals that have failed \code{\link{perSampleQC}} have been written
#' to qcdir/alg.failIDs and removes these individuals when computing missingness
#' rates, HWE p-values and MAF. If the qcdir/alg.failIDs file does not exist, a
#' message is written to stdout but the analyses will continue for all samples
#' in the alg.fam/al.bed/alg.bim dataset.
#' Depicts i) SNP missingness rates (stratified by minor allele
#' frequency) as histograms, ii) p-values of HWE exact test (stratified by all
#' and low p-values) as histograms and iii) the minor allele frequency
#' distribution as a histogram.
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
#' and macTh are specified, macTh is used (macTh = mafTh\*2\*NrSamples)
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_marker) via ggplot2::ggsave(p=p_marker,
#' other_arguments) or pdf(outfile) print(p_marker) dev.off().
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return named [list] with i) fail_list, a named [list] with 1.
#' SNP_missingness, containing SNP IDs failing the missingness threshold
#' lmissTh, 2. hwe, containing SNP IDs failing the HWE exact test threshold
#' hweTh and 3. maf, containing SNPs failing the MAF threshold mafTh/MAC
#' threshold macTh and ii) p_markerQC, a ggplot2-object 'containing' a
#' sub-paneled plot with the QC-plots of \code{\link{check_snp_missingness}},
#' \code{\link{check_hwe}} and \code{\link{check_maf}}, which can
#' be shown by print(p_markerQC).
#' List entries contain NULL if that specific check was not chosen.
#' @details perMarkerQC wraps around the marker QC functions
#' \code{\link{check_snp_missingness}}, \code{\link{check_hwe}} and
#' \code{\link{check_maf}}. For details on the parameters and outputs, check
#' these function documentations.
#' @export
perMarkerQC <- function(qcdir, alg, mafTh=0.01, macTh=20, hweTh=1e-5,
                        lmissTh=0.01, interactive=FALSE, verbose=TRUE,
                        path2plink=NULL, showPlinkError=TRUE) {
    if (verbose) message("Identification of SNPs with high missingness rate")
    fail_snp_missingness <- check_snp_missingness(qcdir=qcdir, alg=alg,
                                                  lmissTh=lmissTh,
                                                  path2plink=path2plink,
                                                  verbose=verbose,
                                                  showPlinkError=showPlinkError)
    if (verbose) message("Identification of SNPs with deviation from HWE")
    fail_hwe <- check_hwe(qcdir=qcdir, alg=alg, hweTh=hweTh,
                          path2plink=path2plink, verbose=verbose,
                          showPlinkError=showPlinkError)
    if (verbose) message("Remove markers with a low minor allele frequency")
    fail_maf <- check_maf(qcdir=qcdir, alg=alg, mafTh=mafTh, macTh=macTh,
                          path2plink=path2plink, verbose=verbose,
                          showPlinkError=showPlinkError)
    fail_list <- list(SNP_missingness=fail_snp_missingness$fail_missingness$SNP,
                      hwe=fail_hwe$fail_hwe$SNP, maf=fail_maf$fail_maf$SNP)
    p_markerQC <- cowplot::plot_grid(fail_snp_missingness$p_lmiss,
                                   fail_hwe$p_hwe, fail_maf$p_maf,
                                   nrow=3, labels=LETTERS[1:3])
    if(interactive) print(p_markerQC)
    return(list(fail_list=fail_list, p_markerQC=p_markerQC))
}

#' Overview of per marker QC
#'
#' overviewPerMarkerQC depicts results of perMarkerQC as intersection plot (via
#' \code{\link[UpSetR]{upset}}) and returns dataframe indicating which QC
#' checks were failed or passed.
#'
#'
#' @param results_perMarkerQC [list] Output of \code{\link{perSampleQC}} i.e.
#' named [list] with i) fail_list, a named [list] with 1.
#' SNP_missingness, containing SNP IDs failing the missingness threshold
#' lmissTh, 2. hwe, containing SNP IDs failing the HWE exact test threshold
#' hweTh and 3. maf, containing SNPs failing the MAF threshold mafTh/MAC
#' threshold macTh and ii) p_markerQC, a ggplot2-object 'containing' a
#' sub-paneled plot with the QC-plots of \code{\link{check_snp_missingness}},
#' \code{\link{check_hwe}} and \code{\link{check_maf}}
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_overview) via ggplot2::ggsave(p=p_overview,
#' other_arguments) or pdf(outfile) print(p_overview) dev.off().
#' @return named [list] with i) nr_fail_markers: total number of markers
#' [integer] failing perMarkerQC, ii) fail_QC containing a [data.frame] with
#' markers that failed QC steps: marker rsIDs in rows,
#' columns are all QC steps applied by perMarkerQC (max=3), with entries=0 if
#' passing the QC and entries=1 if failing that particular QC and ii)
#' p_overview, a ggplot2-object 'containing' the overview plots of the QC
#' failures which can be shown by print(p_overview).
#' @export
overviewPerMarkerQC <- function(results_perMarkerQC, interactive=FALSE) {
    list2counts <- function(element, all_names) {
        all_names[!(all_names %in% element)] <- 0
        all_names[all_names %in% element] <- 1
        return(as.numeric(all_names))
    }
    if (length(perMarkerQC) == 2 &&
        !all(names(results_perMarkerQC) == c("fail_list", "p_markerQC"))) {
        stop("results_perMarkerQC not direct output of perMarkerQC")
    }
    fail_list <- results_perMarkerQC$fail_list
    # Remove null elements
    fail_list <- fail_list[!sapply(fail_list, is.null)]

    unique_samples_fail_all <- unique(unlist(fail_list))

    # overview QC fails
    unique_markers_fail <- unique(unlist(fail_list))
    fail_counts <- sapply(fail_list, list2counts, unique_markers_fail)
    rownames(fail_counts) <-  unique_markers_fail

    p_overview <- UpSetR::upset(UpSetR::fromList(fail_list),
                          order.by = "freq",
                          empty.intersections = "on", text.scale=1.2,
                          # Include when UpSetR v1.4.1 is released
                          # title="Overview quality control failures",
                          mainbar.y.label="Markers failing multiple QC checks",
                          sets.x.label="Marker fails per QC check",
                          main.bar.color="#1b9e77", matrix.color="#1b9e77",
                          sets.bar.color="#d95f02")
    if (interactive) {
        print(p_overview)
    }
    nr_fail_markers <- length(unique_markers_fail)
    return(list(nr_fail_samples=nr_fail_markers,
                fail_QC=fail_counts,
                p_overview=p_overview))
}




#' Identification of SNPs with high missingness rates
#'
#' Uses plink --remove alg.fail.IDs --missing --freq to calculate the rates of
#' missing genotype calls and frequency for all variants in the individuals that
#' passed the perSampleQC. Depicts SNP missingness rates (stratified by
#' minor allele frequency) as histograms.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_lmiss) via ggplot2::ggsave(p=p_lmiss,
#' other_arguments) or pdf(outfile) print(p_lmiss) dev.off().
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return named list with i) fail_missingness containing a [data.frame] with
#' CHR (Chromosome code), SNP (Variant identifier), CLST (Cluster identifier.
#' Only present with --within/--family), N_MISS (Number of missing genotype
#' call(s), not counting obligatory missings), N_CLST  (Cluster size; does not
#' include nonmales on Ychr; Only present with --within/--family), N_GENO
#' (Number of potentially valid call(s)), F_MISS (Missing call rate) for all
#' SNPs failing the lmissTh and ii) p_lmiss, a ggplot2-object 'containing' the
#' SNP missingness histogram which can be shown by (print(p_lmiss)).
#' @export
check_snp_missingness <- function(qcdir, alg, lmissTh=0.01, interactive=FALSE,
                                  path2plink=NULL, verbose=FALSE,
                                  showPlinkError=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    checkPlink(path2plink)
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_SNP_missingness for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!verbose, ignore.stderr=!showPlinkError)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=verbose, ignore.stderr=!showPlinkError)
    }

    lmiss <- read.table(paste(qcdir, "/", alg, suffix, ".lmiss",sep=""),
                       as.is=TRUE, header=TRUE)

    frq <- read.table(paste(qcdir, "/", alg, suffix, ".frq", sep=""),
                  header=TRUE, as.is=TRUE)
    lmiss_frq <- merge(lmiss, frq)
    lmiss_frq$MAF_bin <- ifelse(lmiss_frq$MAF < 0.05, 1, 0)
    p_highMAF <- ggplot(dplyr::filter_(lmiss_frq, ~MAF_bin == 0),
                        aes_string('F_MISS'))
    p_highMAF <- p_highMAF + geom_histogram(binwidth = 0.005,
                                            fill="#66a61e") +
        ylab("Number of SNPs") +
        xlab("Proportion of missing data") +
        ggtitle("SNPs with MAF > 0.05") +
        geom_vline(xintercept=lmissTh, lty=2, col="red") +
        theme_bw() +
        theme(title=element_text(size=10))
    p_lowMAF <- ggplot(dplyr::filter_(lmiss_frq, ~MAF_bin == 1),
                       aes_string('F_MISS'))
    p_lowMAF <- p_lowMAF + geom_histogram(binwidth = 0.005,
                                          fill="#e6ab02") +
        ylab("Number of SNPs") +
        xlab("Proportion of missing data") +
        ggtitle("SNPs with MAF < 0.05") +
        geom_vline(xintercept=lmissTh, lty=2, col="red") +
        theme_bw() +
        theme(title=element_text(size=10))
    p_histo <- cowplot::plot_grid(p_lowMAF, p_highMAF)
    title <- cowplot::ggdraw() +
        cowplot::draw_label("Marker missingness rate")
    p_lmiss <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                rel_heights = c(0.1, 1))
    if (interactive) {
        print(p_lmiss)
    }
    fail_missingness <- lmiss[lmiss$F_MISS > lmissTh,]
    return(list(fail_missingness=fail_missingness, p_lmiss=p_lmiss))
}

#' Identification of SNPs showing a significant deviation from Hardy-Weinberg-
#' equilibrium (HWE)
#'
#' Uses plink --remove alg.fail.IDs --hwe to calculate the observed and expected
#' heterozygote frequencies per SNP in the individuals that
#' passed the perSampleQC and compute the deviation of the frequencies from
#' Hardy-Weinberg equilibrium (HWE) by HWE exact test.
#' The p-values of the HWE exact test are displayed as histograms (stratified by
#' all and low p-values), where the hweTh is used to depict the QC cut-off for
#' SNPs.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_hwe) via ggplot2::ggsave(p=p_hwe,
#' other_arguments) or pdf(outfile) print(p_hwe) dev.off().
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return named list with i) fail_hwe containing a [data.frame] with CHR
#' (Chromosome code), SNP (Variant identifier), TEST (Type of test: one of
#' {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}), A1 (Allele 1; usually minor),
#' A2 (Allele 2; usually major), GENO ('/'-separated genotype counts: A1 hom,
#' het, A2 hom), O(HET) (Observed heterozygote frequency E(HET) (Expected
#' heterozygote frequency), P (Hardy-Weinberg equilibrium exact test p-value)
#' for all SNPs that failed the hweTh and ii) p_hwe, a ggplot2-object
#' 'containing' the HWE p-value distribution histogram which can be shown by
#' (print(p_hwe)).
#' @export
check_hwe <- function(qcdir, alg, hweTh=1e-5, interactive=FALSE,
                      path2plink=NULL, verbose=FALSE, showPlinkError=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    checkPlink(path2plink)
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perSampleQC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_HWE for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --hardy",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!verbose, ignore.stderr=!showPlinkError)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir, "/", alg, ".fail.IDs --hardy",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!verbose, ignore.stderr=!showPlinkError)
    }
    hwe <- read.table(paste(qcdir, "/", alg, suffix, ".hwe", sep=""),
                       header=TRUE, as.is=TRUE)
    hwe <- hwe[grepl("ALL", hwe$TEST),]
    hwe$P_bin <- ifelse(hwe$P < 0.01, 1, 0)
    hwe$minus_log10P <- -log10(hwe$P)
    p_allP <- ggplot(hwe, aes_string('minus_log10P'))
    p_allP <- p_allP + geom_histogram(binwidth = 0.5,
                                      fill="#66a61e") +
        ylab("Number of SNPs") +
        xlab(expression(-log[10](HWE~exact~test~p-value))) +
        ggtitle(expression(All~p-value[HWE])) +
        geom_vline(xintercept=-log10(hweTh), lty=2, col="red") +
        theme_bw() +
        theme(title=element_text(size=10))
    p_lowP <- ggplot(dplyr::filter_(hwe, ~P_bin == 1),
                     aes_string('minus_log10P'))
    p_lowP <- p_lowP + geom_histogram(binwidth = 0.5,
                                            fill="#e6ab02") +
        ylab("Number of SNPs") +
        xlab(expression(-log[10](HWE~exact~test~p-value))) +
        ggtitle(expression(p-value[HWE]<0.01)) +
        geom_vline(xintercept=-log10(hweTh), lty=2, col="red") +
        theme_bw() +
        theme(title=element_text(size=10))
    p_histo <- cowplot::plot_grid(p_allP, p_lowP)
    title <- cowplot::ggdraw() +
        cowplot::draw_label(expression(Distribution~of~-log[10](p-value[HWE])))
    p_hwe <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                  rel_heights = c(0.1, 1))
    if (interactive) {
        print(p_hwe)
    }
    fail_hwe <- hwe[hwe$P < hweTh,]
    return(list(fail_hwe=fail_hwe, p_hwe=p_hwe))
}

#' Identification of SNPs with a very low minor allele frequency
#'
#' Uses plink --remove alg.fail.IDs --freq to calculate the minor allele
#' frequencies for all variants in the individuals that passed the perSampleQC.
#' The minor allele frequency distributions is plotted as a histogram.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specied, macTh is used (macTh = mafTh\*2\*NrSamples)
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_hwe) via ggplot2::ggsave(p=p_maf,
#' other_arguments) or pdf(outfile) print(p_maf) dev.off().
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return named list with i) fail_maf containing a [data.frame] with CHR
#' (Chromosome code), SNP (Variant identifier), A1 (Allele 1; usually minor), A2
#' (Allele 2; usually major), MAF (Allele 1 frequency), NCHROBS (Number of
#' allele observations) for all SNPs that failed the mafTh/macTh and ii) p_maf,
#' a ggplot2-object 'containing' the MAF distribution histogram which can be
#' shown by (print(p_maf)).
#' @export
check_maf <- function(qcdir, alg, mafTh=0.01, macTh=20, verbose=FALSE,
                      interactive=FALSE, path2plink=NULL, showPlinkError=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    checkPlink(path2plink)
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    if (!file.exists(paste(qcdir, "/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perSampleQC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_maf for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --freq ",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!verbose, ignore.stderr=!showPlinkError)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --freq",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!verbose, ignore.stderr=!showPlinkError)
    }
    maf <- read.table(paste(qcdir,"/", alg, suffix, ".frq",sep=""),
                       header=TRUE, as.is=TRUE)

    all_samples <-  R.utils::countLines(paste(qcdir,"/",alg, ".fam",sep=""))
    fail_samples <-  R.utils::countLines(paste(qcdir,"/",alg, ".fail.IDs",
                                               sep=""))
    keep_samples <- all_samples - fail_samples

    if (is.null(mafTh) && is.null(macTh)) {
        stop("Either mafTh or macTh need to be provided")
    }
    if(!is.null(macTh)) {
        mafTh <- macTh/(2*keep_samples)
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
    p_maf <- ggplot(maf, aes_string('MAF'))
    p_maf <- p_maf + geom_histogram(binwidth = 0.01,
                                            fill="#999999") +
        ylab("Number of SNPs") +
        xlab("Minor allele frequency") +
        ggtitle("Minor allele frequency distribution") +
        geom_vline(xintercept=mafTh, lty=2, col="red") +
        theme_bw() +
        theme(title=element_text(size=10))
    if (interactive) {
        print(p_maf)
    }
    fail_maf <- maf[maf$MAF < mafTh,]
    return(list(fail_maf=fail_maf, p_maf=p_maf))
}
