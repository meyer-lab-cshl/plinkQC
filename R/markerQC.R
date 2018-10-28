#' Quality control for all markers in plink-dataset
#'
#' perMarkerQC checks the markers in the plink dataset for their missingness
#' rates across samples, their deviation from Hardy-Weinberg-Equilibrium (HWE)
#' and their minor allele frequencies (MAF). Per default, it assumes that IDs of
#' individuals that have failed \code{\link{perIndividualQC}} have been written
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
#' @param alg [character] prefix of plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param do.check_hwe [logical] If TRUE, run \code{\link{check_hwe}}
#' @param do.check_maf [logical] If TRUE, run \code{\link{check_maf}}.
#' @param do.check_snp_missingness [logical] If TRUE, run
#' \code{\link{check_snp_missingness}}.
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
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @return Named [list] with i) fail_list, a named [list] with 1.
#' SNP_missingness, containing SNP IDs [vector] failing the missingness
#' threshold lmissTh, 2. hwe, containing SNP IDs [vector] failing the HWE exact
#' test threshold hweTh and 3. maf, containing SNPs Ids [vector] failing the MAF
#' threshold mafTh/MAC threshold macTh and ii) p_markerQC, a ggplot2-object
#' 'containing' a sub-paneled plot with the QC-plots of
#' \code{\link{check_snp_missingness}}, \code{\link{check_hwe}} and
#' \code{\link{check_maf}}, which can be shown by print(p_markerQC).
#' List entries contain NULL if that specific check was not chosen.
#' @details perMarkerQC wraps around the marker QC functions
#' \code{\link{check_snp_missingness}}, \code{\link{check_hwe}} and
#' \code{\link{check_maf}}. For details on the parameters and outputs, check
#' these function documentations.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' # All quality control checks
#' \dontrun{
#' fail_markers <- perMarkerQC(qcdir=qcdir, alg=alg, interactive=FALSE,
#' verbose=TRUE, path2plink=path2plink)
#' }
perMarkerQC <- function(qcdir, alg,
                        do.check_snp_missingness=TRUE, lmissTh=0.01,
                        do.check_hwe=TRUE, hweTh=1e-5,
                        do.check_maf=TRUE, mafTh=0.01, macTh=20,
                        interactive=FALSE, verbose=TRUE,
                        path2plink=NULL, showPlinkOutput=TRUE
                        ) {
    if (do.check_snp_missingness) {
        if (verbose) message("Identification of SNPs with high missingness",
                             " rate")
        fail_snp_missingness <- check_snp_missingness(qcdir=qcdir, alg=alg,
                                                      lmissTh=lmissTh,
                                                      path2plink=path2plink,
                                                      verbose=verbose,
                                                      showPlinkOutput=
                                                        showPlinkOutput)
        p_missingness <- fail_snp_missingness$p_lmiss
    } else {
        fail_snp_missingness <- NULL
        p_missingness <- NULL
    }
    if (do.check_hwe) {
        if (verbose) message("Identification of SNPs with deviation from HWE")
        fail_hwe <- check_hwe(qcdir=qcdir, alg=alg, hweTh=hweTh,
                              path2plink=path2plink, verbose=verbose,
                              showPlinkOutput=showPlinkOutput)
        p_hwe <- fail_hwe$p_hwe
    } else {
        fail_hwe <- NULL
        p_hwe <- NULL
    }
    if (do.check_maf) {
        if (verbose) message("Remove markers with a low minor allele frequency")
        fail_maf <- check_maf(qcdir=qcdir, alg=alg, mafTh=mafTh, macTh=macTh,
                              path2plink=path2plink, verbose=verbose,
                              showPlinkOutput=showPlinkOutput)
        p_maf <- fail_maf$p_maf
    } else {
        fail_maf <- NULL
        p_maf <- NULL
    }
    fail_list <- list(SNP_missingness=fail_snp_missingness$fail_missingness$SNP,
                      hwe=fail_hwe$fail_hwe$SNP, maf=fail_maf$fail_maf$SNP)

    plots_markerQC <- list(p_missingness, p_hwe, p_maf)
    plots_markerQC <- plots_markerQC[sapply(plots_markerQC,
                                            function(x) !is.null(x))]
    subplotLabels <- LETTERS[1:length(plots_markerQC)]
    p_markerQC <- cowplot::plot_grid(plotlist=plots_markerQC,
                                     nrow=length(plots_markerQC),
                                     labels=subplotLabels,
                                     rel_heights=c(rep(1,
                                                       length(plots_markerQC)),
                                                       1.5))
    if(interactive) print(p_markerQC)
    return(list(fail_list=fail_list, p_markerQC=p_markerQC))
}

#' Overview of per marker QC
#'
#' overviewPerMarkerQC depicts results of perMarkerQC as intersection plot (via
#' \code{\link[UpSetR]{upset}}) and returns dataframe indicating which QC
#' checks were failed or passed.
#'
#' @param results_perMarkerQC [list] Output of \code{\link{perIndividualQC}} i.e.
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
#' @return Named [list] with i) nr_fail_markers: total number of markers
#' [integer] failing perMarkerQC, ii) fail_QC containing a [data.frame] with
#' markers that failed QC steps: marker rsIDs in rows,
#' columns are all QC steps applied by perMarkerQC (max=3), with entries=0 if
#' passing the QC and entries=1 if failing that particular QC.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' # All quality control checks
#' \dontrun{
#' fail_markers <- perMarkerQC(qcdir=qcdir, alg=alg, interactive=FALSE,
#' verbose=TRUE, path2plink=path2plink)
#' overview <- overviewPerMarkerQC(fail_markers)
#' }
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

    if (interactive) {
        if (length(fail_list) >= 2) {
            UpSetR::upset(UpSetR::fromList(fail_list),
                      order.by = "freq",
                      empty.intersections = "on", text.scale=1.2,
                      # Include when UpSetR v1.4.1 is released
                      # title="Overview quality control failures",
                      mainbar.y.label="Markers failing multiple QC checks",
                      sets.x.label="Marker fails per QC check",
                      main.bar.color="#1b9e77", matrix.color="#1b9e77",
                      sets.bar.color="#d95f02")
        } else {
            message("overviewMarkerQC cannot be displayed with UpSetR: at ",
                    "least two elements in list required, but only ",
                    length(fail_list) ," provided")
        }
        #print(p_overview)
    }
    nr_fail_markers <- length(unique_markers_fail)
    return(list(nr_fail_samples=nr_fail_markers,
                fail_QC=fail_counts))
                #p_overview=p_overview))
}

#' Identification of SNPs with high missingness rate
#'
#' Runs and evaluates results from plink --missing --freq. It calculate the
#' rates of missing genotype calls and frequency for all variants in the
#' individuals that passed the \code{\link{perIndividualQC}}. The SNP
#' missingness rates (stratified by minor allele frequency) as histograms.
#'
#' \code{check_snp_missingness} uses plink --remove alg.fail.IDs --missing
#' --freq to calculate rates of missing genotype calls and frequency per SNP in
#' the individuals that passed the \code{\link{perIndividualQC}}. It does so
#' without generating a new dataset but simply removes the IDs when calculating
#' the statistics.
#'
#' For details on the output data.frame fail_missingness, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#lmiss}.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] Prefix of plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_lmiss) via ggplot2::ggsave(p=p_lmiss,
#' other_arguments) or pdf(outfile) print(p_lmiss) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return Named list with i) fail_missingness containing a [data.frame] with
#' CHR (Chromosome code), SNP (Variant identifier), CLST (Cluster identifier.
#' Only present with --within/--family), N_MISS (Number of missing genotype
#' call(s), not counting obligatory missings), N_CLST  (Cluster size; does not
#' include nonmales on Ychr; Only present with --within/--family), N_GENO
#' (Number of potentially valid call(s)), F_MISS (Missing call rate) for all
#' SNPs failing the lmissTh and ii) p_lmiss, a ggplot2-object 'containing' the
#' SNP missingness histogram which can be shown by (print(p_lmiss)).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' fail_snp_missingness <- check_snp_missingness(qcdir=qcdir, alg=alg,
#' interactive=FALSE, verbose=TRUE, path2plink=path2plink)
#' }
check_snp_missingness <- function(qcdir, alg, lmissTh=0.01, interactive=FALSE,
                                  path2plink=NULL, verbose=FALSE,
                                  showPlinkOutput=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
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
               ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
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
#' Runs and evaluates results from plink --hardy. It calculates the observed and
#' expected heterozygote frequencies for all variants in the individuals that
#' passed the \code{\link{perIndividualQC}} and computes the deviation of the
#' frequencies from Hardy-Weinberg equilibrium (HWE) by HWE exact test. The
#' p-values of the HWE exact test are displayed as histograms (stratified by
#' all and low p-values), where the hweTh is used to depict the quality control
#' cut-off for SNPs.
#'
#' \code{check_hwe} uses plink --remove alg.fail.IDs --hardy to
#' calculate the observed and expected heterozygote frequencies per SNP in the
#' individuals that passed the \code{\link{perIndividualQC}}. It does so
#' without generating a new dataset but simply removes the IDs when calculating
#' the statistics.
#'
#' For details on the output data.frame fail_hwe, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#hwe}.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] Prefix of plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_hwe) via ggplot2::ggsave(p=p_hwe,
#' other_arguments) or pdf(outfile) print(p_hwe) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return Named list with i) fail_hwe containing a [data.frame] with CHR
#' (Chromosome code), SNP (Variant identifier), TEST (Type of test: one of
#' {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}), A1 (Allele 1; usually minor),
#' A2 (Allele 2; usually major), GENO ('/'-separated genotype counts: A1 hom,
#' het, A2 hom), O(HET) (Observed heterozygote frequency E(HET) (Expected
#' heterozygote frequency), P (Hardy-Weinberg equilibrium exact test p-value)
#' for all SNPs that failed the hweTh and ii) p_hwe, a ggplot2-object
#' 'containing' the HWE p-value distribution histogram which can be shown by
#' (print(p_hwe)).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' fail_hwe <- check_hwe(qcdir=qcdir, alg=alg, interactive=FALSE, verbose=TRUE,
#' path2plink=path2plink)
#' }
check_hwe <- function(qcdir, alg, hweTh=1e-5, interactive=FALSE,
                      path2plink=NULL, verbose=FALSE, showPlinkOutput=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (!file.exists(paste(qcdir,"/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividualQC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_HWE for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --hardy",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir, "/", alg, ".fail.IDs --hardy",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
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

#' Identification of SNPs with low minor allele frequency
#'
#' Runs and evaluates results from plink --freq. It calculates the minor allele
#' frequencies for all variants in the individuals that passed the
#' \code{\link{perIndividualQC}}. The minor allele frequency distributions is
#' plotted as a histogram.
#'
#' \code{check_maf} uses plink --remove alg.fail.IDs --freq to calculate the
#' minor allele frequencies for all variants in the individuals that passed the
#' \code{\link{perIndividualQC}}. It does so without generating a new dataset
#' but simply removes the IDs when calculating the statistics.
#'
#' For details on the output data.frame fail_maf, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#frq}.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] Prefix of plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specied, macTh is used (macTh = mafTh\*2\*NrSamples)
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_hwe) via ggplot2::ggsave(p=p_maf,
#' other_arguments) or pdf(outfile) print(p_maf) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out
#' and specifically, if TRUE, plink log will be displayed.
#' @return Named list with i) fail_maf containing a [data.frame] with CHR
#' (Chromosome code), SNP (Variant identifier), A1 (Allele 1; usually minor), A2
#' (Allele 2; usually major), MAF (Allele 1 frequency), NCHROBS (Number of
#' allele observations) for all SNPs that failed the mafTh/macTh and ii) p_maf,
#' a ggplot2-object 'containing' the MAF distribution histogram which can be
#' shown by (print(p_maf)).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' alg <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' fail_maf <- check_maf(qcdir=qcdir, alg=alg, macTh=15, interactive=FALSE,
#' verbose=TRUE, path2plink=path2plink)
#' }
check_maf <- function(qcdir, alg, mafTh=0.01, macTh=20, verbose=FALSE,
                      interactive=FALSE, path2plink=NULL,
                      showPlinkOutput=TRUE) {
    if (!file.exists(paste(qcdir,"/", alg, ".fam",sep=""))){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bim",sep=""))){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(paste(qcdir,"/", alg, ".bed",sep=""))){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (!file.exists(paste(qcdir, "/", alg, ".fail.IDs",sep=""))){
        message("File with individuals that failed perIndividualQC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_maf for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --freq ",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
    } else {
        suffix <- ".no_failIDs"
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --freq",
                     " --out ", qcdir, "/", alg, suffix, sep=""),
               ignore.stdout=!showPlinkOutput, ignore.stderr=!showPlinkOutput)
    }
    maf <- read.table(paste(qcdir,"/", alg, suffix, ".frq",sep=""),
                       header=TRUE, as.is=TRUE)

    all_samples <-  R.utils::countLines(paste(qcdir,"/",alg, ".fam",sep=""))
    fail_samples <-  R.utils::countLines(paste(qcdir,"/",alg, ".fail.IDs",
                                               sep=""))
    keep_samples <- as.numeric(all_samples - fail_samples)

    if (is.null(mafTh) && is.null(macTh)) {
        stop("Either mafTh or macTh need to be provided")
    }
    if(!is.null(macTh)) {
        mafTh <- macTh/(2*keep_samples)
    }
    if (verbose) {
        if (!is.null(mafTh) && !is.null(macTh)) {
            message("Both mafTh and macTh provided, macTh=", macTh,
                    " is used (corresponds to mafTh=", round(mafTh, 6), ")")
        } else if (!is.null(mafTh)) {
            message("The mafTh is ", round(mafTh, 6))
        } else {
            message("The macTh is ", macTh," which corresponds to a mafTh=",
                    round(mafTh, 6), ")")
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
