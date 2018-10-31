#' Quality control for all individuals in plink-dataset
#'
#' perIndividualQC checks the samples in the plink dataset for their total
#' missingness and heterozygosity rates, the concordance of their assigned sex
#' to their SNP sex, their relatedness to other study individuals and their
#' genetic ancestry.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where results will be saved.
#' Per default, qcdir=indir. If do.evaluate_[analysis] is set to TRUE and
#' do.run_[analysis] is FALSE, \code{perIndividualQC} expects the
#' analysis-specific plink output files in qcdir i.e. do.check_sex expects
#' name.sexcheck, do.evaluate_check_het_and_miss expects name.het and name.imiss,
#' do.evauluate_check_relatedness expects name.genome and name.imiss and
#' do.evaluate_check_ancestry expects prefixMergeData.eigenvec. If these files
#' are not present \code{perIndividualQC} will fail with missing file error.
#' Setting do.run_[analysis] TRUE will execute the checks and create the
#' required files. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param dont.check_sex [logical] If TRUE, no sex check will be conducted;
#' short for do.run_check_sex=FALSE and do.evaluate_check_sex=FALSE.
#' Takes precendence over do.run_check_sex and do.evaluate_check_sex.
#' @param dont.check_het_and_miss [logical] If TRUE, no heterozygosity and
#' missingness check will be conducted; short for
#' do.run_check_heterozygosity=FALSE, do.run_check_missingness=FALSE and
#' do.evaluate_check_het_and_miss=FALSE. Takes precendence over
#' do.run_check_heterozygosity, do.run_check_missingness and
#' do.evaluate_check_het_and_miss.
#' @param dont.check_relatedness [logical] If TRUE, no relatedness check will be
#' conducted; short for do.run_check_relatedness=FALSE and
#' do.evaluate_check_relatedness=FALSE. Takes precendence over
#' do.run_check_relatedness and do.evaluate_check_relatedness.
#' @param dont.check_ancestry [logical] If TRUE, no ancestry check will be
#' conducted; short for do.run_check_ancestry=FALSE and
#' do.evaluate_check_ancestry=FALSE. Takes precendence over
#' do.run_check_ancestry and do.evaluate_check_ancestry.
#' @param do.run_check_sex [logical] If TRUE, run \code{\link{run_check_sex}}
#' @param do.run_check_het_and_miss [logical] If TRUE, run
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}}
#' @param do.run_check_relatedness [logical] If TRUE, run
#' \code{\link{run_check_relatedness}}.
#' @param do.run_check_ancestry [logical] If TRUE, run
#' \code{\link{run_check_ancestry}}.
#' @param do.evaluate_check_sex [logical] If TRUE, run
#' \code{\link{evaluate_check_sex}}
#' @param do.evaluate_check_het_and_miss [logical] If TRUE, run
#' \code{\link{evaluate_check_het_and_miss}}.
#' @param do.evaluate_check_relatedness [logical] If TRUE, run
#' \code{\link{evaluate_check_relatedness}}.
#' @param do.evaluate_check_ancestry [logical] If TRUE, run
#' \code{\link{evaluate_check_ancestry}}.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param externalSex [data.frame, optional] Dataframe with sample IDs
#' [externalSexID] and sex [externalSexSex] to double check if external and
#' PEDSEX data (often processed at different centers) match.
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @param externalSexSex [character] Column identifier for column containing sex
#' information in externalSex.
#' @param externalSexID [character] Column identifier for column containing ID
#' information in externalSex.
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex while Sex==SNPSEX automatically corrected: this will directly
#' change the name.bim/.bed/.fam files!
#' @param imissTh [double] Threshold for acceptable missing genotype rate per
#' individual.
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity per individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param prefixMergedDataset [character] Prefix of merged dataset (study and
#' reference samples) with prefixMergedDataset.bed, prefixMergedDataset.bim
#' and prefixMergedDataset.fam, which was used in plink --pca, resulting in
#' prefixMergedDataset.eigenvec
#' @param europeanTh [double] Scaling factor of radius to be drawn around center
#' of reference European samples, with study samples inside this radius
#' considered to be of European descent and samples outside this radius of
#' non-European descent. The radius is computed as the maximum Euclidean distance
#' of reference Europeans to the centre of reference European samples.
#' @param refSamples [data.frame] Dataframe with sample identifiers
#' [refSamplesIID] corresponding to IIDs in prefixMergedDataset.eigenvec and
#' population identifier [refSamplesPop] corresponding to population IDs
#' [refColorsPop] in refColorsfile/refColors. Either refSamples or
#' refSamplesFile have to be specified.
#' @param refColors [data.frame, optional] Dataframe with population IDs in
#' column [refColorsPop] and corresponding color-code for PCA plot in column
#' [refColorsColor]. If not provided and is.null(refColorsFile) default colors
#' are used.
#' @param refSamplesFile [character] /path/to/File/with/reference samples. Needs
#' columns with sample identifiers [refSamplesIID] corresponding to IIDs in
#' prefixMergedDataset.eigenvec and population identifier [refSamplesPop]
#' corresponding to population IDs [refColorsPop] in refColorsfile/refColors.
#' @param refColorsFile [character, optional]
#' /path/to/File/with/Population/Colors cotaining population IDs in column
#' [refColorsPop] and corresponding color-code for PCA plot in column
#' [refColorsColor].If not provided and is.null(refColors) default colors for
#' are used.
#' @param refSamplesIID [character] Column name of reference sample IDs in
#' refSamples/refSamplesFile.
#' @param refSamplesPop [character] Column name of reference sample population
#' IDs in refSamples/refSamplesFile.
#' @param refColorsColor [character] Column name of population colors in
#' refColors/refColorsFile
#' @param refColorsPop [character] Column name of reference sample population
#' IDs in refColors/refColorsFile.
#' @param studyColor [character] Color to be used for study population.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sampleQC) via ggplot2::ggsave(p=p_sampleQC,
#' other_arguments) or pdf(outfile) print(p_sampleQC) dev.off().
#' If TRUE, i) depicts the X-chromosomal heterozygosity (SNPSEX) of the samples
#' split by their PEDSEX (if do.evaluate_check_sex is TRUE), ii) creates a
#' scatter plot
#' with samples' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis (if do.evaluate_check_het_and_miss is TRUE), iii) depicts
#' all pair-wise IBD-estimates as histogram (if do.evaluate_check_relatedness is
#' TRUE)
#' and iv) creates a scatter plot of PC1 versus PC2 color-coded for samples of
#' reference populations and study population (if do.check_ancestry is set to
#' TRUE).
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_list, a named [list] with 1.
#' sample_missingness containing a [vector] with sample IIDs failing the
#' missingness threshold imissTh, 2. highIBD containing a [vector] with sample
#' IIDs failing the relatedness threshold highIBDTh, 3. outlying_heterozygosity
#' containing a [vector] with sample IIDs failing the heterozygosity threshold
#' hetTh, 4. mismatched_sex containing a [vector] with the sample IIDs failing
#' the sexcheck based on SNPSEX and femaleTh/maleTh and 5. ancestry containing
#' a vector with sample IIDs failing the ancestry check based on europeanTh and
#' ii) p_sampleQC, a ggplot2-object 'containing' a sub-panelled plot with the
#' QC-plots of \code{\link{check_sex}},
#' \code{\link{check_het_and_miss}},
#' \code{\link{check_relatedness}} and \code{\link{check_ancestry}}, which can
#' be shown by print(p_sampleQC).
#' List entries contain NULL if that specific check was not chosen.
#' @details perIndividualQC wraps around the individual QC functions
#' \code{\link{check_sex}}, \code{\link{check_het_and_miss}},
#' \code{\link{check_relatedness}} and \code{\link{check_ancestry}}. For details
#' on the parameters and outputs, check these function documentations. For
#' detailed output for fail IIDs (instead of simple IID lists), run each
#' function individually.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir, 'extdata')
#' qcdir <- tempdir()
#' name <- "data"
#' # All quality control checks
#' # In this examples, run_check* already conducted and outcome files present
#' # in qcdir, hence dont.check_* all set to FALSE
#' \dontrun{
#' fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
#' do.run_check_sex=FALSE, do.run_check_ancestry=FALSE)
#'
#' # Only check sex and missingness/heterozygosity
#' fail_sex_het_miss <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' dont.check_ancestry=TRUE, dont.check_relatedness=TRUE,
#' interactive=FALSE, verbose=FALSE, do.run_check_het_and_miss=FALSE,
#' do.run_check_sex=FALSE)
#' }
perIndividualQC <- function(indir, name, qcdir=indir,
                            dont.check_sex=FALSE,
                            do.run_check_sex=TRUE, do.evaluate_check_sex=TRUE,
                            maleTh=0.8, femaleTh=0.2,
                            externalSex=NULL, externalMale="M",
                            externalSexSex="Sex", externalSexID="IID",
                            externalFemale="F", fixMixup=FALSE,
                            dont.check_het_and_miss=FALSE,
                            do.run_check_het_and_miss=TRUE,
                            do.evaluate_check_het_and_miss=TRUE,
                            imissTh=0.03, hetTh=3,
                            dont.check_relatedness=FALSE,
                            do.run_check_relatedness=TRUE,
                            do.evaluate_check_relatedness=TRUE,
                            highIBDTh=0.1875,
                            dont.check_ancestry=FALSE,
                            do.run_check_ancestry=TRUE,
                            do.evaluate_check_ancestry=TRUE,
                            prefixMergedDataset, europeanTh=1.5,
                            refSamples=NULL, refColors=NULL,
                            refSamplesFile=NULL, refColorsFile=NULL,
                            refSamplesIID="IID", refSamplesPop="Pop",
                            refColorsColor="Color", refColorsPop="Pop",
                            studyColor="#2c7bb6",
                            interactive=FALSE, verbose=TRUE,
                            path2plink=NULL, showPlinkOutput=TRUE) {
    fail_sex <- NULL
    fail_het_imiss <- NULL
    fail_relatedness <- NULL
    fail_ancestry <- NULL

    p_sexcheck <- NULL
    p_het_imiss <- NULL
    p_relatedness <- NULL
    p_ancestry <- NULL

    prefix <- paste(qcdir, "/", name, sep="")
    if (!dont.check_sex) {
        if (do.run_check_sex) {
            run <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
                                 path2plink=path2plink,
                                 showPlinkOutput=showPlinkOutput,
                                 verbose=verbose)
        }
        if (do.evaluate_check_sex) {
            if (verbose) {
                message("Identification of individuals with discordant sex ",
                        "information")
            }
            fail_sex <- evaluate_check_sex(qcdir=qcdir, indir=indir, name=name,
                                           maleTh=maleTh, femaleTh=femaleTh,
                                           externalSex=externalSex,
                                           externalMale=externalMale,
                                           externalFemale=externalFemale,
                                           externalSexSex=externalSexSex,
                                           externalSexID=externalSexID,
                                           verbose=verbose,
                                           path2plink=path2plink,
                                           showPlinkOutput=showPlinkOutput,
                                           fixMixup=fixMixup, interactive=FALSE)
            if (!is.null(fail_sex$fail_sex)) {
                write.table(fail_sex$fail_sex[,1:2],
                            file=paste(prefix, ".fail-sexcheck.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            if (!is.null(fail_sex$mixup)) {
                write.table(fail_sex$mixup[,1:2],
                            file=paste(prefix, ".sexcheck_mixup.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            p_sexcheck <- fail_sex$p_sexcheck
        }
    }
    if (!dont.check_het_and_miss) {
        if (do.run_check_het_and_miss) {
            run_miss <- run_check_missingness(qcdir=qcdir, indir=indir, name=name,
                                              path2plink=path2plink,
                                              showPlinkOutput=showPlinkOutput,
                                              verbose=verbose)
            run_het <- run_check_heterozygosity(qcdir=qcdir, indir=indir,
                                                name=name,
                                                path2plink=path2plink,
                                                showPlinkOutput=showPlinkOutput,
                                                verbose=verbose)
        }
        if (do.evaluate_check_het_and_miss) {
            if (verbose) {
                message("Identification of individuals with outlying missing ",
                        "genotype or heterozygosity rates")
            }
            fail_het_imiss <-
                evaluate_check_het_and_miss(qcdir=qcdir, name=name,
                                            imissTh=imissTh,
                                            hetTh=hetTh,
                                            interactive=FALSE)
            if (!is.null(fail_het_imiss$fail_imiss)) {
                write.table(fail_het_imiss$fail_imiss[,1:2],
                            file=paste(prefix, ".fail-imiss.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            if (!is.null(fail_het_imiss$fail_het)) {
                write.table(fail_het_imiss$fail_het[,1:2],
                            file=paste(prefix, ".fail-het.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            p_het_imiss <- fail_het_imiss$p_het_imiss
        }
    }
    if (!dont.check_relatedness) {
        if (do.run_check_relatedness) {
            run <- run_check_relatedness(qcdir=qcdir, indir=indir, name=name,
                                          path2plink=path2plink,
                                          showPlinkOutput=showPlinkOutput,
                                          verbose=verbose)
        }
        if (do.evaluate_check_relatedness) {
            if (verbose) message("Identification of related individuals")
            fail_relatedness <- evaluate_check_relatedness(qcdir=qcdir, name=name,
                                                  imissTh=imissTh,
                                                  highIBDTh=highIBDTh,
                                                  interactive=FALSE)
            if (!is.null(fail_relatedness$failIDs)) {
               write.table(fail_relatedness$failIDs,
                           file=paste(prefix, ".fail-IBD.IDs", sep=""),
                           row.names=FALSE, quote=FALSE, col.names=TRUE,
                           sep="\t")
            }
            p_relatedness <- fail_relatedness$p_IBD
        }
    }
    if (!dont.check_ancestry) {
        if (do.run_check_ancestry) {
            run <- run_check_ancestry(qcdir=qcdir, indir=indir,
                                      prefixMergedDataset=prefixMergedDataset,
                                      path2plink=path2plink,
                                      showPlinkOutput=showPlinkOutput,
                                      verbose=verbose)
        }
        if (do.evaluate_check_ancestry) {
            if (verbose) {
                message("Identification of individuals of divergent ancestry")
            }
            fail_ancestry <- evaluate_check_ancestry(qcdir=qcdir, indir=indir,
                                                     name=name,
                                                     prefixMergedDataset=
                                                        prefixMergedDataset,
                                                     europeanTh=europeanTh,
                                                     refSamples=refSamples,
                                                     refColors=refColors,
                                                     refSamplesFile=
                                                         refSamplesFile,
                                                     refColorsFile=
                                                         refColorsFile,
                                                     refSamplesIID=
                                                         refSamplesIID,
                                                     refSamplesPop=
                                                         refSamplesPop,
                                                     refColorsColor=
                                                         refColorsColor,
                                                     refColorsPop=
                                                         refColorsPop,
                                                     studyColor=studyColor,
                                                     interactive=FALSE)
            if (!is.null(fail_ancestry$fail_ancestry)) {
                write.table(fail_ancestry$fail_ancestry,
                            file=paste(prefix, ".fail-ancestry.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            p_ancestry <- fail_ancestry$p_ancestry
        }
    }

    fail_list <- list(missing_genotype=as.vector(fail_het_imiss$fail_imiss$IID),
                      highIBD=as.vector(fail_relatedness$failIDs$IID),
                      outlying_heterozygosity=
                          as.vector(fail_het_imiss$fail_het$IID),
                      mismatched_sex=as.vector(fail_sex$fail_sex$IID),
                      ancestry=as.vector(fail_ancestry$fail_ancestry$IID))

    if(verbose) message(paste("Combine fail IDs into ", prefix,
                               ".fail.IDs", sep=""))
    uniqueFails <- unique(unlist(fail_list))
    write.table(cbind(uniqueFails, uniqueFails),
                        file=paste(prefix, ".fail.IDs",sep=""),
                        quote=FALSE, row.names=FALSE, col.names=FALSE)

    plots_sampleQC <- list(p_sexcheck, p_het_imiss, p_relatedness, p_ancestry)
    plots_sampleQC <- plots_sampleQC[sapply(plots_sampleQC,
                                            function(x) !is.null(x))]
    subplotLabels <- LETTERS[1:length(plots_sampleQC)]
    p_sampleQC <- cowplot::plot_grid(plotlist=plots_sampleQC,
                                     nrow=length(plots_sampleQC),
                                     labels=subplotLabels,
                                     rel_heights=c(rep(1,
                                                       length(plots_sampleQC)),
                                                       1.5))
    if (interactive) {
        print(p_sampleQC)
    }
    return(list(fail_list=fail_list, p_sampleQC=p_sampleQC))
}

#' Overview of per sample QC
#'
#' \code{overviewPerIndividualQC} depicts results of
#' \code{\link{perIndividualQC}} as intersection plots
#' (via \code{\link[UpSetR]{upset}}) and returns dataframes indicating
#' which QC checks individuals failed or passed.
#'
#'
#' @param results_perIndividualQC [list] Output of \code{\link{perIndividualQC}}
#' i.e. named [list]
#' with i) sample_missingness containing a [vector] with sample IIDs failing
#' the selected missingness threshold imissTh, ii) highIBD containing
#' a [vector] with sample IIDs failing the selected relatedness threshold
#' highIBDTh, iii) outlying_heterozygosity containing a [vector] with sample
#' IIDs failing selected the heterozygosity threshold hetTh, iv) mismatched_sex
#' containing a [vector] with the sample IIDs failing the sexcheck based on
#' SNPSEX and selected femaleTh/maleTh, v) ancestry containing a vector with
#' sample IIDs failing the ancestry check based on the selected europeanTh and
#' vi) p_sampleQC, a ggplot2-object 'containing' a sub-paneled plot with the
#' QC-plots of \code{\link{check_sex}},
#' \code{\link{check_het_and_miss}},
#' \code{\link{check_relatedness}} and \code{\link{check_ancestry}}.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_overview) via ggplot2::ggsave(p=p_overview,
#' other_arguments) or pdf(outfile) print(p_overview) dev.off().
#' @return Named [list] with i) nr_fail_samples: total number of samples
#' [integer] failing perIndividualQC, ii) fail_QC containing a [data.frame] with
#' samples that failed QC steps (excluding ancestry): samples IIDs in rows,
#' columns are all QC steps applied by perIndividualQC (max=4), with entries=0 if
#' passing the QC and entries=1 if failing that particular QC and iii)
#' fail_QC_and_ancestry containing a [data.frame] with samples that failed
#' ancestry and QC checks: samples IIDs in rows, columns are QC_fail and
#' Ancestry_fail, with entries=0 if passing and entries=1 if failing that check.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir, 'extdata')
#' qcdir <- tempdir()
#' name <- "data"
#' \dontrun{
#' fail_individuals <- perIndividualQC(qcdir=qcdir, indir=indir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
#' do.run_check_sex=FALSE, do.run_check_ancestry=FALSE)
#'
#' overview <- overviewPerIndividualQC(fail_individuals)
#' }

overviewPerIndividualQC <- function(results_perIndividualQC, interactive=FALSE) {
    list2counts <- function(element, all_names) {
        all_names[!(all_names %in% element)] <- 0
        all_names[all_names %in% element] <- 1
        return(as.numeric(all_names))
    }
    if (length(perIndividualQC) == 2 &&
        !all(names(results_perIndividualQC) == c("fail_list", "p_sampleQC"))) {
        stop("results_perIndividualQC not direct output of perIndividualQC")
    }
    fail_list <- results_perIndividualQC$fail_list
    # Remove null elements
    fail_list <- fail_list[!sapply(fail_list, is.null)]

    unique_samples_fail_all <- unique(unlist(fail_list))

    # a) overview QC fails independent of ethnicity
    fail_list_wo_ancestry <- fail_list[!names(fail_list) == "ancestry"]
    unique_samples_fail_wo_ancestry <- unique(unlist(fail_list_wo_ancestry))

    fail_counts_wo_ancestry <- sapply(fail_list_wo_ancestry, list2counts,
                                      unique_samples_fail_wo_ancestry)
    rownames(fail_counts_wo_ancestry) <- unique_samples_fail_wo_ancestry

    if (interactive) {
        if (length(fail_list_wo_ancestry) >= 2) {
            UpSetR::upset(UpSetR::fromList(fail_list_wo_ancestry),
                  order.by = "freq",
                  empty.intersections = "on", text.scale=1.2,
                  # Include when UpSetR v1.4.1 is released
                  # title="Overview quality control failures",
                  mainbar.y.label="Samples failing multiple QC checks",
                  sets.x.label="Sample fails per QC check",
                  main.bar.color="#1b9e77", matrix.color="#1b9e77",
                  sets.bar.color="#d95f02")
        } else {
            message("overviewSampleQC for QC fails cannot be displayed with ",
                    "UpSetR: at least two elements in list required, but only ",
                    length(fail_list_wo_ancestry) ," provided")
        }
    }

    if ("ancestry" %in% names(fail_list)) {
        # b) overview of QC and ancestry fails
        fail_all <- list(QC_fail=unique_samples_fail_wo_ancestry,
                    Ancestry_fail=fail_list$ancestry)
        fail_counts_all <- sapply(fail_all, list2counts,
                                  unique_samples_fail_all)
        rownames(fail_counts_all) <- unique_samples_fail_all
        if (interactive) {
            if (length(fail_all) >= 2) {
                UpSetR::upset(UpSetR::fromList(fail_all),
                              order.by = "freq",
                            # Include when UpSetR v1.4.1 is released
                      # title="Intersection between QC and ancestry failures",
                      mainbar.y.label="Samples failing QC and ancestry checks",
                      sets.x.label="Sample fails per QC check",
                      empty.intersections = "on", text.scale=1.2,
                      main.bar.color="#7570b3", matrix.color="#7570b3",
                    sets.bar.color="#e7298a" )
            } else {
                message("overviewSampleQC for QC fails and ancestry cannot be ",
                        "displayed with UpSetR: as no samples are present in ",
                        "QC fails")
            }
        }
        #p_overview <- cowplot::plot_grid(p_qc, p_qc_ancestry, nrow=2)
    } else {
        #p_overview <- p_qc
        fail_counts_all <- NULL
    }
    nr_fail_samples <- length(unique_samples_fail_all)
    return(list(nr_fail_samples=nr_fail_samples,
                fail_QC=fail_counts_wo_ancestry,
                fail_QC_and_ancestry=fail_counts_all))
                #p_overview=p_overview))
}

#' Identification of individuals with discordant sex information
#'
#' Runs and evaluates results from plink --check-sex.
#' \code{\link{check_sex}} returns IIDs for individuals whose SNPSEX != PEDSEX
#' (where the SNPSEX is determined by the heterozygosity rate across
#' X-chromosomal variants).
#' Mismatching SNPSEX and PEDSEX IDs can indicate plating errors, sample-mixup
#' or generally samples with poor genotyping. In the latter case, these IDs are
#' likely to fail other QC steps as well.
#' Optionally, an extra data.frame (externalSex) with sample IDs and sex can be
#' provided to double check if external and PEDSEX data (often processed at
#' different centers) match. If a mismatch between PEDSEX and SNPSEX was
#' detected, while SNPSEX == Sex, PEDSEX of these individuals can optionally be
#' updated (fixMixup=TRUE).
#' \code{\link{check_sex}} depicts the X-chromosomal heterozygosity (SNPSEX) of
#' the individuals split by their (PEDSEX).
#'
#' \code{\link{check_sex}} wraps around \code{\link{run_check_sex}}  and
#' \code{\link{evaluate_check_sex}} . If run.check_sex is TRUE,
#' \code{\link{run_check_sex} } is excuted; otherwise it is assumed that plink
#' --check-sex has been run externally and qcdir/name.sexcheck exists.
#' \code{\link{check_sex}}  will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_sex, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#sexcheck}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.sexcheck as returned
#' by plink --check-sex. Per default qcdir=indir. If run.check_sex is FALSE, it
#' is assumed that plink --check-sex has been run and qcdir/name.sexcheck is
#' present. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam and name.sexcheck.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param run.check_sex [logical] Should plink --check-sex be run? if set to
#' FALSE, it is assumed that plink --check-sex has been run and
#' qcdir/name.sexcheck is present; \code{\link{check_sex}}  will fail with
#' missing file error otherwise.
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex while Sex==SNPSEX automatically corrected: this will directly
#' change the name.bim/.bed/.fam files!
#' @param externalSex [data.frame, optional] Dataframe with sample IDs
#' [externalSexID] and sex [externalSexSex] to double check if external and
#' PEDSEX data (often processed at different centers) match.
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @param externalSexSex [character] Column identifier for column containing sex
#' information in externalSex.
#' @param externalSexID [character] Column identifier for column containing ID
#' information in externalSex.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sexcheck) via ggplot2::ggsave(p=p_sexcheck,
#' other_arguments) or pdf(outfile) print(p_sexcheck) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named list with i) fail_sex: [data.frame] with FID, IID, PEDSEX,
#' SNPSEX and Sex (if externalSex was provided) of individuals failing sex
#' check, ii) mixup: dataframe with FID, IID, PEDSEX, SNPSEX and Sex (if
#' externalSex was provided) of individuals whose PEDSEX != Sex and Sex ==
#' SNPSEX and iii) p_sexcheck, a ggplot2-object 'containing' a scatter plot of
#' the X-chromosomal heterozygosity (SNPSEX) of the sample split by their
#' (PEDSEX), which can be shown by print(p_sexcheck).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_sex <- check_sex(indir=qcdir, name=name, run.check_sex=FALSE,
#' interactive=FALSE, verbose=FALSE)
check_sex <- function(indir, name, qcdir=indir, maleTh=0.8, femaleTh=0.2,
                      run.check_sex=TRUE,
                      externalSex=NULL,
                      externalFemale="F", externalMale="M",
                      externalSexSex="Sex", externalSexID="IID",
                      fixMixup=FALSE,
                      interactive=FALSE, verbose=FALSE,
                      path2plink=NULL, showPlinkOutput=TRUE) {
    if (run.check_sex) {
        run_sex <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
                                 verbose=verbose,
                                 path2plink=path2plink,
                                 showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_sex(qcdir=qcdir, name=name, externalSex=externalSex,
                               maleTh=maleTh, femaleTh=femaleTh,
                               interactive=interactive, fixMixup=fixMixup,
                               indir=indir,
                               externalFemale=externalFemale,
                               externalMale=externalMale,
                               externalSexSex=externalSexSex,
                               externalSexID=externalSexID,
                               verbose=verbose, path2plink=path2plink,
                               showPlinkOutput=showPlinkOutput)
    return(fail)
}


#' Identification of individuals with outlying missing genotype or
#' heterozygosity rates
#'
#' Runs and evaluates results from plink --missing (missing genotype rates
#' per individual) and plink --het (heterozygosity rates per individual).
#' Non-systematic failures in genotyping and outlying heterozygosity (hz) rates
#' per individual are often proxies for DNA sample quality. Larger than expected
#' heterozygosity can indicate possible DNA contamination.
#' The mean heterozygosity in PLINK is computed as hz_mean = (N-O)/N, where
#' N: number of non-missing genotypes and O:observed number of homozygous
#' genotypes for a given individual.
#' Mean heterozygosity can differ between populations and SNP genotyping panels.
#' Within a population and genotyping panel, a reduced heterozygosity rate can
#' indicate inbreeding - these individuals will then likely be returned by
#' \code{\link{check_relatedness}} as individuals that fail the relatedness
#' filters. \code{check_het_and_miss} creates a scatter plot with the
#' individuals' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis.
#'
#' \code{\link{check_het_and_miss}} wraps around
#' \code{\link{run_check_missingness}},
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{evaluate_check_het_and_miss}}.
#' If run.check_het_and_miss is TRUE, \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}} are excuted; otherwise it is assumed
#' that plink --missing and plink --het have been run externally and
#' qcdir/name.het and qcdir/name.imiss exist.  \code{\link{check_het_and_miss}}
#' will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_imiss and fail_het, check the
#' original description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#imiss} and
#' \url{https://www.cog-genomics.org/plink/1.9/formats#het}
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where name.het as returned by
#' plink --het and name.imiss as returned by plink --missing will be saved. Per
#' default qcdir=indir. If run.check_het_and_miss is FALSE, it is assumed that
#' plink --missing and plink --het have been run and qcdir/name.imiss and
#' qcdir/name.het are present. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.het and name.imiss.
#' @param run.check_het_and_miss [logical] Should plink --missing and plink
#' --het be run to determine genotype missingness and heterozygosity rates; if
#' FALSE, it is assumed that plink --missing and plink --het have been run and
#' qcdir/name.imiss and qcdir/name.het are present;
#' \code{\link{check_het_and_miss}} will fail with missing file error otherwise.
#' @param imissTh [double] Threshold for acceptable missing genotype rate per
#' individual; has to be proportion between (0,1)
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity per individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check; has to be
#' larger than 0.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_het_imiss) via ggplot2::ggsave(p=p_het_imiss
#' , other_arguments) or pdf(outfile) print(p_het_imiss) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_imiss [data.frame] containing FID (Family
#' ID), IID (Within-family ID), MISS_PHENO (Phenotype missing? (Y/N)), N_MISS
#' (Number of missing genotype call(s), not including obligatory missings),
#' N_GENO (Number of potentially valid call(s)), F_MISS (Missing call rate) of
#' individuals failing missing genotype check and ii) fail_het [data.frame]
#' containing FID (Family ID), IID (Within-family ID), O(HOM) (Observed number
#' of homozygotes), E(HOM) (Expected number of homozygotes), N(NM) (Number of
#' non-missing autosomal genotypes), F (Method-of-moments F coefficient
#' estimate) of individuals failing  outlying heterozygosity check and iii)
#' p_het_imiss, a ggplot2-object 'containing' a scatter plot with the samples'
#' missingness rates on x-axis and their heterozygosity rates on the y-axis,
#' which can be shown by print(p_het_imiss).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_het_miss <- check_het_and_miss(indir=qcdir, name=name,
#' run.check_het_and_miss=FALSE, interactive=FALSE)
check_het_and_miss <- function(indir, name, qcdir=indir, imissTh=0.03, hetTh=3,
                               run.check_het_and_miss=TRUE,
                               interactive=FALSE, verbose=FALSE,
                               path2plink=NULL, showPlinkOutput=TRUE) {
    if (run.check_het_and_miss) {
        run_het <- run_check_heterozygosity(indir=indir,qcdir=qcdir, name=name,
                                            verbose=verbose,
                                            path2plink=path2plink,
                                            showPlinkOutput=showPlinkOutput)
        run_miss <- run_check_missingness(indir=indir, qcdir=qcdir, name=name,
                                          verbose=verbose,
                                          path2plink=path2plink,
                                          showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,  hetTh=hetTh,
                                 imissTh=imissTh, interactive=interactive)
    return(fail)
}

#' Identification of related individuals
#'
#' Runs and evaluates results from plink --genome.
#' plink --genome calculates identity by state (IBS) for each pair of
#' individuals based on the average proportion of alleles shared at genotyped
#' SNPs. The degree of recent shared ancestry, i.e. the identity by descent
#' (IBD) can be estimated from the genome-wide IBS. The proportion of IBD
#' between two individuals is returned by plink --genome as PI_HAT.
#' check_relatedness finds pairs of samples whose proportion of IBD is larger
#' than the specified highIBDTh. Subsequently, for pairs of individuals that do
#' not have additional relatives in the dataset, the individual with the greater
#' genotype missingness rate is selected and returned as the individual failing
#' the relatedness check. For more complex family structures, the unrelated
#' individuals per family are selected (e.g. in a parents-offspring trio, the
#' offspring will be marked as fail, while the parents will be kept in the
#' analysis).
#' \code{check_relatedness} depicts all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT.
#'
#' \code{\link{check_relatedness}} wraps around
#' \code{\link{run_check_relatedness}} and
#' \code{\link{evaluate_check_relatedness}}. If run.check_relatedness is TRUE,
#' \code{\link{run_check_relatedness}} is excuted; otherwise it is assumed that
#' plink --genome has been run externally and qcdir/name.genome exists.
#' \code{\link{check_relatedness}}  will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_high_IBD, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#genome}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to where name.genome as returned
#' by plink --genome will be saved.  Per default qcdir=indir. If
#' run.check_relatedness is FALSE, it is assumed that plink
#' --missing and plink --genome have been run and qcdir/name.imiss and
#' qcdir/name.genome exist. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.genome and name.imiss.
#' @param run.check_relatedness [logical] Should plink --genome be run to
#' determine pairwise IBD of individuals; if FALSE, it is assumed that
#' plink --genome and plink --missing have been run and qcdir/name.imiss and
#' qcdir/name.genome are present;
#' \code{\link{check_relatedness}} will fail with missing file error otherwise.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_IBD() via ggplot2::ggsave(p=p_IBD,
#' other_arguments) or pdf(outfile) print(p_IBD) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_high_IBD containing a [data.frame] of
#' IIDs and FIDs of individuals who fail the IBDTh in columns
#' FID1 and IID1. In addition, the following columns are returned (as originally
#' obtained by plink --genome):
#' FID2 (Family ID for second sample), IID2 (Individual ID for second sample),
#' RT (Relationship type inferred from .fam/.ped file), EZ (IBD sharing expected
#' value, based on just .fam/.ped relationship), Z0 (P(IBD=0)), Z1 (P(IBD=1)),
#' Z2 (P(IBD=2)), PI_HAT (Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)), PHE
#' (Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)),
#' DST (IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)), PPC (IBS
#' binomial test), RATIO (HETHET : IBS0 SNP ratio (expected value 2)).
#' and ii) failIDs containing a [data.frame] with individual IDs [IID] and
#' family IDs [FID] of individuals failing the highIBDTh iii) p_IBD, a
#' ggplot2-object 'containing' all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT, which can be
#' shown by print(p_IBD).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' relatednessQC <- check_relatedness(indir=qcdir, name=name, interactive=FALSE,
#' run.check_relatedness=FALSE)
check_relatedness <- function(indir, name, qcdir=indir, highIBDTh=0.1875,
                              imissTh=0.03, run.check_relatedness=TRUE,
                              interactive=FALSE, verbose=FALSE,
                              path2plink=NULL, showPlinkOutput=TRUE) {
    if (run.check_relatedness) {
        run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name,
                                     verbose=verbose,
                                     path2plink=path2plink,
                                     showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_relatedness(qcdir=qcdir, name=name,
                                       highIBDTh=highIBDTh,
                                       imissTh=imissTh, interactive=interactive,
                                       verbose=verbose)
    return(fail)
}



#' Identification of individuals of divergent ancestry
#'
#' Runs and evaluates results of plink --pca on merged genotypes from
#' individuals to be QCed and individuals of reference population of known
#' genotypes. Currently, check ancestry only supports automatic selection of
#' individuals of European descent. It uses information from principal
#' components 1 and 2 returned by plink --pca to find the center of the European
#' reference samples (mean(PC1_europeanRef), mean(PC2_europeanRef). It then
#' computes the maximum Euclidean distance (maxDist) of the European reference
#' samples from this centre. All study samples whose Euclidean distance from the
#' centre falls outside the circle described by the radius r=europeanTh* maxDist
#' are considered non-European and their IDs are returned as failing the
#' ancestry check.
#' \code{check_ancestry} creates a scatter plot of PC1 versus PC2 colour-coded
#' for samples of the reference populations and the study population.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where
#' prefixMergedDataset.eigenvec results as returned by plink --pca should be
#' saved. Per default qcdir=indir. If run.check_ancestry is FALSE, it is assumed
#' that plink --pca prefixMergedDataset has been run and
#' qcdir/prefixMergedDataset.eigenvec exists.User needs writing permission to
#' qcdir.
#' @param name [character] prefix of plink files, i.e. name.bed, name.bim,
#' name.fam.
#' @param prefixMergedDataset [character] Prefix of merged dataset (study and
#' reference samples) used in plink --pca, resulting in
#' prefixMergedDataset.eigenvec.
#' @param europeanTh [double] Scaling factor of radius to be drawn around center
#' of European reference samples, with study samples inside this radius
#' considered to be of European descent and samples outside this radius of
#' non-European descent. The radius is computed as the maximum Euclidean
#' distance of European reference samples to the centre of European reference
#' samples.
#' @param refSamples [data.frame] Dataframe with sample identifiers
#' [refSamplesIID] corresponding to IIDs in prefixMergedDataset.eigenvec and
#' population identifier [refSamplesPop] corresponding to population IDs
#' [refColorsPop] in refColorsfile/refColors. Either refSamples or
#' refSamplesFile have to be specified.
#' @param refColors [data.frame, optional] Dataframe with population IDs in
#' column [refColorsPop] and corresponding colour-code for PCA plot in column
#' [refColorsColor]. If not provided and is.null(refColorsFile) default colors
#' are used.
#' @param refSamplesFile [character] /path/to/File/with/reference samples. Needs
#' columns with sample identifiers [refSamplesIID] corresponding to IIDs in
#' prefixMergedDataset.eigenvec and population identifier [refSamplesPop]
#' corresponding to population IDs [refColorsPop] in refColorsfile/refColors.
#' @param refColorsFile [character, optional]
#' /path/to/File/with/Population/Colors cotaining population IDs in column
#' [refColorsPop] and corresponding colour-code for PCA plot in column
#' [refColorsColor].If not provided and is.null(refColors) default colors for
#' are used.
#' @param refSamplesIID [character] Column name of reference sample IDs in
#' refSamples/refSamplesFile.
#' @param refSamplesPop [character] Column name of reference sample population
#' IDs in refSamples/refSamplesFile.
#' @param refColorsColor [character] Column name of population colors in
#' refColors/refColorsFile
#' @param refColorsPop [character] Column name of reference sample population
#' IDs in refColors/refColorsFile.
#' @param studyColor [character] Colour to be used for study population.
#' @param run.check_ancestry [logical] Should plink --pca be run to
#' determine principal components of merged dataset; if FALSE, it is assumed
#' that plink --pca has been run successfuly and
#' qcdir/prefixMergedDataset.eigenvec is present;
#' \code{\link{check_ancestry}} will fail with missing file error otherwise.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_ancestry) via ggplot2::ggsave(p=p_ancestry,
#' other_arguments) or pdf(outfile) print(p_ancestry) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_ancestry, containing a [data.frame] with
#' FID and IID of non-European individuals and ii) p_ancestry, a ggplot2-object
#' 'containing' a scatter plot of PC1 versus PC2 colour-coded for samples of the
#' reference populations and the study population, which can be shown by
#' print(p_ancestry).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_ancestry <- check_ancestry(indir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE,
#' run.check_ancestry=FALSE)

check_ancestry <- function(indir, name, qcdir=indir, prefixMergedDataset,
                           europeanTh=1.5,
                           refSamples=NULL, refColors=NULL,
                           refSamplesFile=NULL, refColorsFile=NULL,
                           refSamplesIID="IID", refSamplesPop="Pop",
                           refColorsColor="Color", refColorsPop="Pop",
                           studyColor="#2c7bb6",
                           run.check_ancestry=TRUE,
                           interactive=FALSE, verbose=verbose,
                           path2plink=NULL, showPlinkOutput=TRUE) {
        if (run.check_ancestry) {
            run <- run_check_ancestry(indir=indir, qcdir=qcdir,
                                      prefixMergedDataset=prefixMergedDataset,
                                      verbose=verbose,
                                      path2plink=path2plink,
                                      showPlinkOutput=showPlinkOutput)
        }
        fail <- evaluate_check_ancestry(qcdir=qcdir,indir=indir, name=name,
                                        prefixMergedDataset=prefixMergedDataset,
                                        europeanTh=europeanTh,
                                        refSamples=refSamples,
                                        refColors=refColors,
                                        refSamplesFile=refSamplesFile,
                                        refColorsFile=refColorsFile,
                                        refSamplesIID=refSamplesIID,
                                        refSamplesPop=refSamplesPop,
                                        refColorsColor=refColorsColor,
                                        refColorsPop=refColorsPop,
                                        studyColor=studyColor,
                                        interactive=interactive)
        return(fail)
}


#' Run PLINK sexcheck
#'
#' Run  plink --sexcheck to calculate the heterozygosity rate
#' across X-chromosomal variants.
#'
#' Both \code{\link{run_check_sex}} and its evaluation
#' \code{\link{evaluate_check_sex}} can simply be invoked by
#' \code{\link{check_sex}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.sexcheck as returned
#' by plink --check-sex. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' qcdir <- tempdir()
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run <- run_check_sex(indir=indir, qcdir=qcdir, name=name)
#' }
run_check_sex <- function(indir, name, qcdir=indir, verbose=FALSE,
                          path2plink=NULL,showPlinkOutput=TRUE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (verbose) message("Run check_sex via plink --check-sex")
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --check-sex",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
}

#' Evaluate results from PLINK sex check.
#'
#' Evaluates and depicts results from plink --check-sex (via
#' \code{\link{run_check_sex}} or externally conducted sex check).
#' Takes file qcdir/name.sexcheck and returns IIDs for samples whose
#' SNPSEX != PEDSEX (where the SNPSEX is determined by the heterozygosity rate
#' across X-chromosomal variants).
#' Mismatching SNPSEX and PEDSEX IDs can indicate plating errors, sample-mixup
#' or generally samples with poor genotyping. In the latter case, these IDs are
#' likely to fail other QC steps as well.
#' Optionally, an extra data.frame (externalSex) with sample IDs and sex can be
#' provided to double check if external and PEDSEX data (often processed at
#' different centers) match. If a mismatch between PEDSEX and SNPSEX was
#' detected while SNPSEX == Sex, PEDSEX of these individuals can optionally be
#' updated (fixMixup=TRUE).
#' \code{evaluate_check_sex} depicts the X-chromosomal heterozygosity (SNPSEX)
#' of the samples split by their (PEDSEX).
#'
#' Both \code{\link{run_check_sex}} and \code{\link{evaluate_check_sex}} can
#' simply be invoked by \code{\link{check_sex}}.
#'
#' For details on the output data.frame fail_sex, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#sexcheck}.
#'
#' @param qcdir [character] /path/to/directory containing name.sexcheck as
#' returned by plink --check-sex.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam and name.sexcheck.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex, with Sex==SNPSEX automatically corrected: this will directly
#' change the name.bim/.bed/.fam files!
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files; only required of fixMixup==TRUE. User
#' needs writing permission to indir.
#' @param externalSex [data.frame, optional] with sample IDs [externalSexID] and
#' sex [externalSexSex] to double check if external and PEDSEX data (often
#' processed at different centers) match.
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @param externalSexSex [character] Column identifier for column containing sex
#' information in externalSex.
#' @param externalSexID [character] Column identifier for column containing ID
#' information in externalSex.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sexcheck) via ggplot2::ggsave(p=p_sexcheck,
#' other_arguments) or pdf(outfile) print(p_sexcheck) dev.off().
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return named list with i) fail_sex: dataframe with FID, IID, PEDSEX, SNPSEX
#' and Sex (if externalSex was provided) of individuals failing sex check,
#' ii) mixup: dataframe with FID, IID, PEDSEX, SNPSEX and Sex (if externalSex
#' was provided) of individuals whose PEDSEX != Sex and Sex == SNPSEX and iii)
#' p_sexcheck, a ggplot2-object 'containing' a scatter plot of the X-chromosomal
#' heterozygosity (SNPSEX) of the individuals split by their (PEDSEX), which can
#' be shown by print(p_sexcheck).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_sex <- evaluate_check_sex(qcdir=qcdir, name=name, interactive=FALSE,
#' verbose=FALSE)
evaluate_check_sex <- function(qcdir, name, maleTh=0.8,
                               femaleTh=0.2, externalSex=NULL,
                               fixMixup=FALSE, indir=qcdir,
                               externalFemale="F", externalMale="M",
                               externalSexSex="Sex", externalSexID="IID",
                               verbose=FALSE,
                               path2plink=NULL,
                               showPlinkOutput=TRUE,
                               interactive=FALSE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(out, ".sexcheck", sep=""))){
        stop("plink --check-sex results file: ", out,
             ".sexcheck does not exist.")
    }
    testNumerics(numbers=c(maleTh, femaleTh), positives=c(maleTh, femaleTh),
                 proportions=c(maleTh, femaleTh))
    sexcheck <- read.table(paste(out, ".sexcheck",sep=""),
                           header=TRUE, stringsAsFactors=FALSE)

    names_sexcheck <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
    if (!all(names_sexcheck == names(sexcheck))) {
        stop("Header of", out, ".sexcheck is not correct. Was your
             file generated with plink --check-sex?")
    }
    if (is.null(externalSex)) {
        fail_sex <- sexcheck[sexcheck$STATUS == "PROBLEM",]
        if (nrow(fail_sex) == 0) fail_sex <- NULL
        mixup_geno_pheno <- NULL
    } else {
        if (!(externalSexSex %in% names(externalSex))) {
            stop("Column ", externalSexSex, " not found in externalSex!")
        }
        if (!(externalSexID %in% names(externalSex))) {
            stop("Column ", externalSexID, " not found in externalSex!")
        }
        names(externalSex)[names(externalSex) == externalSexSex] <- "Sex"
        names(externalSex)[names(externalSex) == externalSexID] <- "IID"

        sexcheck_fuse <- merge(sexcheck, externalSex, by="IID")
        sex_mismatch <-
            apply(dplyr::select_(sexcheck_fuse, ~Sex, ~PEDSEX, ~SNPSEX), 1,
                  function(ind) {
                      # Mismatch Sex in pheno and geno files
                      if (ind[1] == externalFemale && ind[2] %in% c(0, 1)) {
                          return(ifelse(ind[3] == 2, FALSE, TRUE))
                      }
                      if (ind[1] == externalMale && ind[2] %in% c(0, 2)) {
                          return(ifelse(ind[3] == 1, FALSE, TRUE))
                      }
                      # Mismatch Sex assigned and estimated
                      if (ind[1] == externalFemale && ind[2] == 2) {
                          return(ifelse(ind[3] == 1, TRUE, NA))
                      }
                      if (ind[1] == externalMale && ind[2] == 1) {
                          return(ifelse(ind[3] == 2, TRUE, NA))
                      }
                  })
        # SNPSEX != (Sex in pheno file, PEDSEX)
        fail_sex <-
            dplyr::select_(sexcheck_fuse, ~FID, ~IID, ~Sex, ~PEDSEX,
                           ~SNPSEX, ~F)[which(sex_mismatch),]
        if (nrow(fail_sex) == 0) {
            fail_sex <- NULL
            mixup_geno_pheno <- NULL
        } else {
            # SNPSEX == Sex in pheno file and PEDSEX != (SNPSEX, Sex)
            mixup_geno_pheno <-
                dplyr::select_(sexcheck_fuse, ~FID, ~IID, ~Sex, ~PEDSEX,
                               ~SNPSEX, ~F)[which(!sex_mismatch),]
            # Fix mismatch between PEDSEX and sex
            if (fixMixup) {
                if (!is.null(path2plink)) {
                    path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
                }
                checkPlink(path2plink)
                if (nrow(mixup_geno_pheno) != 0) {
                    file_mixup <- paste(out, ".mismatched_sex_geno_pheno",
                                        sep="")
                    write.table(dplyr::select_(mixup_geno_pheno, ~FID, ~IID,
                                               ~SNPSEX),
                                file=file_mixup,
                                row.names=FALSE, quote=FALSE, col.names=FALSE)
                    if (!file.exists(paste(prefix, ".fam", sep=""))){
                        stop("plink family file: ", prefix,
                             ".fam does not exist.")
                    }
                    if (!file.exists(paste(prefix, ".bim", sep=""))){
                        stop("plink snp file: ", prefix,
                             ".bim does not exist.")
                    }
                    if (!file.exists(paste(prefix, ".bed", sep=""))){
                        stop("plink binary file: ", prefix,
                             ".bed does not exist.")
                    }
                    system(paste(path2plink, "plink --bfile ", prefix,
                                 " --update-sex ", file_mixup,
                                 " --make-bed ",
                                 " --out ", prefix, sep=""),
                           ignore.stdout=showPlinkOutput,
                           ignore.stderr=!showPlinkOutput)
                } else {
                    if (verbose) {
                        message("All assigned genotype sexes (PEDSEX) match",
                                " external sex assignment (Sex)")
                    }
                    mixup_geno_pheno <- NULL
                }
            } else {
                # Append mismatch to fail_sex, ie treat as fail IDs
                if (nrow(mixup_geno_pheno) != 0) {
                    fail_sex <- rbind(fail_sex, mixup_geno_pheno)
                } else {
                    mixup_geno_pheno <- NULL
                }
            }
        }
    }
    sexcheck$PEDSEX <- as.factor(sexcheck$PEDSEX)
    p_sexcheck <- ggplot()
    p_sexcheck <- p_sexcheck + geom_point(data=sexcheck,
                                          aes_string(x='PEDSEX', y='F',
                                                     color='PEDSEX')) +
        ggtitle("Check assigned sex versus SNP sex") +
        xlab("Reported Sex (PEDSEX)") +
        ylab("ChrX heterozygosity") +
        ggrepel::geom_label_repel(data=data.frame(x=fail_sex$PEDSEX,
                                                  y=fail_sex$F,
                                                  label=fail_sex$IID),
                                  aes_string(x='x',y='y', label='label'),
                                  size=2) +
        geom_segment(data=data.frame(x=0.8, xend=1.2, y=maleTh,
                                     yend=maleTh),
                     aes_string(x='x', xend='xend', y='y', yend='yend'), lty=2,
                     color="#e7298a") +
        geom_segment(data=data.frame(x=1.8, xend=2.2, y=femaleTh,
                                     yend=femaleTh), lty=2,
                     aes_string(x='x', xend='xend', y='y', yend='yend'),
                     color="#e7298a") +
        theme_bw()
    if (length(unique(sexcheck$PEDSEX)) == 2) {
        p_sexcheck <- p_sexcheck +
            scale_color_manual(values=c("#377eb8", "#e41a1c"),
                               guide=FALSE) +
            scale_x_discrete(labels=c("Male", "Female"))
    } else if (length(unique(sexcheck$PEDSEX)) == 3 ) {
        p_sexcheck <- p_sexcheck +
            scale_color_manual(values=c("#999999", "#377eb8", "#e41a1c"),
                               guide=FALSE) +
            scale_x_discrete(labels=c("Unassigned", "Male", "Female"))
    }
    if (interactive) {
        print(p_sexcheck)
    }
    return(list(fail_sex=fail_sex, mixup=mixup_geno_pheno,
                p_sexcheck=p_sexcheck))
}

#' Run PLINK heterozygosity rate calculation
#'
#' Run plink --het to calculate heterozygosity rates per individual.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and their evaluation by
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.het as returned by
#' plink --het. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' qcdir <- tempdir()
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run <- run_check_heterozygosity(indir=indir, qcdir=qcdir, name=name)
#' }
run_check_heterozygosity <- function(indir, name, qcdir=indir, verbose=FALSE,
                                     path2plink=NULL,
                                     showPlinkOutput=TRUE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (verbose) message("Run check_heterozygosity via plink --het")
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --het",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
}

#' Run PLINK missingness rate calculation
#'
#' Run  plink --missing to calculate missing genotype rates
#' per individual.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and their evaluation by
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.imiss as returned
#' by plink --missing. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' qcdir <- tempdir()
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run <- run_check_missingnessness(indir=indir, qcdir=qcdir, name=name)
#' }
run_check_missingness <- function(indir, name, qcdir=indir, verbose=FALSE,
                                  path2plink=NULL,
                                  showPlinkOutput=TRUE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (verbose) message("Run check_missingness via plink --missing")
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --missing",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
}

#' Evaluate results from PLINK missing genotype and heterozygosity rate check.
#'
#' Evaluates and depicts results from plink --missing (missing genotype rates
#' per individual) and plink --het (heterozygosity rates per individual) via
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}} or externally conducted check.)
#' Non-systematic failures in genotyping and outlying heterozygosity (hz) rates
#' per individual are often proxies for DNA sample quality. Larger than expected
#' heterozygosity can indicate possible DNA contamination.
#' The mean heterozygosity in Plink is computed as hz_mean = (N-O)/N, where
#' N: number of non-missing genotypes and O:observed number of homozygous
#' genotypes for a given individual.
#' Mean heterozygosity can differ between populations and SNP genotyping panels.
#' Within a population and genotyping panel, a reduced heterozygosity rate can
#' indicate inbreeding - these individuals will then be returned by
#' check_relatedness as individuals that fail the relatedness filters.
#' \code{evaluate_check_het_and_miss} creates a scatter plot with the
#' individuals' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' For details on the output data.frame fail_imiss and fail_het, check the
#' original description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#imiss} and
#' \url{https://www.cog-genomics.org/plink/1.9/formats#het}
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' name.imiss and name.het results as returned by plink --missing and plink
#' --het.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.het and name.imiss.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity in any individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check; has to be
#' larger than 0.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_het_imiss) via ggplot2::ggsave(p=p_het_imiss
#' , other_arguments) or pdf(outfile) print(p_het_imiss) dev.off().
#' @return named [list] with i) fail_imiss dataframe containing FID (Family ID),
#' IID (Within-family ID), MISS_PHENO (Phenotype missing? (Y/N)), N_MISS (Number
#' of missing genotype call(s), not including obligatory missings), N_GENO (
#' Number of potentially valid call(s)), F_MISS (Missing call rate) of
#' individuals failing missing genotype check and ii) fail_het dataframe
#' containing FID (Family ID), IID (Within-family ID), O(HOM) (Observed number
#' of homozygotes), E(HOM) (Expected number of homozygotes), N(NM) (Number of
#' non-missing autosomal genotypes), F (Method-of-moments F coefficient
#' estimate) of individuals failing  outlying heterozygosity check and iii)
#' p_het_imiss, a ggplot2-object 'containing' a scatter plot with the samples'
#' missingness rates on x-axis and their heterozygosity rates on the y-axis,
#' which can be shown by print(p_het_imiss).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_het_miss <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,
#' interactive=FALSE)
evaluate_check_het_and_miss <- function(qcdir, name, imissTh=0.03,
                                  hetTh=3, interactive=FALSE) {
    prefix <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".imiss",sep=""))){
        stop("plink --missing output file: ", prefix,
             ".imiss does not exist.")
    }
    if (!file.exists(paste(prefix, ".het",sep=""))){
        stop("plink --het output file: ", prefix,
             ".het does not exist.")
    }
    testNumerics(numbers=c(imissTh, hetTh), positives=c(imissTh, hetTh),
                 proportions=imissTh)
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_imiss == names(imiss))) {
        stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    fail_imiss <- imiss[imiss$F_MISS > imissTh,]

    names_het <- c("FID", "IID", "O.HOM.", "E.HOM.", "N.NM.", "F")
    het <- read.table(paste(prefix, ".het", sep=""), header=TRUE)
    if (!all(names_het == names(het))) {
        stop("Header of ", prefix, ".het is not correct. Was your
             file generated with plink --het?")
    }
    fail_het <- het[het$F < (mean(het$F)  - hetTh*sd(het$F)) |
                        het$F > (mean(het$F) + hetTh*sd(het$F)),]

    nr_samples <- nrow(imiss)
    imiss$logF_MISS <- log10(imiss$F_MISS)
    het_imiss <- merge(imiss, het, by="IID")
    fail_het_imiss <- het_imiss[which(het_imiss$IID %in%
                                          union(fail_het$IID, fail_imiss$IID)),]
    het_imiss$type <- 1
    het_imiss$type[het_imiss$IID %in% fail_het$IID] <- 2
    het_imiss$type[het_imiss$IID %in% fail_imiss$IID] <- 3
    het_imiss$type[het_imiss$IID %in%
                       intersect(fail_het$IID, fail_imiss$IID)] <- 4
    het_imiss$type <- as.factor(het_imiss$type)
    minus_sd <- mean(het_imiss$F) - 1:5*(sd(het_imiss$F))
    plus_sd <- mean(het_imiss$F) + 1:5*(sd(het_imiss$F))
    p_het_imiss <- ggplot()
    p_het_imiss <- p_het_imiss + geom_point(data=het_imiss,
                                            aes_string(x='logF_MISS', y='F',
                                                       color='type')) +
        scale_color_manual(values=c("#666666", "#1b9e77", "#d95f02",
                                    "#7570b3"), guide=FALSE) +
        xlab("Proportion of missing SNPs") +
        ylab("Heterozygosity rate (and sd)") +
        ggtitle("Heterozygosity by Missingness across samples") +
        geom_hline(yintercept=c(minus_sd[1:3], plus_sd[1:3]), lty=2,
                   col="azure4") +
        scale_y_continuous(labels=c("-5", "-4", "-3" ,"+3", "+4", "+5"),
                           breaks=c(minus_sd[3:5], plus_sd[3:5])) +
        scale_x_continuous(labels=c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1),
                           breaks=c(-4,-3,-2, log10(0.03), log10(0.05),-1,0)) +
        geom_hline(yintercept=mean(het_imiss$F) - (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_hline(yintercept=mean(het_imiss$F) + (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_vline(xintercept=log10(imissTh), col="#e7298a", lty=2) +
        ggrepel::geom_label_repel(data=data.frame(x=fail_het_imiss$logF_MISS,
                                                  y=fail_het_imiss$F,
                                                  label=fail_het_imiss$IID),
                                  aes_string(x='x',y='y', label='label'),
                                  size=2) +
        theme_bw()
    if (interactive) {
        print(p_het_imiss)
    }
    return(list(fail_imiss=fail_imiss, fail_het=fail_het,
                p_het_imiss=p_het_imiss))
}

#' Run PLINK IBD estimation
#'
#' Run LD pruning on dataset with plink --exclude range highldfile
#' --indep-pairwise 50 5 0.2, where highldfile contains regions of high LD as
#' provided by Anderson et (2010) Nature Protocols. Subsequently, plink
#' --genome is run on the LD pruned data. plink --genome calculates identity
#' by state (IBS) for each pair of individuals based on the average proportion
#' of alleles shared at genotyped SNPs. The degree of recent shared ancestry,
#' i.e. the identity by descent (IBD) can be estimated from the genome-wide IBS.
#' The proportion of IBD between two individuals is returned by --genome as
#' PI_HAT.
#'
#' Both \code{\link{run_check_relatedness}} and its evaluation via
#' \code{\link{evaluate_check_relatedness}} can simply be invoked by
#' \code{\link{check_relatedness}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.genome as returned
#' by plink --genome. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' qcdir <- tempdir()
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name)
#' }
run_check_relatedness <- function(indir, name, qcdir=indir, path2plink=NULL,
                                  showPlinkOutput=TRUE, verbose=FALSE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    package.dir <- find.package('plinkQC')
    highld <- file.path(package.dir,'extdata', 'high-LD-regions.txt')
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (verbose) message(paste("Prune", prefix, "for relatedness estimation"))
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --exclude range ", highld,
                 " --indep-pairwise 50 5 0.2",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
    if (verbose) message("Run check_relatedness via plink --genome")
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --extract ", out, ".prune.in",
                 " --maf 0.1 --genome",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
    if (! file.exists(paste(prefix, ".imiss", sep=""))) {
        system(paste(path2plink, "plink --bfile ", prefix,
                     " --missing",
                     " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
    }
}

#' Evaluate results from PLINK IBD estimation.
#'
#' Evaluates and depicts results from plink --genome on the LD pruned dataset
#' (via \code{\link{run_check_relatedness}} or externally conducted IBD
#' estimation). plink --genome calculates identity by state (IBS) for each pair
#' of individuals based on the average proportion of alleles shared at genotyped
#' SNPs. The degree of recent shared ancestry, i.e. the identity by descent
#' (IBD) can be estimated from the genome-wide IBS. The proportion of IBD
#' between two individuals is returned by --genome as PI_HAT.
#' \code{evaluate_check_relatedness} finds pairs of samples whose proportion of
#' IBD is larger than the specified highIBDTh. Subsequently, for pairs of
#' individual that do not have additional relatives in the dataset, the
#' individual with the greater genotype missingness rate is selected and
#' returned as the individual failing the relatedness check. For more complex
#' family structures, the unrelated individuals per family are selected (e.g. in
#' a parents-offspring trio, the offspring will be marked as fail, while the
#' parents will be kept in the analysis).
#' \code{evaluate_check_relatedness} depicts all pair-wise IBD-estimates as
#' histograms stratified by value of PI_HAT.
#'
#' Both \code{\link{run_check_relatedness}} and
#' \code{\link{evaluate_check_relatedness}} can simply be invoked by
#' \code{\link{check_relatedness}}.
#'
#' For details on the output data.frame fail_high_IBD, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#genome}.
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' name.imiss and name.genome results as returned by plink --missing and plink
#' --genome.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.genome and name.imiss.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_IBD() via ggplot2::ggsave(p=p_IBD,
#' other_arguments) or pdf(outfile) print(p_IBD) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return a named [list] with i) fail_high_IBD containing a [data.frame] of
#' IIDs and FIDs of individuals who fail the IBDTh in columns
#' FID1 and IID1. In addition, the following columns are returned (as originally
#' obtained by plink --genome):
#' FID2 (Family ID for second sample), IID2 (Individual ID for second sample),
#' RT (Relationship type inferred from .fam/.ped file), EZ (IBD sharing expected
#' value, based on just .fam/.ped relationship), Z0 (P(IBD=0)), Z1 (P(IBD=1)),
#' Z2 (P(IBD=2)), PI_HAT (Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)), PHE
#' (Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)),
#' DST (IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)), PPC (IBS
#' binomial test), RATIO (HETHET : IBS0 SNP ratio (expected value 2)).
#' and ii) failIDs containing a [data.frame] with individual IDs [IID] and
#' family IDs [FID] of individuals failing the highIBDTh iii) p_IBD, a
#' ggplot2-object 'containing' all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT, which can be
#' shown by print(p_IBD).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' qcdir <- file.path(package.dir,'extdata')
#' name <- 'data'
#' relatednessQC <- evaluate_check_relatedness(qcdir=qcdir, name=name,
#' interactive=FALSE)
evaluate_check_relatedness <- function(qcdir, name, highIBDTh=0.1875,
                                       imissTh=0.03, interactive=FALSE,
                                       verbose=FALSE) {
    prefix <- paste(qcdir, "/", name, sep="")
    if (!file.exists(paste(prefix, ".imiss", sep=""))){
        stop("plink --missing output file: ", prefix,
             ".imiss does not exist.")
    }
    if (!file.exists(paste(prefix, ".genome",sep=""))){
        stop("plink --genome output file: ", prefix,
             ".genome does not exist.")
    }
    testNumerics(numbers=highIBDTh, positives=highIBDTh, proportions=highIBDTh)
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_imiss == names(imiss))) {
        stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    names_genome <- c("FID1", "IID1", "FID2", "IID2", "RT", "EZ", "Z0", "Z1",
                      "Z2", "PI_HAT", "PHE", "DST", "PPC", "RATIO")
    genome <- read.table(paste(prefix, ".genome", sep=""), header=TRUE,
                         as.is=TRUE)
    if (!all(names_genome == names(genome))) {
        stop("Header of ", prefix, ".genome is not correct. Was your
             file generated with plink --genome?")
    }

    fail_highIBD <- relatednessFilter(relatedness=genome,  otherCriterion=imiss,
                                      relatednessTh=highIBDTh,
                                      relatednessFID1="FID1",
                                      relatednessFID2="FID2",
                                      otherCriterionTh=imissTh,
                                      otherCriterionThDirection="gt",
                                      otherCriterionMeasure="F_MISS" )

    genome$PI_HAT_bin <- ifelse(genome$PI_HAT > 0.05, 0, 1)
    p_allPI_HAT <- ggplot(genome, aes_string('PI_HAT'))
    p_allPI_HAT <- p_allPI_HAT + geom_histogram(binwidth = 0.005,
                                                fill="#66a61e") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for all sample pairs") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(title=element_text(size=8))
    p_highPI_HAT <- ggplot(dplyr::filter_(genome, ~PI_HAT_bin == 0),
                           aes_string('PI_HAT'))
    p_highPI_HAT <- p_highPI_HAT + geom_histogram(binwidth = 0.005,
                                                  fill="#e6ab02") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for sample pairs with PI_HAT >0.1") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(title=element_text(size=8))
    p_histo <- cowplot::plot_grid(p_allPI_HAT, p_highPI_HAT)
    title <- cowplot::ggdraw() +
        cowplot::draw_label("Relatedness estimated as pairwise IBD (PI_HAT)")
    p_IBD <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                rel_heights = c(0.1, 1))
    if (interactive) {
        print(p_IBD)
    }
    return(list(fail_highIBD=fail_highIBD$relatednessFails,
                failIDs=fail_highIBD$failIDs, p_IBD=p_IBD))
}

#' Run PLINK principal component analysis
#'
#' Run plink --pca to calculate the principal components on merged genotypes
#' of the study and reference dataset.
#'
#' Both, \code{\link{run_check_ancestry}} and its evaluation by
#' \code{\link{evaluate_check_ancestry}} can simply be invoked by
#' \code{\link{check_ancestry}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files prefixMergedDataset.bim,prefixMergedDataset.fam and
#' prefixMergedDataset.bed.
#' @param qcdir [character] /path/to/directory to save
#' prefixMergedDataset.eigenvec as returned by plink --pca. User needs writing
#' permission to qcdir. Per default qcdir=indir.
#' @param prefixMergedDataset [character] Prefix of merged study and
#' reference data files, i.e. prefixMergedDataset.bed, prefixMergedDataset.bim,
#' prefixMergedDataset.fam.
#' @param path2plink [character] Absolute path to directory where external plink
#' software \url{https://www.cog-genomics.org/plink/1.9/} can be found, i.e.
#' plink should be accesible as path2plink/plink -h. If not
#' provided, assumed that PATH set-up works and plink will be found by
#' system("plink").
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir,'extdata')
#' qcdir <- tempdir()
#' prefixMergedDataset <- 'data.HapMapIII'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run <- run_check_ancestry(indir=indir, qcdir=qcdir, prefixMergedDataset)
#' }
run_check_ancestry <- function(indir, prefixMergedDataset,
                               qcdir=indir,
                               verbose=FALSE, path2plink=NULL,
                               showPlinkOutput=TRUE) {
    prefix <- paste(indir, "/", prefixMergedDataset, sep="")
    out <- paste(qcdir, "/", prefixMergedDataset, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
    if (!is.null(path2plink)) {
        path2plink <- paste(gsub("/$", "", path2plink), "/", sep="")
    }
    checkPlink(path2plink)
    if (verbose) message("Run check_ancestry via plink --pca")
    system(paste(path2plink, "plink --bfile ", prefix,
                 " --pca",
                 " --out ", out, sep=""),
           ignore.stdout=showPlinkOutput, ignore.stderr=!showPlinkOutput)
}

#' Evaluate results from PLINK PCA on combined study and reference data
#'
#' Evaluates and depicts results of plink --pca (via
#' \code{\link{run_check_ancestry}} or externally conducted pca) on merged
#' genotypes from individuals to be QCed and individuals of reference population
#' of known genotypes. Currently, check ancestry only supports automatic
#' selection of individuals of European descent. It uses information from
#' principal components 1 and 2 returned by plink --pca to find the center of
#' the European reference samples (mean(PC1_europeanRef), mean(PC2_europeanRef).
#' It computes the maximum Euclidean distance (maxDist) of the European
#' reference samples from this centre. All study samples whose Euclidean
#' distance from the centre falls outside the circle described by the radius
#' r=europeanTh* maxDist are considered non-European and their IDs are returned
#' as failing the ancestry check.
#' check_ancestry creates a scatter plot of PC1 versus PC2 colour-coded for
#' samples of the reference populations and the study population.
#'
#' Both \code{\link{run_check_ancestry}} and
#' \code{\link{evaluate_check_ancestry}} can simply be invoked by
#' \code{\link{check_ancestry}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' prefixMergedDataset.eigenvec results as returned by plink --pca. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param prefixMergedDataset [character] Prefix of merged dataset (study and
#' reference samples) used in plink --pca, resulting in
#' prefixMergedDataset.eigenvec.
#' @param europeanTh [double] Scaling factor of radius to be drawn around center
#' of European reference samples, with study samples inside this radius
#' considered to be of European descent and samples outside this radius of
#' non-European descent. The radius is computed as the maximum Euclidean distance
#' of European reference samples to the centre of European reference samples.
#' @param refSamples [data.frame] Dataframe with sample identifiers
#' [refSamplesIID] corresponding to IIDs in prefixMergedDataset.eigenvec and
#' population identifier [refSamplesPop] corresponding to population IDs
#' [refColorsPop] in refColorsfile/refColors. Either refSamples or
#' refSamplesFile have to be specified.
#' @param refColors [data.frame, optional] Dataframe with population IDs in
#' column [refColorsPop] and corresponding colour-code for PCA plot in column
#' [refColorsColor]. If not provided and is.null(refColorsFile) default colors
#' are used.
#' @param refSamplesFile [character] /path/to/File/with/reference samples. Needs
#' columns with sample identifiers [refSamplesIID] corresponding to IIDs in
#' prefixMergedDataset.eigenvec and population identifier [refSamplesPop]
#' corresponding to population IDs [refColorsPop] in refColorsfile/refColors.
#' @param refColorsFile [character, optional]
#' /path/to/File/with/Population/Colors cotaining population IDs in column
#' [refColorsPop] and corresponding colour-code for PCA plot in column
#' [refColorsColor].If not provided and is.null(refColors) default colors for
#' are used.
#' @param refSamplesIID [character] Column name of reference sample IDs in
#' refSamples/refSamplesFile.
#' @param refSamplesPop [character] Column name of reference sample population
#' IDs in refSamples/refSamplesFile.
#' @param refColorsColor [character] Column name of population colors in
#' refColors/refColorsFile
#' @param refColorsPop [character] Column name of reference sample population
#' IDs in refColors/refColorsFile.
#' @param studyColor [character] Colour to be used for study population if plot
#' is TRUE.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_ancestry) via ggplot2::ggsave(p=p_ancestry,
#' other_arguments) or pdf(outfile) print(p_ancestry) dev.off().
#' @return Named [list] with i) fail_ancestry, containing a [data.frame] with
#' FID and IID of non-European individuals and ii) p_ancestry, a ggplot2-object
#' 'containing' a scatter plot of PC1 versus PC2 colour-coded for samples of the
#' reference populations and the study population, which can be shown by
#' print(p_ancestry).
#' @export
#' @examples
#' package.dir <- find.package('plinkQC')
#' indir <- file.path(package.dir, 'extdata')
#' name <- "data"
#' fail_ancestry <- evaluate_check_ancestry(indir=indir, name=name,
#' refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(indir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE)

evaluate_check_ancestry <- function(indir, name, prefixMergedDataset,
                                    qcdir=indir,
                                    europeanTh=1.5,
                                    refSamples=NULL, refColors=NULL,
                                    refSamplesFile=NULL, refColorsFile=NULL,
                                    refSamplesIID="IID", refSamplesPop="Pop",
                                    refColorsColor="Color", refColorsPop="Pop",
                                    studyColor="#2c7bb6",
                                    interactive=FALSE) {
    prefix <- paste(indir, "/", name, sep="")
    out <- paste(qcdir, "/", prefixMergedDataset, sep="")
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    samples <- data.table::fread(paste(prefix, ".fam", sep=""),
                                 header=FALSE, stringsAsFactors=FALSE,
                                 data.table=FALSE)[,1:2]
    colnames(samples) <- c("FID", "IID")
    if (!file.exists(paste(out, ".eigenvec", sep=""))){
        stop("plink --pca output file: ", out, ".eigenvec does not exist.")
    }
    testNumerics(numbers=europeanTh, positives=europeanTh)
    pca_data <- data.table::fread(paste(out, ".eigenvec", sep=""),
                                  stringsAsFactors=FALSE, data.table=FALSE)
    colnames(pca_data) <- c("FID", "IID", paste("PC",1:(ncol(pca_data)-2),
                                                sep=""))

    if (is.null(refSamples) && !file.exists(refSamplesFile)) {
        stop("refSamples are not specified and refSamplesFile file",
             refSamplesFile, "does not exist.")
    }
    if (!is.null(refSamplesFile)) {
        refSamples <- read.table(refSamplesFile, header=TRUE,
                                 stringsAsFactors=FALSE)
    }
    if (!(refSamplesIID  %in% names(refSamples))) {
        stop(paste("Column", refSamplesIID, "not found in refSamples."))
    }
    if (!(refSamplesPop %in% names(refSamples))) {
        stop(paste("Column", refSamplesPop, "not found in refSamples."))
    }
    names(refSamples)[names(refSamples) == refSamplesIID] <- "IID"
    names(refSamples)[names(refSamples) == refSamplesPop] <- "Pop"
    refSamples <- dplyr::select_(refSamples, ~IID, ~Pop)

    if (!is.null(refColorsFile) && !file.exists(refColorsFile)) {
        stop("refColorsFile file", refColorsFile, "does not exist.")
    }
    if (!is.null(refColorsFile)) {
        refColors <- read.table(refColorsFile, header=TRUE,
                                stringsAsFactors=FALSE)
    }
    if (!is.null(refColors)) {
        if (!(refColorsColor  %in% names(refColors))) {
            stop(paste("Column", refColorsColor, "not found in refColors."))
        }
        if (!(refColorsPop %in% names(refColors))) {
            stop(paste("Column", refColorsPop, "not found in refColors."))
        }
        names(refColors)[names(refColors) == refColorsColor] <- "Color"
        names(refColors)[names(refColors) == refColorsPop] <- "Pop"
        refColors <- dplyr::select_(refColors, ~Pop, ~Color)
    } else {
        refColors <- data.frame(Pop=unique(as.character(refSamples$Pop)))
        refColors$Color <- 1:nrow(refColors)
    }
    refSamples <- merge(refSamples, refColors, by="Pop", all.X=TRUE)

    ## Combine pca data and population information ####
    data_all <- merge(pca_data, refSamples, by="IID", all.x=TRUE)
    data_all$Pop[is.na(data_all$Pop)] <- name
    data_all$Color[is.na(data_all$Color)] <- studyColor
    data_all <- data_all[order(data_all$Pop, decreasing=FALSE),]

    refColors <- rbind(refColors, c(name, studyColor))
    data_all$Color <- as.factor(data_all$Color)
    data_all$Pop <- factor(data_all$Pop, levels=refColors$Pop)

    ## Find mean coordinates and distances of reference Europeans ####
    all_european <- dplyr::filter_(data_all, ~Pop %in% c("CEU", "TSI"))
    euro_pc1_mean <- mean(all_european$PC1)
    euro_pc2_mean <- mean(all_european$PC2)

    all_european$euclid_dist <- sqrt((all_european$PC1 - euro_pc1_mean)^2 +
                                         (all_european$PC2 - euro_pc2_mean)^2)

    max_euclid_dist <- max(all_european$euclid_dist)

    ## Find samples' distances to reference Europeans ####
    data_name <- dplyr::filter_(data_all, ~Pop == name)
    data_name$euclid_dist <- sqrt((data_name$PC1 - euro_pc1_mean)^2 +
                                     (data_name$PC2 - euro_pc2_mean)^2)
    non_europeans <- dplyr::filter_(data_name, ~euclid_dist >
                                        (max_euclid_dist * europeanTh))
    fail_ancestry <- dplyr::select_(non_europeans, ~FID, ~IID)
    p_ancestry <- ggplot()
    p_ancestry <- p_ancestry + geom_point(data=data_all,
                                          aes_string(x='PC1', y='PC2',
                                                     color='Pop')) +
        geom_point(data=dplyr::filter_(data_all, ~Pop != name),
                   aes_string(x='PC1', y='PC2',
                              color='Pop'),
                   size=1) +
        scale_color_manual(values=refColors$Color,
                           name="Population") +
        guides(color=guide_legend(nrow=2, byrow=TRUE)) +
        ggforce::geom_circle(aes(x0=euro_pc1_mean, y0=euro_pc2_mean,
                                 r=(max_euclid_dist * europeanTh))) +
        ggtitle("PCA on combined reference and study genotypes") +
        theme_bw() +
        theme(legend.position='bottom')
    if (interactive) {
        print(p_ancestry)
    }
    return(list(fail_ancestry=fail_ancestry, p_ancestry=p_ancestry))
}


