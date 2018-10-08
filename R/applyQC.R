#' Create plink dataset with samples and markers passing quality control
#'
#' Remove samples that fail per individual QC
#' Remove markers that fail per marker QC
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param filterRelated [logical]
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specified, macTh is used (macTh = mafTh\*2\*NrSamples)
#' @param path2plink [character] Absolute path to where external plink software
#' \url{https://www.cog-genomics.org/plink/1.9/} can be found. If not provided,
#' assumed that PATH set-up works and plink will be found by system("plink").
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export

clean_data <- function(qcdir, alg, lmissTh=0.01, hweTh=1e-5, macTh=20,
                       mafTh=0.01, filterRelated=FALSE, path2plink=NULL,
                       verbose=FALSE) {
    checkPlink(path2plink)
    if (!is.null(path2plink)) paste(gsub("/$", "", path2plink), "/", sep="")
    if (filterRelated) {
        # i) with sample filtering based on all fail-IDs (including relatedness)
        # and variant filter for missing genotypes, MAF and HWE
        system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir,"/", alg, ".fail.IDs",
                     " --hwe ", hweTh,
                     " --geno ", lmissTh,
                     " --maf ", mafTh ,
                     " --make-bed --out ", qcdir, "/", alg, ".clean", sep=""))
    }
    # ii) with sample filtering based on fail-IDs not including failed IBD and
    # variant filtering for missing genotypes and HWE (keep related samples and
    # MAF since only part of the cohort; smart relatedness filtering when
    # combining cohorts)
    system(paste("cat ", qcdir ,"/",  alg, ".fail-ancestry.IDs  ",
                 qcdir ,"/",  alg, ".fail-het.IDs ",
                 qcdir ,"/",  alg, ".fail-sexcheck.IDs  ",
                 qcdir ,"/",  alg, ".fail-imiss.IDs | sort | uniq >",
                 qcdir, "/", alg, ".noIBDfail.IDs", sep=""), wait=TRUE)
    system(paste(path2plink, "plink --bfile ", qcdir, "/", alg,
                 " --remove ", qcdir,"/", alg, ".noIBDfail.IDs",
                 " --hwe ", hweTh,
                 " --geno ", lmissTh,
                 " --make-bed --out ", qcdir, "/", alg, ".clean.related",
                 sep=""))
}
