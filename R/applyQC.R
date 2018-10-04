#' Create plink dataset with samples and markers passing quality control
#'
#' Remove samples that fail per individual QC
#' Remove markers that fail per marker QC
clean_data <- function(qcdir, alg, macTh=20, mafTh=0.01, highIBDTh=NULL) {
    if (!is.null(highIBDTh)) {
        # i) with sample filtering based on all fail-IDs (including relatedness)
        # and variant filter for missing genotypes, MAF and HWE
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir,"/", alg, ".fail.IDs ",
                     "--hwe ", hweTh,
                     " --geno ", lmissTh,
                     " --maf ", MAF ,
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
    system(paste("plink --bfile ", qcdir, "/", alg,
                 " --remove ", qcdir,"/", alg, ".noIBDfail.IDs ",
                 "--hwe ", hweTh,
                 " --geno ", lmissTh,
                 " --make-bed --out ", qcdir, "/", alg, ".clean.related",
                 sep=""))
    return(MAF)
}
