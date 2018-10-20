context('Test sampleQC functions')
# package.dir <- find.package('plinkQC')
package.dir <- "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/plinkQC"
qcdir <- file.path(package.dir, 'extdata')
alg <- 'data'


fail_sexIDs <- read.table(paste(qcdir, "/", alg, ".fail-sexcheck.IDs", sep=""))
fail_highIBDIDs <- read.table(paste(qcdir, "/", alg, ".fail-IBD.IDs", sep=""))
fail_imissIDs <- read.table(paste(qcdir, "/", alg, ".fail-imiss.IDs", sep=""))
fail_hetIDs <- read.table(paste(qcdir, "/", alg, ".fail-het.IDs", sep=""))

context('Test check_sex')
test_that('check_sex throws file error',{
    expect_error(check_sex(qcdir, "nodata", verbose=FALSE),
                 "plink --check-sex results file")
})
test_that('check_sex fails with femaleTh/maleTh numbers error', {
              expect_error(check_sex(qcdir, alg, maleTh=2, femaleTh=3),
                           "Proportions have to be specified")
          })
#test_that('check_sex throws path2plink error',{
#    expect_error(check_sex(qcdir, alg, path2plink="~/test", verbose=FALSE),
#                 "PLINK software required for running this ")
#})
#
test_that('check_sex returns correct output type',{
    fail_sex <- check_sex(qcdir, alg, verbose=FALSE)
    expect_true(is.list(fail_sex))
})

test_that('check_sex returns correct output length',{
    fail_sex <- check_sex(qcdir, alg, verbose=FALSE)
    expect_equal(length(fail_sex), 3)
})

test_that('check_sex returns correct output names',{
    fail_sex <- check_sex(qcdir, alg, verbose=FALSE)
    expect_equal(names(fail_sex), c("fail_sex", "mixup", "p_sexcheck"))
})

test_that('check_sex returns correct fail IDs for example data',{
    fail_sex <- check_sex(qcdir, alg, verbose=FALSE)
    expect_true(all(fail_sex$fail_sex$IID %in% fail_sexIDs[,1]))
})


context('Test check_heterozygosity_and_missingness')
test_that('check_heterozygosity_and_missingness fails with hetTh numbers error',
          {
    expect_error(check_heterozygosity_and_missingness(qcdir, alg,
                                                              hetTh=-1),
                 "has/have to be greater than zero")
})

test_that('check_heterozygosity_and_missingness fails with imissTh error',{
    expect_error(check_heterozygosity_and_missingness(qcdir, alg,
                                                      imissTh=2),
                 "Proportions have to be specified")
})

test_that('check_heterozygosity_and_missingness returns correct output type',{
    fail_het_imiss <- check_heterozygosity_and_missingness(qcdir, alg)
    expect_true(is.list(fail_het_imiss))
})

test_that('check_heterozygosity_and_missingness returns correct output length',{
    fail_het_imiss <- check_heterozygosity_and_missingness(qcdir, alg)
    expect_equal(length(fail_het_imiss), 3)
})

test_that('check_heterozygosity_and_missingness returns correct output names',{
    fail_het_imiss <- check_heterozygosity_and_missingness(qcdir, alg)
    expect_equal(names(fail_het_imiss), c("fail_imiss", "fail_het",
                                          "p_het_imiss"))
})

test_that('check_heterozygosity_and_missingness returns correct fail IDs',{
    fail_het_imiss <- check_heterozygosity_and_missingness(qcdir, alg)
    expect_true(all(fail_het_imiss$fail_het$IID %in% fail_hetIDs[,1]))
    expect_true(all(fail_het_imiss$fail_imiss$IID %in% fail_imissIDs[,1]))
})

context('Test check_relatedness')

test_that('check_relatedness throws file error',{
    expect_error(check_relatedness(qcdir, "nodata"),
                 "plink --missing output file")
})
test_that('check_relatedness fails with highIBDTh numbers error', {
    expect_error(check_relatedness(qcdir, alg, highIBDTh=2),
                 "Proportions have to be specified")
})

test_that('check_relatedness returns correct output type',{
    fail_relatedness <- check_relatedness(qcdir, alg, verbose=FALSE)
    expect_true(is.list(fail_relatedness))
})

test_that('check_relatedness returns correct output length',{
    fail_relatedness <- check_relatedness(qcdir, alg, verbose=FALSE)
    expect_equal(length(fail_relatedness), 3)
})

test_that('check_relatedness returns correct output names',{
    fail_relatedness <- check_relatedness(qcdir, alg, verbose=FALSE)
    expect_equal(names(fail_relatedness), c("fail_highIBD", "failIDs", "p_IBD"))
})

test_that('check_relatedness returns correct fail IDs for example data',{
    fail_relatedness <- check_relatedness(qcdir, alg, verbose=FALSE)
    expect_true(all(fail_relatedness$fail_highIBD$IID %in%
                        fail_highIBDIDs[,1]))
})

test_that('overviewPerSampleQC returns QC fails message', {
    results_perSampleQC <- list(fail_list=list(a=1:10), p_SampleQC=NULL)
    expect_message(overviewPerSampleQC(results_perSampleQC, interactive=TRUE),
                   "overviewSampleQC for QC fails cannot be displayed with")
})
