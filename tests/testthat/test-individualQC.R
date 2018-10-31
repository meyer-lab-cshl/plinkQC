context('Test individualQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'


fail_sexIDs <- read.table(paste(name, ".fail-sexcheck.IDs", sep=""))
fail_highIBDIDs <- read.table(paste(name, ".fail-IBD.IDs", sep=""))
fail_imissIDs <- read.table(paste(name, ".fail-imiss.IDs", sep=""))
fail_hetIDs <- read.table(paste(name, ".fail-het.IDs", sep=""))

context('Test check_sex')
test_that('check_sex throws file error',{
    expect_error(evaluate_check_sex(qcdir, "nodata", verbose=FALSE),
                 "plink --check-sex results file")
})
test_that('evaluate_check_sex fails with femaleTh/maleTh numbers error', {
              expect_error(evaluate_check_sex(qcdir, name, maleTh=2, femaleTh=3),
                           "Proportions have to be specified")
          })

test_that('evaluate_check_sex returns correct output type',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_true(is.list(fail_sex))
})

test_that('evaluate_check_sex returns correct output length',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_equal(length(fail_sex), 3)
})

test_that('evaluate_check_sex returns correct output names',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_equal(names(fail_sex), c("fail_sex", "mixup", "p_sexcheck"))
})

test_that('evaluate_check_sex returns correct fail IDs for example data',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_true(all(fail_sex$fail_sex$IID %in% fail_sexIDs[,1]))
})


context('Test evaluate_check_het_and_miss')
test_that('evaluate_check_het_and_miss fails with hetTh numbers error',
          {
    expect_error(evaluate_check_het_and_miss(qcdir, name, hetTh=-1),
                 "has/have to be greater than zero")
})

test_that('evaluate_check_het_and_miss fails with imissTh error',{
    expect_error(evaluate_check_het_and_miss(qcdir, name, imissTh=2),
                 "Proportions have to be specified")
})

test_that('evaluate_check_het_and_miss returns correct output type',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_true(is.list(fail_het_imiss))
})

test_that('evaluate_check_het_and_miss returns correct output length',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_equal(length(fail_het_imiss), 3)
})

test_that('evaluate_check_het_and_miss returns correct output names',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_equal(names(fail_het_imiss), c("fail_imiss", "fail_het",
                                          "p_het_imiss"))
})

test_that('evaluate_check_het_and_miss returns correct fail IDs',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_true(all(fail_het_imiss$fail_het$IID %in% fail_hetIDs[,1]))
    expect_true(all(fail_het_imiss$fail_imiss$IID %in% fail_imissIDs[,1]))
})

context('Test evaluate_check_relatedness')

test_that('evaluate_check_relatedness throws file error',{
    expect_error(evaluate_check_relatedness(qcdir, "nodata"),
                 "plink --missing output file")
})
test_that('evaluate_check_relatedness fails with highIBDTh numbers error', {
    expect_error(evaluate_check_relatedness(qcdir, name, highIBDTh=2),
                 "Proportions have to be specified")
})

test_that('evaluate_check_relatedness returns correct output type',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_true(is.list(fail_relatedness))
})

test_that('evaluate_check_relatedness returns correct output length',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_equal(length(fail_relatedness), 3)
})

test_that('evaluate_check_relatedness returns correct output names',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_equal(names(fail_relatedness), c("fail_highIBD", "failIDs", "p_IBD"))
})

test_that('evaluate_check_relatedness returns correct fail IDs for example data',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_true(all(fail_relatedness$fail_highIBD$IID %in%
                        fail_highIBDIDs[,1]))
})

test_that('overviewPerIndividualQC returns QC fails message', {
    results_perSampleQC <- list(fail_list=list(a=1:10), p_SampleQC=NULL)
    expect_message(overviewPerIndividualQC(results_perSampleQC,
                                           interactive=TRUE),
                   "overviewSampleQC for QC fails cannot be displayed with")
})
