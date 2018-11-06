context('Test individualQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'
prefix <- 'data.HapMapIII'
refSamplesFile <-  "./HapMap_ID2Pop.txt"
refColorsFile <-  "./HapMap_PopColors.txt"

fail_sexIDs <- read.table(paste(qcdir, '/', name, ".fail-sexcheck.IDs", sep=""))
fail_highIBDIDs <- read.table(paste(qcdir, '/', name, ".fail-IBD.IDs", sep=""))
fail_imissIDs <- read.table(paste(qcdir, '/', name, ".fail-imiss.IDs", sep=""))
fail_hetIDs <- read.table(paste(qcdir, '/', name, ".fail-het.IDs", sep=""))
fail_ancestryIDs <- read.table(paste(qcdir, '/', name, ".fail-ancestry.IDs",
                                     sep=""))
refSamples <- read.table(paste(qcdir, '/', refSamplesFile, sep=""), header=TRUE,
                         stringsAsFactors=FALSE)
refColors <- read.table(paste(qcdir, '/', refColorsFile, sep=""), header=TRUE,
                         stringsAsFactors=FALSE)

context('Test evaluate_check_sex')
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

context('Test evaluate_check_ancestry')
test_that('evaluate_check_ancestry throws sample file error',{
    expect_error(evaluate_check_ancestry(qcdir, "nodata",
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples),
                 "plink family file:")
})

test_that('evaluate_check_ancestry throws input file error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset='nodata',
                                         refSamples=refSamples),
                 "plink --pca output file:")
})

test_that('evaluate_check_ancestry throws refSamples missing error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix),
                 "Neither refSamples nor refSamplesFile are specified")
})

test_that('evaluate_check_ancestry throws refSamples file error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamplesFile="nofile"),
                 "refSamplesFile file")
})

test_that('evaluate_check_ancestry throws refSamples IID column error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refSamplesIID="none"),
                 "Column")
})

test_that('evaluate_check_ancestry throws refSamples Pop column error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refSamplesPop="none"),
                 "Column")
})

test_that('evaluate_check_ancestry throws refColors file error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refColorsFile="nofile"),
                 "refColorsFile file")
})

test_that('evaluate_check_ancestry throws refColors Pop column error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refColors=refColors,
                                         refColorsPop="none"),
                 "Column")
})

test_that('evaluate_check_ancestry throws refColors Color column error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refColors=refColors,
                                         refColorsColor="none"),
                 "Column")
})


test_that('evaluate_check_ancestry throws missing population code error',{
    testColors <- refColors[1:4,]
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refColors=testColors),
                 "Not all refSamples populations")
})

test_that('evaluate_check_ancestry throws missing missing ID error',{
    testRefs <- refSamples[1:4,]
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=testRefs),
                 "There are samples in the prefixMergedDataset")
})

test_that('evaluate_check_ancestry fails with numbers error', {
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         europeanTh=-1,
                                         refSamples=refSamples),
                 "has/have to be greater than zero")
})

test_that('evaluate_check_ancestry returns correct output type',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                                    refSamples=refSamples)
    expect_true(is.list(fail))
})

test_that('evaluate_check_ancestry returns correct output length',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                                        refSamples=refSamples)
    expect_equal(length(fail), 2)
})

test_that('evaluate_check_ancestry returns correct output names',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                               refSamples=refSamples)
    expect_equal(names(fail), c("fail_ancestry", "p_ancestry"))
})

test_that('evaluate_check_ancestry returns correct fail IDs for example data',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                                    refSamples=refSamples)
    expect_true(all(fail$fail_ancestry$IID %in% fail_ancestryIDs[,1]))
})

context('Test overviewPerIndividualQC')
test_that('overviewPerIndividualQC returns QC fails message', {
    results_perSampleQC <- list(fail_list=list(a=1:10), p_SampleQC=NULL)
    expect_message(overviewPerIndividualQC(results_perSampleQC,
                                           interactive=TRUE),
                   "overviewSampleQC for QC fails cannot be displayed with")
})
