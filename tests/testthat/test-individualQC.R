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
fail_ancestryIDs <- read.table(paste(qcdir, '/', name, ".exclude-ancestry.IDs",
                                     sep=""))
refSamples <- read.table(paste(qcdir, '/', refSamplesFile, sep=""), header=TRUE,
                         stringsAsFactors=FALSE)
refColors <- read.table(paste(qcdir, '/', refColorsFile, sep=""), header=TRUE,
                         stringsAsFactors=FALSE)

highlight_samples <- paste0("ID_", 1:10)

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
    expect_equal(length(fail_sex), 4)
})

test_that('evaluate_check_sex returns correct output names',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_equal(names(fail_sex), c("fail_sex", "mixup", "p_sexcheck",
                                    "plot_data"))
})

test_that('evaluate_check_sex returns correct fail IDs for example data',{
    fail_sex <- evaluate_check_sex(qcdir, name, verbose=FALSE)
    expect_true(all(fail_sex$fail_sex$IID %in% fail_sexIDs[,1]))
})

test_that('evaluate_check_sex runs successfully if no sample fails',{
    fail_sex <- evaluate_check_sex(qcdir, name="data_all_passing",
                                   verbose=FALSE, interactive=TRUE)
    expect_true(is.null(fail_sex$fail_sex))
})

test_that('evaluate_check_sex fails with additional sample in highlighting',{
    expect_error(evaluate_check_sex(qcdir, name=name,
                                    highlight_samples = "XBC",
                                    verbose=FALSE, interactive=FALSE
                                    ),
    "Not all samples to be highlighted")
})

test_that('evaluate_check_sex fails with text and label highlighting',{
    expect_error(evaluate_check_sex(qcdir, name=name,
                                    highlight_samples = highlight_samples,
                                    highlight_type = c("text", "label"),
                                    verbose=FALSE, interactive=FALSE
    ),
    "Only one of text or label")
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
    expect_equal(length(fail_het_imiss), 4)
})

test_that('evaluate_check_het_and_miss returns correct output names',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_equal(names(fail_het_imiss), c("fail_imiss", "fail_het",
                                          "p_het_imiss", "plot_data"))
})

test_that('evaluate_check_het_and_miss returns correct fail IDs',{
    fail_het_imiss <- evaluate_check_het_and_miss(qcdir, name)
    expect_true(all(fail_het_imiss$fail_het$IID %in% fail_hetIDs[,1]))
    expect_true(all(fail_het_imiss$fail_imiss$IID %in% fail_imissIDs[,1]))
})

test_that('evaluate_check_het_and_miss runs successfully if no sample fails',{
    fail_check_het_and_miss <-
        evaluate_check_het_and_miss(qcdir, name="data_all_passing",
                                    interactive=TRUE)
    expect_true(is.null(fail_check_het_and_miss$fail_het_and_miss))
})

test_that('evaluate_check_het_and_miss fails with additional sample in highlighting',{
    expect_error(evaluate_check_het_and_miss(qcdir, name=name,
                                    highlight_samples = "XBC", interactive=FALSE
    ),
    "Not all samples to be highlighted")
})

test_that('evaluate_check_het_and_miss fails with text and label highlighting',{
    expect_error(evaluate_check_het_and_miss(qcdir, name=name,
                                    highlight_samples = highlight_samples,
                                    highlight_type = c("text", "label"),
                                    interactive=FALSE
    ),
    "Only one of text or label")
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
    expect_equal(length(fail_relatedness), 4)
})

test_that('evaluate_check_relatedness returns correct output names',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_equal(names(fail_relatedness), c("fail_highIBD", "failIDs", "p_IBD",
                                            "plot_data"))
})

test_that('evaluate_check_relatedness returns correct fail IDs for example data',{
    fail_relatedness <- evaluate_check_relatedness(qcdir, name, verbose=FALSE)
    expect_true(all(fail_relatedness$fail_highIBD$IID %in%
                        fail_highIBDIDs[,1]))
})


test_that('perIndividualQC works if all samples pass', {
    fail_individuals <-
        perIndividualQC(indir=indir,
                        qcdir=qcdir, name="data.clean",
                        interactive=FALSE, verbose=FALSE,
                        do.run_check_sex = FALSE,
                        do.run_check_relatedness = FALSE,
                        do.run_check_het_and_miss = FALSE,
                        do.run_superpop_classification = FALSE)
    expect_equal(unlist(fail_individuals$fail_list), NULL)
})
