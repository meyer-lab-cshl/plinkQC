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

test_that('evaluate_check_ancestry uses HapMap as default reference',{
    expect_message(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         verbose=TRUE),
                   "Using HapMap as reference samples.")
})

test_that('evaluate_check_ancestry uses HapMap as default reference',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                           prefixMergedDataset=prefix,
                                           verbose=TRUE,
                                           defaultRefSamples = "Genomes"),
                   "defaultRefSamples should be one of 'HapMap'.")
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

test_that('evaluate_check_ancestry throws missing refPopulation code error',{
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         refPopulation = c("CCC")),
                 "Not all refPopulation populations")
})

test_that('evaluate_check_ancestry throws missing missing ID error',{
    testRefs <- refSamples[84:88,]
    expect_error(evaluate_check_ancestry(qcdir, name,
                                         prefixMergedDataset=prefix,
                                         refSamples=testRefs,
                                         refPopulation = "CEU"),
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
    expect_equal(length(fail), 3)
})

test_that('evaluate_check_ancestry returns correct output names',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                               refSamples=refSamples)
    expect_equal(names(fail), c("fail_ancestry", "p_ancestry", "plot_data"))
})

test_that('evaluate_check_ancestry returns correct fail IDs for example data',{
    fail <- evaluate_check_ancestry(qcdir, name, prefixMergedDataset=prefix,
                                    refSamples=refSamples)
    expect_true(all(fail$fail_ancestry$IID %in% fail_ancestryIDs[,1]))
})

test_that('evaluate_check_ancestry fails with all missing sample error',{
    expect_error(evaluate_check_ancestry(
        qcdir, name, prefixMergedDataset="data.HapMapIII_no_samples",
        refSamples=refSamples),
        "There are no")
})

test_that('evaluate_check_ancestry fails with some missing sample error',{
    expect_error(evaluate_check_ancestry(
        qcdir, name, prefixMergedDataset="data.HapMapIII_some_samples",
        refSamples=refSamples),
        "Not all")
})

test_that('evaluate_check_ancestry fails with additional sample in highlighting',{
    expect_error(evaluate_check_ancestry(qcdir, name=name,
                                         highlight_samples = "XBC",
                                         interactive=FALSE,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples
    ),
    "Not all samples to be highlighted")
})

test_that('evaluate_check_ancestry fails with text and label highlighting',{
    expect_error(evaluate_check_ancestry(qcdir, name=name,
                                         prefixMergedDataset=prefix,
                                         refSamples=refSamples,
                                         highlight_samples = highlight_samples,
                                         highlight_type = c("text", "label"),
                                         interactive=FALSE
    ),
    "Only one of text or label")
})


test_that('perIndividualQC works if all samples pass', {
    fail_individuals <-
        perIndividualQC(indir=indir,
                        qcdir=qcdir, name="data.clean",
                        refSamplesFile=paste0(qcdir, "/HapMap_ID2Pop.txt"),
                        refColorsFile=paste0(qcdir, "/HapMap_PopColors.txt"),
                        prefixMergedDataset="data.clean.HapMapIII",
                        interactive=FALSE, verbose=FALSE,
                        dont.check_ancestry = TRUE,
                        do.run_check_sex = FALSE,
                        do.run_check_relatedness = FALSE,
                        do.run_check_het_and_miss = FALSE)
    expect_equal(unlist(fail_individuals$fail_list), NULL)
})
