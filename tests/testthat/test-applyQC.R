context('Test applyQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'

fail_individuals <-
    perIndividualQC(qcdir=qcdir, indir=indir, name=name,
                    refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt", sep=""),
                    refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
                    do.run_check_sex=FALSE, do.run_check_ancestry=FALSE,
                    do.run_check_het_and_miss=FALSE,
                    do.run_check_relatedness=FALSE,
                    prefixMergedDataset="data.HapMapIII",
                    interactive=FALSE, verbose=FALSE)

context('Test cleanData')
test_that('cleanData throws file error',{
    expect_error(cleanData(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})
test_that('cleanData fails with missing check error', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterAncestry=FALSE,
                           filterRelated=FALSE,
                           filterSex=FALSE, filterHeterozygosity=FALSE,
                           filterSampleMissingness=FALSE,
                           filterSNPMissingness=FALSE, filterHWE=FALSE,
                           filterMAF=FALSE),
                 "No per-sample and per-marker filters chosen")
          })

test_that('cleanData fails with missing hweTh error', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterHWE=TRUE,
                           hweTh=NULL),
                 "filterHWE is TRUE but hweTh not specified")
})

test_that('cleanData fails with missing mafTh and macTh error', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterMAF=TRUE,
                           mafTh=NULL, macTh=NULL),
                 "filterMAF is TRUE but neither mafTh or macTh")
})

test_that('cleanData fails with missing lmissTh error', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterSNPMissingness=TRUE,
                           lmissTh=NULL),
                 "filterSNPMissingness is TRUE but lmissTh not specified")
})

test_that('cleanData returns message with missing sample check', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterAncestry=FALSE,
                           filterRelated=FALSE, filterSampleMissingness=FALSE,
                           filterSex=FALSE, filterHeterozygosity=FALSE,
                           filterMAF=FALSE, filterSNPMissingness=FALSE,
                           filterHWE=FALSE),
                 "No per-sample and per-marker filters chosen")
})



