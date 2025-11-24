context('Test applyQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'

fail_individuals <-
    perIndividualQC(qcdir=qcdir, indir=indir, name=name,
                    do.run_check_sex=FALSE,
                    do.run_check_het_and_miss=FALSE,
                    do.run_check_relatedness=FALSE,
                    interactive=FALSE, verbose=FALSE,
                    dont.ancestry_prediction = TRUE)

context('Test cleanData')
test_that('cleanData throws file error',{
    expect_error(cleanData(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})
test_that('cleanData fails with missing check error', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterRelated=FALSE,
                           filterSex=FALSE, filterHeterozygosity=FALSE,
                           filterAncestry=FALSE,
                           filterSampleMissingness=FALSE,
                           filterSNPMissingness=FALSE, filterHWE=FALSE,
                           filterMAF=FALSE),
                 "No per-sample and per-marker filters chosen")
          })


test_that('cleanData returns message with missing sample check', {
    expect_error(cleanData(qcdir=qcdir, indir=indir, name=name,
                           filterRelated=FALSE, filterAncestry=FALSE,
                           filterSampleMissingness=FALSE,
                           filterSex=FALSE, filterHeterozygosity=FALSE,
                           filterMAF=FALSE, filterSNPMissingness=FALSE,
                           filterHWE=FALSE),
                 "No per-sample and per-marker filters chosen")
})



