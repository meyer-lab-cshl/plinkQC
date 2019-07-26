context('Test markerQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'


context('Test check_snp_missingness')
test_that('check_snp_missingness throws file error',{
    expect_error(check_snp_missingness(qcdir=qcdir, indir=indir, name="nodata",
                                       verbose=FALSE),
                 "plink family file")
})

context('Test check_hwe')
test_that('check_hwe throws file error',{
    expect_error(check_hwe(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})

context('Test check_maf')
test_that('check_maf throws file error',{
    expect_error(check_maf(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})
