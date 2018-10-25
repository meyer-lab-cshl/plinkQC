context('Test markerQC functions')
package.dir <- find.package('plinkQC')
#package.dir <- "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/plinkQC"
qcdir <- file.path(package.dir, 'extdata')
alg <- 'data'


context('Test check_snp_missingness')
test_that('check_snp_missingness throws file error',{
    expect_error(check_snp_missingness(qcdir, "nodata", verbose=FALSE),
                 "plink family file")
})

context('Test check_hwe')
test_that('check_hwe throws file error',{
    expect_error(check_hwe(qcdir, "nodata", verbose=FALSE),
                 "plink family file")
})

context('Test check_maf')
test_that('check_maf throws file error',{
    expect_error(check_maf(qcdir, "nodata", verbose=FALSE),
                 "plink family file")
})

test_that('overviewPerMarkerQC returns QC fails message', {
    results_perMarkerQC <- list(fail_list=list(a=1:10), p_markerQC=NULL)
    expect_message(overviewPerMarkerQC(results_perMarkerQC, interactive=TRUE),
                    "overviewMarkerQC cannot be displayed with UpSetR")
})
