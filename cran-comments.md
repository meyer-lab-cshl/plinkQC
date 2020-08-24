## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release, addressing some minor bug fixes and adding some more
functionality.

# plinkQC 0.3.3
## major changes
* Output all IBD estimates in `run_check_relatedness` not just fail IBD
estimates [d009da90](https://github.com/meyer-lab/plinkQC/commit/d009da90434a83cb4867881e3772b518d17dc2dd)
* Add midp modifier for --hwe and --hardy as recommended on [plink hwe documentation](https://www.cog-genomics.org/plink/1.9/filter#hwe) and discussed in issue [#27](https://github.com/meyer-lab-cshl/plinkQC/issues/27); [43230a6c](https://github.com/meyer-lab/plinkQC/commit/43230a6ce5f079b46c85b1a7a526213177db1f1b)

## minor changes
* Add checks in `cleanData` that can distinguish between miss-sepcified filters and non-failing samples [7432b1ee](https://github.com/meyer-lab/plinkQC/commit/7432b1ee880ad14c7a078f3a0379f7fdcff1684b); addresses issue [#31](https://github.com/meyer-lab-cshl/plinkQC/issues/31) 
* Add checks and tests in `check_snp_missingness`, `check_maf` and `check_hwe` 
in cases where all samples pass `perIndividualQC`.
[31ed6f65](https://github.com/meyer-lab/plinkQC/commit/31ed6f653902857a0a78ae723e786dd84d240b75); related to issue [#31](https://github.com/meyer-lab-cshl/plinkQC/issues/31) 
* Add checks and tests in `evaluate_check_sex` and `evaluate_check_het_imiss`
for all passing samples and accordingly, no fail samples label in plot.
[d071a43](https://github.com/meyer-lab/plinkQC/commit/d071a4341b0d4c40220b42f81c412da4d320b8ec)
