## Resubmission
# plinkQC 0.3.2
## minor changes
* Add checks and tests in `evaluate_check_ancestry` for missing non-reference samples
[5c03971](https://github.com/meyer-lab/plinkQC/commit/5c03971b92341b5a8cb26a1d5f2e37d637f5cbd7)
* Clarified tutorials:
  * R internal file.copy instead of systems `cp` command to make run on windows: [6c91ef3](https://github.com/meyer-lab/plinkQC/commit/6c91ef3368f601658967f88bc20db94a28307fce),
  fixes [#26](https://github.com/meyer-lab/plinkQC/issues/26) 
  * Add additional filter for A->T, C->G variants as suggested in [#24](https://github.com/meyer-lab/plinkQC/issues/24); [11e0375](https://github.com/meyer-lab/plinkQC/commit/11e03756852215915310a836d2e73dbcedf12bfc)
* Fixed plotting issues in PCA plot of ancestry check:
  * allow for supplying names of European reference population [ae09e64](https://github.com/meyer-lab/plinkQC/commit/ae09e64ee2eaeb8337b60b24b126a9c049d6315f)
  * provide argument to specify number of populations per legend row [c7fe85d](https://github.com/meyer-lab/plinkQC/commit/c7fe85d555c256bf8be3843803ff1cc445f7898f)
  * Ensure correct ordering of population colors when reference population is not HapMap [86275b7](https://github.com/meyer-lab/plinkQC/commit/86275b7c281835681f38ab16672db61c288cfed9)

# plinkQC 0.3.1
## minor changes

* Fixed dead links in vignettes (caused by migration of repository): [da987d8](https://github.com/meyer-lab/plinkQC/commit/da987d8f225aa6aca0596b9c4f6a2484b102bdb6)
* Added note about chrY in Hapmap data (vignette): [e8afbb9](https://github.com/meyer-lab/plinkQC/commit/e8afbb9842ed9421461a8114ac0a00f7955cf0c0)
* Added note about recommended use of plink1.9 (vs 2.0):
[b69d3d7](https://github.com/meyer-lab/plinkQC/commit/b69d3d71d23e9b161176a635fcb2b5a2b524591f)

## Test environments
* local R installation, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

## Comments
* This is a minor release
* Mostly updates to vignettes and documentation
* Some additional checks and tests on ancestry check
