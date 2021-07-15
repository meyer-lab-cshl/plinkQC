# plinkQC 0.3.4
## major changes
* Add flag for removing high LD regions and add option for user provided file,
allowing removal of high LD regions for non-human organisms; fix missing flag
[19bc569](https://github.com/meyer-lab-cshl/plinkQC/commit/19bc5697f80c4766d91c24d773337b8bf88a3e6f) and
[#42](https://github.com/meyer-lab-cshl/plinkQC/issues/42)
* Fixed failing data.no_failIDs.lmiss file in check_snp_missingness function
when failIDs were present
[86adc39](https://github.com/meyer-lab-cshl/plinkQC/commit/86adc39764ed6271e947d12063229be0e0189d12) and
[#43](https://github.com/meyer-lab-cshl/plinkQC/issues/43)

# plinkQC 0.3.3
## major changes
* Output all IBD estimates in `run_check_relatedness` not just fail IBD
estimates [d009da90](https://github.com/meyer-lab-cshl/plinkQC/commit/d009da90434a83cb4867881e3772b518d17dc2dd)
* Add midp modifier for --hwe and --hardy as recommended on [plink hwe documentation](https://www.cog-genomics.org/plink/1.9/filter#hwe) and discussed in issue [#27](https://github.com/meyer-lab-cshl/plinkQC/issues/27); [43230a6c](https://github.com/meyer-lab-cshl/plinkQC/commit/43230a6ce5f079b46c85b1a7a526213177db1f1b)
* Provide HapMap and 1000Genomes individual and population identifiers as default
for ancestry check functions; include optional color mapping for these reference
populations; raised in [#37](https://github.com/meyer-lab-cshl/plinkQC/issues/37), incorporated in [65cfcba](https://github.com/meyer-lab-cshl/plinkQC/commit/65cfcbaf5b0af3b1d3d4b710a8710331a39e0f6e)
* for all individual QC functions return data.frame with FID, IID and measure
used for filtering [#5](https://github.com/meyer-lab-cshl/plinkQC/issues/5)
incorporated in [295505b](https://github.com/meyer-lab-cshl/plinkQC/commit/295505bc2d2a6278539e1d078c338bcf3a128338)
* allow for user-supplied individual and marker filtering in perIndividualQC,
perMarkerQC and cleanData [#33](https://github.com/meyer-lab-cshl/plinkQC/issues/33)
* individual highlighting in perIndividualQC plots; several highlight modalities
(text/label, shape, color) [#38](https://github.com/meyer-lab-cshl/plinkQC/issues/38)


## minor changes
* Add checks in `cleanData` that can distinguish between miss-specified filters and non-failing samples [7432b1ee](https://github.com/meyer-lab-cshl/plinkQC/commit/7432b1ee880ad14c7a078f3a0379f7fdcff1684b); addresses issue [#31](https://github.com/meyer-lab-cshl/plinkQC/issues/31) 
* Add checks and tests in `check_snp_missingness`, `check_maf` and `check_hwe` 
in cases where all samples pass `perIndividualQC`.
[31ed6f65](https://github.com/meyer-lab-cshl/plinkQC/commit/31ed6f653902857a0a78ae723e786dd84d240b75); related to issue [#31](https://github.com/meyer-lab-cshl/plinkQC/issues/31) 
* Add checks and tests in `evaluate_check_sex` and `evaluate_check_het_imiss`
for all passing samples and accordingly, no fail samples label in plot.
[d071a43](https://github.com/meyer-lab-cshl/plinkQC/commit/d071a4341b0d4c40220b42f81c412da4d320b8ec)
* Update dropbox links and add note to vignette [9bae988](https://github.com/meyer-lab-cshl/plinkQC/commit/9bae9886dc637569698d883654fcf71486b6634e)
* additional plotting parameters for `evaluate_*` functions, controlling 
text and title sizes ([b38a348](https://github.com/meyer-lab-cshl/plinkQC/commit/b38a348da0a4dd60860543b2d53c4de1df058a54),
[f25b8d0](https://github.com/meyer-lab-cshl/plinkQC/commit/f25b8d05bc3462228afb26bb1caec8d829bf7fca))


# plinkQC 0.3.2
## minor changes
* Add checks and tests in `evaluate_check_ancestry` for missing non-reference samples
[5c03971](https://github.com/meyer-lab-cshl/plinkQC/commit/5c03971b92341b5a8cb26a1d5f2e37d637f5cbd7)
* Clarified tutorials:
  * R internal file.copy instead of systems `cp` command to make run on windows: [6c91ef3](https://github.com/meyer-lab-cshl/plinkQC/commit/6c91ef3368f601658967f88bc20db94a28307fce),
  fixes [#26](https://github.com/meyer-lab-cshl/plinkQC/issues/26) 
  * Add additional filter for A->T, C->G variants as suggested in [#24](https://github.com/meyer-lab-cshl/plinkQC/issues/24); [11e0375](https://github.com/meyer-lab-cshl/plinkQC/commit/11e03756852215915310a836d2e73dbcedf12bfc)
* Fixed plotting issues in PCA plot of ancestry check:
  * allow for supplying names of European reference population [ae09e64](https://github.com/meyer-lab-cshl/plinkQC/commit/ae09e64ee2eaeb8337b60b24b126a9c049d6315f)
  * provide argument to specify number of populations per legend row [c7fe85d](https://github.com/meyer-lab-cshl/plinkQC/commit/c7fe85d555c256bf8be3843803ff1cc445f7898f)
  * Ensure correct ordering of population colors when reference population is not HapMap [86275b7](https://github.com/meyer-lab-cshl/plinkQC/commit/86275b7c281835681f38ab16672db61c288cfed9)

# plinkQC 0.3.1
## minor changes

* Fixed dead links in vignettes (caused by migration of repository): [da987d8](https://github.com/meyer-lab-cshl/plinkQC/commit/da987d8f225aa6aca0596b9c4f6a2484b102bdb6)
* Added note about chrY in Hapmap data (vignette): [e8afbb9](https://github.com/meyer-lab-cshl/plinkQC/commit/e8afbb9842ed9421461a8114ac0a00f7955cf0c0)
* Added note about recommended use of plink1.9 (vs 2.0):
[b69d3d7](https://github.com/meyer-lab-cshl/plinkQC/commit/b69d3d71d23e9b161176a635fcb2b5a2b524591f)

# plinkQC 0.3.0
## major changes
* Relationship filter can handle more complicated relationship scenarios as
  observed in plant genotype sets (fixed [#11](https://github.com/meyer-lab-cshl/plinkQC/issues/11))
* code moved to the meyer-lab repository: https://github.com/meyer-lab-cshl/plinkQC

# plinkQC 0.2.3
## major changes
* Enable return of overview plots as ggplot object
  (fixes  [#6](https://github.com/meyer-lab-cshl/plinkQC/issues/6) in [ab2e840](https://github.com/meyer-lab-cshl/plinkQC/commit/ab2e840f0f22ccdeb5317475698aa58e7eecd345) and
  [101e74e](https://github.com/meyer-lab-cshl/plinkQC/commit/101e74e318752b2f8038330245d875f6afd5b57b))
* Relationship filter now deals with more complicated relationship scenarios as
  observed in plant genotype sets (addresses [#11](https://github.com/meyer-lab-cshl/plinkQC/issues/11))

## minor changes
* give user option to choose maf threshold for relatedness filtering (relates to
   [#3](https://github.com/meyer-lab-cshl/plinkQC/issues/3))

## bug fixes
* Include check for zero related individuals fixing  [#3](https://github.com/meyer-lab-cshl/plinkQC/issues/3) in
  [1445a88](https://github.com/meyer-lab-cshl/plinkQC/commit/1445a88d41d985e73ffdd161144229683f8352cd)
* Include check in case all samples fail perIndividual QC in
   [894acc1](https://github.com/meyer-lab-cshl/plinkQC/commit/894acc1fa03dadfe0ad2028888142171bcc641eb) and
   [0464224](https://github.com/meyer-lab-cshl/plinkQC/commit/04642246d18ed4eaac5b9d9a6931d1ecb08308e8))
* Include checks for diagonal derived relationship estimates, and estimate data
  containing only related individuals; addresses  [#11](https://github.com/meyer-lab-cshl/plinkQC/issues/11)
* Fix command for genotype conversion in 1000Genomes vignette, addressing issue
   [#10](https://github.com/meyer-lab-cshl/plinkQC/issues/10)
* fix missing rownames error for overviewPerIndividualQC, when relatedness check
  was included (issue [#16](https://github.com/meyer-lab-cshl/plinkQC/issues/16), [fc7a38b](https://github.com/meyer-lab-cshl/plinkQC/commit/fc7a38b1f2b345d9c6c5d69f5dcf0bc57a857f62))
* fix vignette mismatch (issue [#16](https://github.com/meyer-lab-cshl/plinkQC/issues/16), [09dcd59](https://github.com/meyer-lab-cshl/plinkQC/commit/09dcd59e77178b35747aae81a5c1988712e20de9))

# plinkQC 0.2.2
## minor changes
* Fix IDs written to fail.IDs file: previous versions wrote IID,IID, now fixed
to FID, IID (fixes [#2](https://github.com/meyer-lab-cshl/plinkQC/issues/2)).

# plinkQC 0.2.1
## minor changes
* Fix path check bug in checkPlink
* Include test data in build!

# plinkQC 0.2.0
## major changes
* All system calls to plink are conducted with sys::exec_wait - this should
solve platform dependent issues, mainly compatibility with windows.
* Make path construction compatible with windows.
* path2plink now requires full path to plink executable, no tilde expansion or
simple pointer to directory supported.
* Fix bug in return of cleanData function: list now contains keep and fail IDs.
* Fix bug in return of maf computation: if fail.IDs does not exist, set
  fail_samples to zero.

## minor changes
* IBD-fail.IDs now saved without column names to be consistent with other
xxx-fail.IDs files.
* Include additional progress messages in cleanData()
* Remove default double-specification of mafTh and macTh
* use checkPlink to return correct path2plink and export to make checkPlink
directly accesible to user.


# plinkQC 0.1.1
## major changes
* run_check_relatedness will only save IBD estimates of individuals whose
estimates are higher than the threshold.

## minor changes
* Fix examples in vignettes 1000 Genomes and HapMap III reference.
* Change file access in function examples
* Add additional checks in check_ancestry and fix missing refSamplesFile test


# plinkQC 0.1.0

* Added a `NEWS.md` file to track changes to the package.



