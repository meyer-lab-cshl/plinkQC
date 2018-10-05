###############################################################################
###                                                                         ###
### Genotype Quality Control                                                ###
###                                                                         ###
###     * numbers in QC indicate steps described in the protocol            ###
###      "Anderson et al. (2010) Data quality control in genetic            ###
###       case-control association studies." Nature protocols               ###
###       5(9):1564-73                                                      ###
###                                                                         ###
###                                                                         ###
###     by Hannah Meyer                                                     ###
###                                                                         ###
###############################################################################

####################################
### 1. input parameters and data ###
####################################

### a) parameters
alg=$1
qcdir=$2
highld=$3

mkdir -p $qcdir/plink_log

##############
### 2. QC  ###
##############

### a) Identification of Individuals with discordant sex information
# 4.  calculate mean homozygosity rate across X for each individual
plink --bfile $qcdir/$alg --check-sex --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.sexcheck
grep 'PROBLEM' $qcdir/$alg.sexcheck > $qcdir/$alg.failsex

### b) Identification of individuals with elevated missing data rates or
### outlying heterozygosity rate
# 7.  creates files .imiss (4th column: number of missing SNPs N_MISS,
# 6th column: missing SNPs %/individum F_miss) and .lmiss
plink --bfile  $qcdir/$alg --missing --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.missing

# 8.  creates the file .het (3rd column: number of homozygous gt O(Hom),
# 5th column: number of nonmissing genotypes per individual N(NM)
plink --bfile  $qcdir/$alg --het --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.hr

### c) Identification of duplicated or related individuals
# 11.  minimize computational complexity by reducing the number of SNPs in
# IBS matrix through removal of SNPs within LD (r2 threshold 0.05) and MAF
# threshold results in files: .prune.in --> SNPs kept in the analyses;
plink --bfile $qcdir/$alg --exclude range $highld --indep-pairwise 50 5 0.2 \
      --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.prune

# 12.  create pairwise IBS for all pairs of individuals
plink --bfile  $qcdir/$alg --extract $qcdir/$alg.prune.in --make-bed \
      --out $qcdir/$alg.pruned
mv $qcdir/$alg.pruned.log $qcdir/plink_log/$alg.pruned

plink --bfile $qcdir/$alg.pruned --maf 0.1 --genome --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.pruned.IBS
