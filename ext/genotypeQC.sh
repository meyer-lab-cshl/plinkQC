###############################################################################
###                                                                         ###
### Pipeline (partly) based on "Anderson et al. (2010) Data quality         ###
### control in genetic case-control association studies."                   ###
### Nature protocols 5(9):1564-73  	                                        ###
###                                                                         ###
### Script deals mainly with per-individual QC and passes on parameters for ###
### per-marker qc to the Rscript called at the end                          ###
###                                                                         ###
###                                                                         ###
###     by Hannah Meyer                                                     ###
###                                                                         ###
###############################################################################

####################################
### 1. input parameters and data ###
####################################

### a) parameters
# set default parameters
ethnicityFilter=yes
hapmap=HapMapIII_CGRCh37.good_pos
ensemblNotInUK10K1KG=~/GWAS/data/hmeyer/Supplementary/20160415_GRCh37p13_not_in_UK10K1KG.txt

# source config file
. $1

# change working directory to data dir
cd $qcdir
mkdir -p plink_log

### b)  format raw data and transfer to QC directory ###
if [[ ! -s $qcdir/${alg}.raw.fam ]]; then
    ### Reformat input files if parental information is encoded incorrectly:
    if [[ $center == *"singapore"* ]]; then
        reformat=1
        echo "Reformating Singapore files..."
        perl ~/GWAS/analysis/genotyping/parse_singapore_data.pl \
            --inped $rawdir/$rawdata \
            --samples $rawdir/${rawdata}_sample_sex.txt \
            --outped $rawdir/$alg.raw \
            --incorrect $rawdir/corrected_sample_id.txt
    elif [[ $center == *"sanger"* ]]; then
        if  cut -d " " -f 3 $rawdir/${rawdata}.fam | uniq |  grep -q -e "-9" ; 
        then
            reformat=1
            echo "Reformating Sanger files..."
            # get original non-binary file
            echo "Convert into ped/map file..."
            plink  --bfile $rawdir/$rawdata  --recode --out $rawdir/${alg}.orig

            # PLINK doesn't accept -9 as parental information
            # Rename the samples: rm prefix "urn:wtsi:[\d\w]_A01_"
            echo "Format sample identifiers"
            perl ~/GWAS/analysis/general/Replace.pl $rawdir/$alg.orig.ped  \
                -sub_string "-9 -9" \
                -replace_string "0 0" | perl ~/GWAS/analysis/general/Replace.pl \
                -sub_string "urn:wtsi:[\d_ABCDEFGH]*_" \
                -replace_string "" \
                -out_file $rawdir/$alg.tmp.ped
            mv $rawdir/$alg.orig.map $rawdir/$alg.tmp.map

            # make .bed again
            echo "Convert into binary (bed/bim/fam) format again"
            plink --file $rawdir/${alg}.tmp \
                --make-bed \
                --out $rawdir/${alg}.raw
        else
            reformat=0
            cp $rawdir/$rawdata.bim $qcdir/$alg.raw.bim
            cp $rawdir/$rawdata.fam $qcdir/$alg.raw.fam
            cp $rawdir/$rawdata.bed $qcdir/$alg.raw.bed
        fi
    fi

    if [[ $reformat == 1 ]]; then
        ### Copy reformated data from rawdir to qcdir
        cp $rawdir/$alg.raw.bim $qcdir/$alg.raw.bim
        cp $rawdir/$alg.raw.fam $qcdir/$alg.raw.fam
        cp $rawdir/$alg.raw.bed $qcdir/$alg.raw.bed
    fi
fi

### c) Check for duplicate IDs in dataset ###
cut -d " " -f 1 $qcdir/$alg.raw.fam | sort | uniq -d \
    > $qcdir/$alg.duplicate.IDs
if [[ -s $qcdir/$alg.duplicate.IDs ]]; then
    echo "Remove duplicate ID's..."
    plink --bfile $qcdir/$alg.raw --recode --out $qcdir/$alg.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.recode.log
    egrep -v -f $qcdir/$alg.duplicate.IDs  $alg.raw.ped \
        > $alg.raw.duplicate.remove.ped
    plink --ped $qcdir/$alg.raw.duplicate.remove.ped \
        --map $qcdir/$alg.raw.map \
        --make-bed --out $qcdir/${alg}.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.duplicateremove.log
fi

### d)  Remove IDs (files generated via ~/GWAS/analysis/GWAS/IDmatching.R)
if [[ -s $qcdir/$alg.preQCfail.IDs ]]; then
    plink --noweb --bfile $qcdir/$alg.raw \
        --remove $qcdir/$alg.preQCfail.IDs \
        --make-bed --out $qcdir/$alg.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.removequestionableIDs.log
fi

##################################################
### 2. match input files to reference datasets ###
##################################################

# reference data sets: imputation (UK10K1KG) and ENSEMBL 
# (if not present in UK10K1000Genomes) ###
# (ENSEMBL generated via ~/genotyping/format_20160415_GRCh37p13_human_vcf.sh)

### a) remove variants that are neither ensembl nor UK10K1KG annotated 
# (important for later fusion of genotypes from different batches)
awk 'FNR==NR {a[$2]++; next} a[$2] ' $qcdir/$alg.raw.bim  $UK10K1KGdir/genome.legend.format > $qcdir/${alg}_genome.legend.format
awk 'FNR==NR {a[$2]++; next} a[$3] ' $qcdir/$alg.raw.bim  $ensemblNotInUK10K1KG > $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt

awk 'FNR==NR {a[$2]++; next} !a[$2] ' $UK10K1KGdir/genome.legend.format $qcdir/$alg.raw.bim  > $qcdir/${alg}_not_in_genome.legend.format
awk 'FNR==NR {a[$3]++; next} !a[$2] ' $ensemblNotInUK10K1KG $qcdir/${alg}_not_in_genome.legend.format  > $qcdir/${alg}_not_in_GRCh37p13_or_UK10K1KG

plink --noweb --bfile $qcdir/$alg.raw \
    --exclude <(cut -f 2 $qcdir/${alg}_not_in_GRCh37p13_or_UK10K1KG) \
    --make-bed --out $qcdir/$alg.onlyannotated
mv $qcdir/${alg}.onlyannotated.log $qcdir/plink_log

### b) check if SNPs are concordant with UK10K1KG information
# get summary of SNPs in batch and UK10K1KG
plink1.07 --noweb --id-match $qcdir/$alg.onlyannotated.bim SNP,2 \
    $qcdir/${alg}_genome.legend.format SNP,2 + complete \
    --out $qcdir/${alg}.match_UK10K1KG
mv $qcdir/${alg}.match_UK10K1KG.log $qcdir/plink_log

# check if chromosomes match
awk '$1 != $7 {print $2, "\t", $7}' $qcdir/${alg}.match_UK10K1KG.matched \
    > $qcdir/${alg}.match_UK10K1KG.bad_chr
plink --bfile  $qcdir/$alg.onlyannotated \
    --update-map $qcdir/${alg}.match_UK10K1KG.bad_chr --update-chr \
    --make-bed --out $qcdir/$alg.UK10K1KG.good_chr
mv $qcdir/${alg}.UK10K1KG.good_chr.log $qcdir/plink_log

# check if chromosomal position match
awk '$4 != $10 {print $2, $10}' $qcdir/${alg}.match_UK10K1KG.matched > \
    $qcdir/${alg}.match_UK10K1KG.bad_pos

# check if alleles match
awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' $qcdir/${alg}.match_UK10K1KG.matched  > $qcdir/${alg}.match_UK10K1KG.bad_allele ## FLIP or EXCLUDE
plink --bfile $qcdir/${alg}.UK10K1KG.good_chr \
    --update-map $qcdir/${alg}.match_UK10K1KG.bad_pos \
    --flip $qcdir/${alg}.match_UK10K1KG.bad_allele \
    --make-bed \
    --out $qcdir/$alg.UK10K1KG.flipped
mv $qcdir/${alg}.UK10K1KG.flipped.log $qcdir/plink_log

# check if re-mapping was successful
plink1.07 --noweb --id-match $qcdir/$alg.UK10K1KG.flipped.bim SNP,2 \
    $qcdir/${alg}_genome.legend.format SNP,2 + complete \
    --out $qcdir/$alg.flipped_UK10K1KG
mv $qcdir/${alg}.flipped_UK10K1KG.log $qcdir/plink_log
awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' $qcdir/$alg.flipped_UK10K1KG.matched > $qcdir/$alg.flipped.UK10K1KG.bad_allele

# generate UK10K1KG-matched plink file set 
plink --bfile $qcdir/$alg.UK10K1KG.flipped \
    --exclude $qcdir/$alg.flipped.UK10K1KG.bad_allele \
    --make-bed \
    --out $qcdir/$alg.UK10K1KG
mv $qcdir/${alg}.UK10K1KG.log $qcdir/plink_log/${alg}.UK10K1KG.log

### c) check if SNPs are concordant with ENSEMBL reference info (for SNPs not in UK10K1KG reference set) 
plink1.07 --noweb --id-match $qcdir/$alg.UK10K1KG.bim SNP,2 \
    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt SNP,3 + complete \
    --out $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG
mv $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.log $qcdir/plink_log

# check if chromosomes match
awk '$1 != $7 {print $2, "\t", $7}' $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.matched > $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_chr
plink --bfile  $qcdir/$alg.UK10K1KG \
    --update-map $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_chr \
    --update-chr \
    --make-bed \
    --out $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.good_chr
mv $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.good_chr.log $qcdir/plink_log

# check if chromosomal position match
awk '$4 != $8 {print $2, $8}' $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.matched > $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_pos

# check if alleles match
awk '($5$6 != $10$11) && ($5$6 != $11$10) {print $2}' $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.matched  > $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_allele ## FLIP or EXCLUDE
plink --bfile $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.good_chr \
    --update-map $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_pos \
    --flip $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.bad_allele \
    --make-bed --out $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped
mv $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.log $qcdir/plink_log

# check if re-mapping was successful
plink1.07 --noweb \
    --id-match $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.bim SNP,2 \
    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt SNP,3 + complete \
    --out $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped
mv $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.log $qcdir/plink_log
awk '($5$6 != $10$11) && ($5$6 != $11$10) {print $2}' $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.matched > $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.bad_allele

# generate UK10K1KG and ENSEMBL-matched plink file set for QC
plink --bfile $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped \
    --exclude $qcdir/${alg}.GRCh37.p13_not_in_UK10K1KG.flipped.bad_allele \
    --make-bed \
    --out $qcdir/$alg
mv $qcdir/${alg}.log $qcdir/plink_log/${alg}.UK10K1KG_ensembl.annotated.log


##########################################################################
### 3. Quality Control (# according to Anderson et al (2010) protocol) ###
##########################################################################

### a) Identification of Individuals with discordant sex information
#4  calculate mean homozygosity rate across X for each individual
plink --bfile $qcdir/$alg --check-sex --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.sexcheck
grep 'PROBLEM' $qcdir/$alg.sexcheck > $qcdir/$alg.failsex

### b) Identification of individuals with elevated missing data rates or 
###outlying heterozygosity rate
#7  creates files .imiss (4th column: number of missing SNPs N_MISS, 
# 6th column: missing SNPs %/individum F_miss) and .lmiss
plink --bfile  $qcdir/$alg --missing --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.missing

#8  creates the file .het (3rd column: number of homozygous gt O(Hom), 
# 5th column: number of nonmissing genotypes per individual N(NM)
plink --bfile  $qcdir/$alg --het --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.hr

### c) Identification of duplicated or related individuals
#11  minimize computational complexity by reducing the number of SNPs in 
# IBS matrix through removal of SNPs within LD (r2 threshold 0.05) and MAF 
# threshold results in files: .prune.in --> SNPs kept in the analyses; 
plink --bfile $qcdir/$alg --exclude range $highld --indep-pairwise 50 5 0.2 \
    --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.prune
#12  create pairwise IBS for all pairs of individuals
plink --bfile  $qcdir/$alg --extract $qcdir/$alg.prune.in --make-bed \
    --out $qcdir/$alg.pruned
mv $qcdir/$alg.pruned.log $qcdir/plink_log/$alg.pruned
plink --bfile $qcdir/$alg.pruned --maf 0.1 --genome --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.pruned.IBS

### d) Identification of individuals of divergent ancestry 
# 15 merge study genotypes to HapMap3 data
# Filter hapmap3 data for the same SNPset as in study 
plink --noweb --bfile  $hapmapdir/$hapmap --extract $qcdir/$alg.prune.in \
    --make-bed --out $qcdir/$hapmap.pruned
mv  $qcdir/$hapmap.pruned.log $qcdir/plink_log/$hapmap.pruned

# check if SNPs which are in the HapMapIII (position matched to UK10K1KG) 
# reference sets are concordant with HapMapIII info by comparing SNP names·
plink1.07 --noweb --id-match $qcdir/$hapmap.pruned.bim SNP,2  \
    $qcdir/$alg.pruned.bim SNP,2 + complete --out $qcdir/$hapmap.pruned
mv $qcdir/$hapmap.pruned.log $qcdir/plink_log

# check if chromosomes match
awk '$1 != $7 {print $2, "\t", $7}' $qcdir/$hapmap.pruned.matched > $qcdir/$hapmap.pruned.bad_chr
plink --bfile  $qcdir/$hapmap.pruned \
    --update-chr $qcdir/$hapmap.pruned.bad_chr --make-bed \
    --out $qcdir/$hapmap.pruned.good_chr
mv $qcdir/$hapmap.pruned.good_chr.log $qcdir/plink_log

# check if chromosomal position match
awk '$4 != $10 {print $2, $10}' $qcdir/$hapmap.pruned.matched > \
    $qcdir/$hapmap.pruned.bad_pos

# check if alleles match
awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' $qcdir/$hapmap.pruned.matched  > $qcdir/$hapmap.pruned.bad_allele ## FLIP or EXCLUDE

# remove variants with bad positions and flip alleles that didn't match
plink --bfile $qcdir/$hapmap.pruned.good_chr  \
    --exclude $qcdir/$hapmap.pruned.bad_pos \
    --flip $qcdir/$hapmap.pruned.bad_allele \
    --make-bed --out $qcdir/$hapmap.pruned.flipped
mv $qcdir/$hapmap.pruned.flipped.log $qcdir/plink_log

#check if re-mapping was successful
plink1.07 --noweb --id-match $qcdir/$hapmap.pruned.flipped.bim SNP,2 \
    $qcdir/$alg.pruned.bim SNP,2 + complete --out $qcdir/$hapmap.pruned.flipped
mv $qcdir/$hapmap.pruned.flipped.log $qcdir/plink_log
awk '($5$6 != $11$12) && ($5$6 != $12$11) {print $2}' $qcdir/$hapmap.pruned.flipped.matched > $qcdir/$hapmap.pruned.flipped.bad_allele
plink --bfile $qcdir/$hapmap.pruned.flipped \
    --exclude $qcdir/$hapmap.pruned.flipped.bad_allele \
    --make-bed --out $qcdir/$hapmap.clean
mv $qcdir/$hapmap.clean.log $qcdir/plink_log/$hapmap.clean.log

# merge HapMap-mapped genotypes to genotypes of HapMap populations
plink --bfile $qcdir/$alg.pruned  \
    --bmerge $qcdir/$hapmap.clean.bed $qcdir/$hapmap.clean.bim \
    $qcdir/$hapmap.clean.fam  \
    --make-bed --out $qcdir/$alg.merge.$hapmap
mv $qcdir/$alg.merge.$hapmap.log $qcdir/plink_log

#17 Conduct PCA on the merged data
plink --bfile $qcdir/$alg.merge.$hapmap --pca \
    --out $qcdir/$alg.HapMapIII.pruned.pca
mv $qcdir/$alg.HapMapIII.pruned.pca.log $qcdir/plink_log

##############################################################################
### 7.  QC for per-individual and per-marker analysis: Plots and filtering ###
##############################################################################


Rscript ~/GWAS/analysis/genotyping/genotypeQC.R --vanilla \
    --default-packages=R.utils \
    --alg=$alg \
    --qcdir=$qcdir \
    --sexCheckfilter=$sexCheckfilter \
    --maleTh=$maleTh \
    --femaleTh=$femaleTh \
    --imissTh=$imissTh \
    --hetTh=$hetTh  \
    --C1=$C1 \
    --C2=$C2 \
    --C3=$C3 \
    --highIBDfilter=$highIBDfilter \
    --highIBDTh=$highIBDTh \
    --crypticIBDTh=$crypticIBDTh \
    --lmissTh=$lmissTh \
    --hweTh=$hweTh \
    --mafTh=$mafTh \
    --omnixfile=$omnix \
    --ethnicityFilter=$ethnicityFilter \
    --center=$center \
    --sample=$sample

