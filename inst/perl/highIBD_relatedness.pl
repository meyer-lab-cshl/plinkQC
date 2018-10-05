#! /nfs/software/software/prefix/usr/bin/perl

#######################################################
#######################################################
###                                                 ###
###     Script to filter related samples based on   ### 
###         * PI_HAT (plink-generated IBD meassure) ###
###         * SNP missingness rates                 ###
###         * known family structure                ###
###         * patient status                        ###
###                                                 ###
###     by Hannah Meyer                             ###
###                                                 ###
#######################################################
#######################################################

################
### modules  ###
################

use strict;
use diagnostics;
use warnings;

use Getopt::Long;
use List::Util qw(reduce sum min);

#################
### variables ###
#################

# defaults 
my $genoID_field = 1;
my $famID_field = 3;
my $diagnosisID_field = 6;
my $missThr = 0.03;

my ($file, $thres, $famfile, $diag);
my (%imiss, %removed,%pihat, %fam2id, %id2fam, %fail_id, %diagnosis);


############################
### example command line ###
############################
 # /homes/hannah/GWAS/analysis/genotyping/highIBD.pl --file /homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall_all/HVOL/gencall_all --thres 0.125 --famfile /homes/hannah/GWAS/data/phenotype/2Dphenotype/20160209_All_BRU_family_format.txt

###################
### subroutines ###
###################

sub usage {
    print "Unknown option: @_\n" if (@_);
    print "usage: Script to filter related samples based on PI_HAT (plink-generated IBD meassure), SNP missingness rates, known family structure and patient status";
    print "\n
           --file\tPath and prefix to plink-file fow which relatedness should be analysed
           --thres\tHigh IBD threshold (standard: 0.125 for 3rd degree relatives)
           --famfile\tPath to file containing known-family relations between samples \n";
    exit;
}


### returns index of array element with minimmum value
sub argmin {
  	my ($ar) = @_;
   	return reduce { $ar->[$b] < $ar->[$a] ? $b : $a } (0..$#$ar);
}

#reduce {$arr[$a] > $arr[$b] ? $a : $b} (0 .. $#arr);
### returns index of array element with maximum value
sub argmax {
  	my ($ar) = @_;
   	return reduce { $ar->[$b] > $ar->[$a] ? $b : $a } (0..$#$ar);
}

################
### analyses ###
################

##########################################
### 1. input parameters and input data ###
##########################################

### a) Read in command line options
usage() if (@ARGV <  1 or ! GetOptions ("file=s" => \$file, "thres=s" => \$thres, "famfile=s" => \$famfile));#, "diag=s" => \$diag));
GetOptions ("file=s" => \$file, "thres=s" => \$thres, "famfile=s" => \$famfile);# , "diag=s" => \$diag);

### b) read family and patient status info
print "Reading family ID file\n";
open (FAM, '<', $famfile) or die "Cannot open family ID file (".$famfile."): $!\n";
while(defined( my $line = <FAM>)){
    # skip header
	next if $. == 1;
	chomp $line;

    # find individual and family ids (famID = NA if doesnt exist)
    my @array = split "\t", $line; 
	my $famID = $array[$famID_field -1];
	my $genoID = $array[$genoID_field -1];
    $diagnosis{$genoID} = $array[$diagnosisID_field -1];
	
    # create entry in lookup tables if individual part of family :
    # - mapping family Ids to all individaul ids (fam2id)
    # - mapping individual Ids to their family ids (id2fam)
    if ($famID  ne "NA") {
		push @{$fam2id{$famID}}, $genoID;
		$id2fam{$genoID} = $famID;
	}
}
close(FAM);


### c)  read SNP missingness rates per individual (plink output-file of 'plink --bfile --missing --out' with missing data rate per individual)
print "Reading PLINK .imiss file (missing genotyping rate) ".$file.".imiss\n";
open (IMISS, '<', $file.".imiss") or die "Cannot open genotypes file (".$file.".imiss): $!\n";
while(defined( my $line = <IMISS>)){
    # skip header
	next if $. == 1;
    chomp $line;
    
	# remove leading white spaces of first column
    $line =~ s/^\s+//;
    my @array = split /\s+/, $line;
    
    # individual id as hash key and save genotype missing rate as value
    $imiss{$array[0]} = $array[5];
}
close(IMISS);

###########################
### 2. Sample filtering ###
###########################

### a) read pairwise PIHAT measurement as indicator for relatedness (plink output-file of 'plink --bfile --genome --out')
###  i) filtering based on patient status
print "Reading PLINK .genome file (PIHAT measurement) ".$file.".genome\n";
open (GENOME, '<', $file.".genome")  or die "Cannot open genotypes file (".$file.".genome): $!\n";
open (OUT, '>', $file.".fail-IBD.IDs") or die "Cannot write to fail-ID file (".$file.".fail-IBD.IDs): $!\n";        

open (OUT_related, '>', $file.".fail-IBD.relatedIDs") or die "Cannot write to fail-ID file (".$file.".fail-IBD.IDs): $!\n";   
print OUT_related "Excluded ID\tRelated ID\tPIHAT\tExcluded Family ID\tRelated Family ID\tExcluded imiss\tRelated imiss\tExcluded diagnosis\tRealted diagnosis\n";

while(defined(my $line =<GENOME>)){
    # skip header
	next if $. == 1;
    chomp $line;
    
	# remove leading white spaces of first column
    $line =~ s/^\s+//;
    my @array = split /\s+/, $line;
    
    # check if PI_HAT value greater than threshold
    my $both=0;
    if ($array[9] >= $thres) {

        # get famID if individuals are part of registered family
        my ($famID1, $famID2) = ("", "");
        $famID1 = $id2fam{$array[0]} if exists $id2fam{$array[0]};
        $famID2 = $id2fam{$array[2]} if exists $id2fam{$array[2]};
        print $array[0], "\t", $array[2], "\t", $imiss{$array[0]}, "\t", $imiss{$array[2]}, "\t", $famID1, "\t",  $famID2, "\t", $diagnosis{$array[0]}, "\t", $diagnosis{$array[2]},"\n";

		# if either of the individuals failed the missingness threshold, don't consider pair for IBD analysis 
		# create key in fail_id 
		if ( ($imiss{$array[0]} >= $missThr) &&  ($imiss{$array[2]} >= $missThr) ) {
            # only write one as failed IBD to be consistent with rest of analyses (only filter one sample per pair)
            print OUT "$array[0] $array[0]\n" unless exists $fail_id{$array[0]};
            print OUT_related join("\t", ($array[0], $array[2], $array[9], $famID1, $famID2, $imiss{$array[0]}, $imiss{$array[2]},  $diagnosis{$array[0]},  $diagnosis{$array[2]})), "\n" unless exists $fail_id{$array[0]};
            $fail_id{$array[0]} = 1;
            $fail_id{$array[2]} = 1;
            $both = 1;
        # if disease and healthy cohorts are combined (eg HVOL_DCM) remove related control sample
        } elsif ($diagnosis{$array[0]} ne $diagnosis{$array[2]}) {
            my $fail_diag = $diagnosis{$array[0]} eq "HVOL" ? $array[0] : $array[2];
            my $keep_diag = $diagnosis{$array[0]} ne "HVOL" ? $array[0] : $array[2];
            print OUT "$fail_diag $fail_diag\n" unless exists $fail_id{$fail_diag};
            print OUT_related join("\t", ($fail_diag, $keep_diag, $array[9], $famID1, $famID2, $imiss{$array[0]}, $imiss{$array[2]},  $diagnosis{$fail_diag},  $diagnosis{$keep_diag})), "\n" unless exists $fail_id{$fail_diag};
            $fail_id{$fail_diag} = 1;
        } elsif (($imiss{$array[0]} >= $missThr) && ($both == 0)) {
            print OUT "$array[0] $array[0]\n" unless exists $fail_id{$array[0]};
            print OUT_related join("\t", ($array[0], $array[2], $array[9], $famID1, $famID2, $imiss{$array[0]}, $imiss{$array[2]},  $diagnosis{$array[0]},  $diagnosis{$array[2]})), "\n" unless exists $fail_id{$array[0]};
            $fail_id{$array[0]} = 1;
        } elsif (($imiss{$array[2]} >= $missThr)  && ($both == 0)) {
            print OUT "$array[2] $array[2]\n" unless exists $fail_id{$array[2]};
            print OUT_related join("\t", ($array[2], $array[0], $array[9], $famID1, $famID2, $imiss{$array[2]},$imiss{$array[0]},  $diagnosis{$array[2]},  $diagnosis{$array[0]})), "\n" unless exists $fail_id{$array[2]};
            $fail_id{$array[2]} = 1;
        } else {
			# if both individuals pass missingness threshold evaluate which one to keep in analysis
            #print $array[0], "\t", $array[2], "\n";
			$pihat{$array[0]}{$array[2]} = $array[9];
		}
	}
}
close(GENOME);


### b) iterate over all known family relations and filter such that maximal number of individuals per family is retained in study
foreach my $keys (keys %fam2id) {
	my %pi =();
    my %pi_id= ();

    # iterate over all possible individual combinations of the family 
	foreach my $element (@{$fam2id{$keys}}) {
		foreach my $e (@{$fam2id{$keys}}) {
			next if $element eq $e;

            # if combination exists in %pihat, ie IBD > threshold, save Id names (both directions) and PIHAT value
			if (exists $pihat{$element}{$e}) {
				$pi{$element.$e} = $pihat{$element}{$e};
				$pi{$e.$element} = $pihat{$element}{$e};
                $pi_id{$element}++;
                $pi_id{$e}++;
			}
		}
	}
    # print scalar keys %pi_id, "\n";
    # if there are more than 2 individuals of the same family with IBD > threshold, select which one to exclude from future analysis 
    if (scalar keys %pi_id > 2 ) {
        
        # so far only families of three can be easily selected, give warning if larger family structure is detected
        if (scalar keys %pi_id > 3  ){
            print "Family larger than 3 encountered, program doesn't deal with it yet:\n", join("\n",keys %pi_id), "\n";
            exit 1;
        }
        
        # if all individuals of family of three are related (eg parent, two offspring), keep individual with best genotyping rate
        if (sum (values %pi_id) == 6 ) {
            my @missrate = ();
			my @trio = sort keys %pi_id;
            #print "All related: ", join ("\t", @trio), "\n";
			
			# iterate over individuals in trio and get genotyping info
            foreach my $id (@trio) {
                push @missrate, $imiss{$id};
            }

			# remove individual with best genotyping rate from array and map the remaining to keys in fail_id hash
			my $best_genotyping = argmin (\@missrate);
			my $id_keep = splice(@trio, $best_genotyping, 1);
            
			foreach my $id (@trio) {
				print OUT "$id $id\n";
                print OUT_related join("\t", ($id, $id_keep,  $pi{$id.$id_keep}, $id2fam{$id}, $id2fam{$id_keep}, $imiss{$id},$imiss{$id_keep}, $diagnosis{$id},  $diagnosis{$id_keep})), "\n" unless exists $fail_id{$id};
				$fail_id{$id} = 1;
			}
        } else {

        # if only one individual is related to both other individuals (eg both parents one offspring), exclude single individual
		    my @trio = sort keys %pi_id;
			my @numberrel = ();

            #print join("\t", @trio), "\n";
            # iterate over individuals in trio and get number of relatives
            foreach my $id (@trio) {
                push @numberrel, $pi_id{$id};
            }
            my $mostrel = argmax (\@numberrel);
			my $id_discard = splice(@trio, $mostrel, 1);
			print OUT "$id_discard $id_discard\n";
            print OUT_related join("\t", ($id_discard, $trio[0],  $pi{$id_discard.$trio[0]}, $id2fam{$id_discard}, $id2fam{$trio[0]}, $imiss{$id_discard}, $imiss{$trio[0]},  $diagnosis{$id_discard},  $diagnosis{$trio[0]})), "\n" unless exists $fail_id{$id_discard};
            print OUT_related join("\t", ($id_discard, $trio[1],  $pi{$id_discard.$trio[1]}, $id2fam{$id_discard}, $id2fam{$trio[1]}, $imiss{$id_discard}, $imiss{$trio[1]},  $diagnosis{$id_discard},  $diagnosis{$trio[1]})), "\n" unless exists $fail_id{$id_discard};
            $fail_id{$id_discard} = 1;
        }
    }
}

### c) filter samples that passed PI_HAT threshold based on SNP missingness rate
foreach my $id1 (keys %pihat) {
    foreach my $id2 (keys  %{$pihat{$id1}}) {

        next if exists $fail_id{$id1} || exists $fail_id{$id2};
        
        # get famID if individuals are part of registered family
        my ($famID1, $famID2) = ("", "");
        $famID1 = $id2fam{$id1} if exists $id2fam{$id1};
        $famID2 = $id2fam{$id2} if exists $id2fam{$id2};
        
		# check which individual has higher SNP missingness rate
        if ($imiss{$id1} >= $imiss{$id2}) {
            print OUT "$id1 $id1\n";
            print OUT_related join("\t", ($id1, $id2, $pihat{$id1}{$id2}, $famID1, $famID2, $imiss{$id1}, $imiss{$id2},  $diagnosis{$id1},  $diagnosis{$id2})),"\n";
            $fail_id{$id1} = 1;
        } else {
            print OUT "$id2 $id2\n";
            print OUT_related join ("\t", ($id2, $id1, $pihat{$id1}{$id2}, $famID1, $famID2, $imiss{$id2}, $imiss{$id1},  $diagnosis{$id2},  $diagnosis{$id1})),"\n";
            $fail_id{$id2} = 1;
        }
    }
}
close(OUT);

__END__
