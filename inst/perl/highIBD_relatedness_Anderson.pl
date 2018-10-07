#!/usr/bin/perl

use strict;

my %imiss;
my %removed;

my $verbose = $ARGV[1];

open IMISS, '<', $ARGV[0].".imiss"
        or die "Cannot open genotypes file (".$ARGV[0].".genome): $!\n";
if($verbose eq "verbose") {
    print "Reading PLINK .imiss file ".$ARGV[0].".imiss\n";
}
while(<IMISS>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', $ARGV[0].".genome"
        or die "Cannot open genotypes file (".$ARGV[0].".genome): $!\n";
# Adapted output file name
#open OUT, '>', "fail-IBD-QC.txt";
open OUT, '>', $ARGV[0].".fail-IBD.IDs";
if($verbose eq "verbose") {
    print "Reading PLINK .genome file ".$ARGV[0].".genome\n";
}
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
    if($fields[9] > 0.185){
        if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
            unless($removed{$fields[0]}{$fields[1]}){
                print OUT "$fields[0] $fields[1]\n";
                $removed{$fields[0]}{$fields[1]} = 1;
            }
        }
        elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
            unless($removed{$fields[2]}{$fields[3]}){
                print OUT "$fields[2] $fields[3]\n";
                $removed{$fields[2]}{$fields[3]} = 1;
            }
        }
        else{
            unless($removed{$fields[0]}{$fields[1]}){
                print OUT "$fields[0] $fields[1]\n";
                $removed{$fields[0]}{$fields[1]} = 1;
            }
        }
    }
}
    
    

