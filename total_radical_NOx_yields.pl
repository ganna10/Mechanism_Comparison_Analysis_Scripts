#! /usr/bin/env perl 
# Calculate total radical and NOx yields from each mechanism
# Version 0: Jane Coates 29/9/2014

use strict;
use diagnostics;
use KPP;

my $base = "/work/users/jco/MECCA";
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );

my (%families, %weights);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{'radicals'} = [ @radicals ];
    $families{'NOx'} = [ qw(NO NO2 NO3 N2O5) ];
    $weights{'NOx'} = { N2O5 => 2 };
    print "\n$run\n";
    my %yield;
    foreach my $species (qw(radicals NOx)) {
        print "\t$species\n";
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        my $producers = $kpp->producing($species);
        my $producer_yields = $kpp->effect_on($species, $producers);

        my $yield;
        foreach my $i (0..$#$producers) {
            $yield += $producer_yields->[$i];
        }
        $yield{$species} = $yield;
        print "\tTotal $species yield => $yield\n";
    }
    print "Net radical yield => ", $yield{'radicals'} - $yield{'NOx'}, "\n";
}

sub get_species { 
    my ($file) = @_; 

    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my %hash;
    my @lines = <$in>;
    close $in;
    my @separate = map { split /\s/, $_ } @lines;
    $hash{$_} += 1 foreach (@separate);
    return keys %hash;
} 
