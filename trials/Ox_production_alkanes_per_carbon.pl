#! /usr/bin/env perl
# Amount of Ox produced per C during alkane degradation
# Version 0: Jane Coates 14/1/2015

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper/TOPP_plots";
opendir DIR, $base or die "Can't open $base : $!";
my @daily_TOPP_files = grep { $_ =~ /_TOPP_values/ } readdir DIR;
closedir DIR;

foreach my $file (@daily_TOPP_files) {
    my $path = "$base/$file";
    my @lines = split /\n/, read_file($path);
    (my $mechanism = $file) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next unless ($line =~ /C2H6|C3H8|NC|IC|ETH|HC|BIGALK/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $VOC = get_chemical_name($VOC);
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $TOPP{$mechanism}{$VOC} = $TOPPs[0];
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` data = data.frame(Mechanism = as.numeric(0), VOC = as.numeric(0), TOPP.per.C = as.numeric(0)) `); 
foreach my $mechanism (sort keys %TOPP) {
    $R->set('mechanism', $mechanism);
    foreach my $alkane (sort keys %{$TOPP{$mechanism}}) {
        if ($alkane eq "Ethane") {
            $TOPP{$mechanism}{$alkane} /= 2;
        } elsif ($alkane eq "Propane") {
            $TOPP{$mechanism}{$alkane} /= 3;
        } elsif ($alkane eq "Butane" or $alkane eq "2-Methylpropane") {
            $TOPP{$mechanism}{$alkane} /= 4;
        } elsif ($alkane eq "Pentane" or $alkane eq "2-Methylbutane") {
            $TOPP{$mechanism}{$alkane} /= 5;
        } elsif ($alkane eq "Hexane") {
            $TOPP{$mechanism}{$alkane} /= 6;
        } elsif ($alkane eq "Heptane") {
            $TOPP{$mechanism}{$alkane} /= 7;
        } elsif ($alkane eq "Octane") {
            $TOPP{$mechanism}{$alkane} /= 8;
        } else {
            print "No C for $alkane\n";
        }
        $R->set('alkane', $alkane);
        $R->set('topp.per.c', $TOPP{$mechanism}{$alkane});
        $R->run(q` data = rbind(data, c(mechanism, alkane, topp.per.c)) `);
    } 
}
$R->run(q` data = data[-1,] `);
my $out = "Alkanes_TOPP_per_C.txt";
my $p = $R->run(q` print(data) `);
open my $file, ">:encoding(utf-8)", $out or die $!;
print $file $p;
close $file;

$R->stop();

sub read_file {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    close $in;
    return $data;
} 

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
