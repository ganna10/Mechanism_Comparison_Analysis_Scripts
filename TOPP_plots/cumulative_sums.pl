#! /usr/bin/perl
# calculate cumulative sums of TOPP values
# Version 0: Jane Coates 5/1/2015

use strict;
use diagnostics;

my $base = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper/TOPP_plots";
opendir DIR, $base or die "Can't open $base : $!";
my @daily_TOPP_files = grep { $_ =~ /_TOPP_values\.txt/ } readdir DIR;
closedir DIR;
my %TOPP;

foreach my $file (@daily_TOPP_files) {
    my $path = "$base/$file";
    my @lines = split /\n/, read_file($path);
    (my $mechanism = $file) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $VOC = get_chemical_name($VOC);
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $TOPP{$mechanism}{$VOC} = \@TOPPs;
    }
}

#calculate cumulative sums of TOPPs
my %TOPP_cumsum;
foreach my $mechanism (keys %TOPP) {
    foreach my $VOC (sort keys %{$TOPP{$mechanism}}) {
        my @TOPP_sum;
        foreach my $TOPP ($TOPP{$mechanism}{$VOC}) {
            my $sum = 0;
            @TOPP_sum = map { $sum += $_ } @$TOPP;
        }
        $TOPP_cumsum{$mechanism}{$VOC} = \@TOPP_sum;
    }
}

my $output = "_TOPP_cumulative.txt";
foreach my $mechanism (sort keys %TOPP_cumsum) {
    my $file = $mechanism . $output;
    open my $out, '>:encoding(utf-8)', $file or die "Can't open $file : $!";
    foreach my $VOC (sort keys %{$TOPP_cumsum{$mechanism}}) {
        my @topps = @{$TOPP_cumsum{$mechanism}{$VOC}};
        print $out "$VOC => @topps\n";
    }
    close $out;
}

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
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
