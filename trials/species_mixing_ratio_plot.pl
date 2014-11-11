#! /usr/bin/env perl
# Plot mixing ratio time series of species specified in ARGV. Species name is MCM name.
# Version 0: Jane Coates 17/9/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $mcm_species = $ARGV[0];
my $base = "/work/users/jco/MECCA";
my %mixing_ratios;

my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
$times = $times(1:$NTIME-2);
$times -= $times->at(0);
$times /= 86400;
my @time_axis = map { $_ } $times->dog;

my @dirs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
#my @dirs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my $index = 0;

foreach my $run (@dirs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $spc_file = "$base/$run/gas.spc";
    my $mech_species = get_mechanism_species($mcm_species, $mechanisms[$index]);
    my $tagged_mechanism_species = get_tagged_species($mech_species, $spc_file);
    $mixing_ratios{$mechanisms[$index]} = get_mixing_ratios($tagged_mechanism_species, $mecca);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('Time', [@time_axis]);
$R->run(q` data = data.frame(Time) `);

foreach my $run (sort keys %mixing_ratios) {
    $R->set('Mechanism', $run); 
    $R->set('Mixing.ratios', [@{$mixing_ratios{$run}}]);
    $R->run(q` data[Mechanism] = Mixing.ratios `);
}

$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Mechanism", value.name = "Mixing.ratio") `); #arrange data
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->set('plot.title', "${mcm_species} Mixing Ratio Comparison");
$R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);
$R->run(q` plot = ggplot(data = data, aes(x = Time, y = Mixing.ratio, colour = Mechanism)) `, #plotting
        q` plot = plot + geom_line(size = 1) `,
        q` plot = plot + ggtitle(plot.title) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + scale_x_continuous(limits = c(0,7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(0.99, 0.99)) `,
        q` plot = plot + theme(legend.justification = c(0.99, 0.99)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->set('filename', "${mcm_species}_mixing_ratio_comparison.pdf");
$R->run(q` CairoPDF(file = filename, width = 10, height = 7) `, #print to file
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

sub get_mixing_ratios {
    my ($all_species, $mecca) = @_;
    my $mixing_ratio;
    $mixing_ratio += $mecca->tracer($_) foreach (@$all_species);
    $mixing_ratio = $mixing_ratio(1:$NTIME-2);
    $mixing_ratio *= 1e9; # convert to ppbv
    my @mixing_ratios = map { $_ } $mixing_ratio->dog;
    return \@mixing_ratios;
}

sub get_tagged_species {
    my ($species, $spc_file) = @_;
    my $all_species = read_file($spc_file);
    my @tagged_species;
    foreach my $line (@$all_species) {
        next unless ($line =~ /^$species/);
        $line =~ s/\s=.*$//;
        push @tagged_species, $line;
    }
    return \@tagged_species;
}

sub get_mechanism_species {
    my ($mcm_species, $mechanism) = @_;
    my $species;

    if ($mcm_species eq "HCHO") {
        if ($mechanism =~ /MCM|CRI|RADM2|RACM|CBM/) {
            $species = "HCHO";
        } elsif ($mechanism =~ /MOZART/) {
            $species = "CH2O";
        } elsif ($mechanism =~ /CB05/) {
            $species = "FORM";
        }
    } elsif ($mcm_species eq "HO2") {
        $species = "HO2";
    } else {
        die "No mapping for $mcm_species\n";
    }
    return $species;
}

sub read_file {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file for reading : $!";
    chomp(my @all = <$in>);
    close $in;
    return \@all;
} 
