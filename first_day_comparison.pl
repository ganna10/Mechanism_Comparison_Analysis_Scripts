#! /usr/bin/perl
# compare first day TOPP values of all VOCs from all mechanisms to MCM v3.2 TOPP values in a separate facet including y = x line
# Version 0: Jane Coates 14/02/2014
# Version 1: Jane Coates 7/11/2014 updating plot to include title of each facet rather than geom_text
# Version 2: Jane Coates 8/12/2014 script updates for constant emissions runs

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper";
opendir DIR, $base or die "Can't open $base : $!";
my @daily_TOPP_files = grep { $_ =~ /TOPP_values/ } readdir DIR;
closedir DIR;

foreach my $file (@daily_TOPP_files) {
    my @lines = split /\n/, read_file($file);
    (my $mechanism = $file) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $VOC = get_chemical_name($VOC);
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $TOPP{$mechanism}{$VOC} = $TOPPs[0];
    }
}

my $R = Statistics::R->new();
$R->run(q` library(dplyr) `);
$R->run(q` library(tidyr) `);
$R->run(q` library(ggplot2) `);
$R->run(q` library(Cairo) `);

$R->set('VOCs', [sort keys %{$TOPP{"MCMv3.2"}}]);
$R->run(q` data = data.frame(VOC = VOCs) `); 
foreach my $mechanism (sort keys %TOPP) {
    $R->set('mechanism', $mechanism);
    $R->run(q` topp = c() `);
    foreach my $VOC (sort keys %{$TOPP{"MCMv3.2"}}) { 
        $R->set('value', $TOPP{$mechanism}{$VOC});
        $R->run(q` topp = c(topp, value) `);
    }
    $R->run(q` data[mechanism] = topp `);
}
$R->run(q` data = gather(data, Mechanism, TOPP, -VOC, -MCMv3.2) `);
#my $p = $R->run(q` print(label.text) `);
#print "$p\n";

#specifiying colours and names
$R->run(q` my.colours = c( "Ethane " = "#696537", "Propane " = "#f9c600", "Butane " = "#76afca", "2-Methylpropane " = "#dc3522", "Pentane " = "#8c6238", "2-Methylbutane " = "#9bb08f", "Hexane " = "#8b1537", "Heptane " = "#ba8b01", "Octane " = "#0352cb", "Ethene " = "#86b650", "Propene " = "#6c254f", "Butene " = "#ee6738", "2-Methylpropene " = "#58691b", "Isoprene " = "#8ed6d5", "Benzene " = "#f3aa7f", "Toluene " = "#c65d6c", "m-Xylene " = "#888a87", "o-Xylene " = "#0e5c28", "p-Xylene " = "#b569b3", "Ethylbenzene " = "#2c9def" ) `,
        q` data$VOC = factor(data$VOC, levels = c("Ethane ", "Propane ", "Butane ", "2-Methylpropane ", "Pentane ", "2-Methylbutane ", "Hexane ", "Heptane ", "Octane ", "Ethene ", "Propene ", "Butene ", "2-Methylpropene ", "Isoprene ", "Benzene ", "Toluene ", "m-Xylene ", "o-Xylene ", "p-Xylene ", "Ethylbenzene ")) `,
        q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `,
);

$R->run(q` plot = ggplot(data = data, aes(x = MCMv3.2, y = TOPP, colour = VOC, group = VOC)) `,
        q` plot = plot + geom_point(shape = 19) `,
        q` plot = plot + facet_wrap( ~ Mechanism, ncol = 2) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + xlab("MCM v3.2 TOPP (molecules(Ox)/molecules(VOC))") `,
        q` plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC))") `,
        q` plot = plot + geom_abline(intercept = 0, slope = 1) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + scale_colour_manual(values = my.colours, guide = guide_legend(nrow = 3)) `,
);

$R->run(q` CairoPDF(file = "first_day_values.pdf", width = 8, height = 11.3) `,
        q` print(plot) `,
        q` dev.off() `,
);

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
        $chemical_species = 'Ethane ';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane ';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane ';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane ';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane ';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane ';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane ';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane ";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane ";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene ';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene ';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene ";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene ';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene ";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene ";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene ';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene ";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene ';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene ";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene ";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
