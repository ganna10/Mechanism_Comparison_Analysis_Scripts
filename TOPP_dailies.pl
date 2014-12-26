#! /usr/bin/perl
# print out daily TOPP values on one plot
# Version 0: Jane Coates 10/02/2014
# Version 1: Jane Coates 10/11/2014 aesthetic updates
# Version 2: Jane Coates 8/12/2014 script updates after constant emissions runs

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper";
opendir DIR, $base or die "Can't open $base : $!";
my @daily_TOPP_files = grep { $_ =~ /_TOPP_values/ } readdir DIR;
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
        $TOPP{$mechanism}{$VOC} = \@TOPPs;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [1..7]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %TOPP) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $VOC (sort keys %{$TOPP{$mechanism}}) {
        $R->set('voc', $VOC);
        $R->set('topps', [@{$TOPP{$mechanism}{$VOC}}]);
        $R->run(q` pre[voc] = topps `);
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, VOC, TOPP, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `,
        q` data$VOC = factor(data$VOC, levels = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Ethylbenzene")) `,
        q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = TOPP, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ VOC ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("TOPP (molecules(Ox) / molecules(VOC))") `,
        q` plot = plot + scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, 1)) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + theme(plot.title = element_blank()) `,
        q` plot = plot + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) `,
        q` plot = plot + theme(legend.margin = unit(0, "lines")) `,
        q` plot = plot + scale_colour_manual(values = my.colours, guide = guide_legend(nrow = 1)) `,
);

$R->run(q` CairoPDF(file = "TOPP_daily_time_series_all_VOC.pdf", width = 8, height = 5.7) `,
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
