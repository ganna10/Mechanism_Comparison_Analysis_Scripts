#! /usr/bin/perl
# compare first day TOPP values of all VOCs from all mechanisms to MCM v3.2 TOPP values in a separate facet including y = x line
# Version 0: Jane Coates 14/02/2014
# Version 1: Jane Coates 7/11/2014 updating plot to include title of each facet rather than geom_text

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/local/home/coates/MECCA/Mechanism_Comparison/TOPP_values_comparison";
my $mcm3_1_file = "$base/MCMv3.1_TOPP_values.txt";
my $mcm3_2_file = "$base/MCMv3.2_TOPP_values.txt";
my $cri_file = "$base/CRIv2_TOPP_values.txt";
my $mozart_file = "$base/MOZART-4_TOPP_values.txt";
my $radm2_file = "$base/RADM2_TOPP_values.txt";
my $racm_file = "$base/RACM_TOPP_values.txt";
my $racm2_file = "$base/RACM2_TOPP_values.txt";
my $cbm4_file = "$base/CBM-IV_TOPP_values.txt";
my $cb05_file = "$base/CB05_TOPP_values.txt";

my @file_list = ($mcm3_1_file, $mcm3_2_file, $cri_file, $mozart_file, $radm2_file, $racm_file, $racm2_file, $cbm4_file, $cb05_file);
get_TOPPs(@file_list);

#re-assign species to MCM v3.2 species names
$TOPP{'NC5H12_MOZART-4'} = $TOPP{'BIGALK_MOZART-4'};
delete $TOPP{'BIGALK_MOZART-4'};
$TOPP{'BUT1ENE_MOZART-4'} = $TOPP{'BIGENE_MOZART-4'};
delete $TOPP{'BIGENE_MOZART-4'};
$TOPP{'C5H8_MOZART-4'} = $TOPP{'ISOP_MOZART-4'};
delete $TOPP{'ISOP_MOZART-4'};

$TOPP{'C2H6_RADM2'} = $TOPP{'ETH_RADM2'};
delete $TOPP{'ETH_RADM2'};
$TOPP{'C3H8_RADM2'} = $TOPP{'HC3_RADM2'};
delete $TOPP{'HC3_RADM2'};
$TOPP{'NC5H12_RADM2'} = $TOPP{'HC5_RADM2'};
delete $TOPP{'HC5_RADM2'};
$TOPP{'NC8H18_RADM2'} = $TOPP{'HC8_RADM2'};
delete $TOPP{'HC8_RADM2'};
$TOPP{'C2H4_RADM2'} = $TOPP{'OL2_RADM2'};
delete $TOPP{'OL2_RADM2'};
$TOPP{'C3H6_RADM2'} = $TOPP{'OLT_RADM2'};
delete $TOPP{'OLT_RADM2'};
$TOPP{'MEPROPENE_RADM2'} = $TOPP{'OLI_RADM2'};
delete $TOPP{'OLI_RADM2'};
$TOPP{'C5H8_RADM2'} = $TOPP{'ISO_RADM2'};
delete $TOPP{'ISO_RADM2'};
$TOPP{'TOLUENE_RADM2'} = $TOPP{'TOL_RADM2'};
delete $TOPP{'TOL_RADM2'};
$TOPP{'MXYL_RADM2'} = $TOPP{'XYL_RADM2'};
delete $TOPP{'XYL_RADM2'};

$TOPP{'C2H6_RACM'} = $TOPP{'ETH_RACM'};
delete $TOPP{'ETH_RACM'};
$TOPP{'C3H8_RACM'} = $TOPP{'HC3_RACM'};
delete $TOPP{'HC3_RACM'};
$TOPP{'NC5H12_RACM'} = $TOPP{'HC5_RACM'};
delete $TOPP{'HC5_RACM'};
$TOPP{'NC8H18_RACM'} = $TOPP{'HC8_RACM'};
delete $TOPP{'HC8_RACM'};
$TOPP{'C2H4_RACM'} = $TOPP{'ETE_RACM'};
delete $TOPP{'ETE_RACM'};
$TOPP{'C3H6_RACM'} = $TOPP{'OLT_RACM'};
delete $TOPP{'OLT_RACM'};
$TOPP{'MEPROPENE_RACM'} = $TOPP{'OLI_RACM'};
delete $TOPP{'OLI_RACM'};
$TOPP{'C5H8_RACM'} = $TOPP{'ISO_RACM'};
delete $TOPP{'ISO_RACM'};
$TOPP{'TOLUENE_RACM'} = $TOPP{'TOL_RACM'};
delete $TOPP{'TOL_RACM'};
$TOPP{'MXYL_RACM'} = $TOPP{'XYL_RACM'};
delete $TOPP{'XYL_RACM'};

$TOPP{'C2H6_RACM2'} = $TOPP{'ETH_RACM2'};
delete $TOPP{'ETH_RACM2'};
$TOPP{'C3H8_RACM2'} = $TOPP{'HC3_RACM2'};
delete $TOPP{'HC3_RACM2'};
$TOPP{'NC5H12_RACM2'} = $TOPP{'HC5_RACM2'};
delete $TOPP{'HC5_RACM2'};
$TOPP{'NC8H18_RACM2'} = $TOPP{'HC8_RACM2'};
delete $TOPP{'HC8_RACM2'};
$TOPP{'C2H4_RACM2'} = $TOPP{'ETE_RACM2'};
delete $TOPP{'ETE_RACM2'};
$TOPP{'C3H6_RACM2'} = $TOPP{'OLT_RACM2'};
delete $TOPP{'OLT_RACM2'};
$TOPP{'MEPROPENE_RACM2'} = $TOPP{'OLI_RACM2'};
delete $TOPP{'OLI_RACM2'};
$TOPP{'C5H8_RACM2'} = $TOPP{'ISO_RACM2'};
delete $TOPP{'ISO_RACM2'};
$TOPP{'BENZENE_RACM2'} = $TOPP{'BEN_RACM2'};
delete $TOPP{'BEN_RACM2'};
$TOPP{'TOLUENE_RACM2'} = $TOPP{'TOL_RACM2'};
delete $TOPP{'TOL_RACM2'};
$TOPP{'MXYL_RACM2'} = $TOPP{'XYM_RACM2'};
delete $TOPP{'XYM_RACM2'};
$TOPP{'OXYL_RACM2'} = $TOPP{'XYO_RACM2'};
delete $TOPP{'XYO_RACM2'};
$TOPP{'PXYL_RACM2'} = $TOPP{'XYP_RACM2'};
delete $TOPP{'XYP_RACM2'};

my @mechanisms_fill = qw( MOZART-4 RADM2 RACM RACM2 );
my @VOCs = map { local $_ = $_; s/_(.*?)$// ; $_ } grep { $_ =~ /MCMv3\.2/ } sort keys %TOPP;
foreach my $mechanism (@mechanisms_fill) {
    foreach my $VOC (@VOCs) {
        my $label = "${VOC}_$mechanism";
        $TOPP{$label} = 999 if (not defined $TOPP{$label});
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(dplyr) `);
$R->run(q` library(plyr) `);
$R->run(q` library(reshape2) `);
$R->run(q` library(grid) `);
$R->run(q` library(Cairo) `);

$R->set('VOCs', [map { local $_ = $_; s/_(.*?)$// ; $_ } grep { $_ =~ /MCMv3\.2/ } sort keys %TOPP]);
$R->run(q` data = data.frame(NMVOC = VOCs) `);

my @mechanisms = ( "MCMv3.2", "(a) MCMv3.1",  "(b) CRIv2", "(c) RADM2", "(d) RACM", "(e) RACM2", "(f) MOZART-4", "(g) CBM-IV", "(h) CB05" );
foreach my $mechanism (@mechanisms) { 
    (my $true = $mechanism) =~ s/\([a-z]\)\s//;
    $R->set('TOPPs', [ map { $TOPP{$_} } grep { $_ =~ /$true$/ } sort keys %TOPP]);
    $R->set('mechanism', $mechanism);
    $R->run(q` data[mechanism] = TOPPs `);
}

$R->run(q` data = melt(data, id.vars = c("NMVOC", "MCMv3.2"), variable.name = "Mechanism", value.name = "TOPP") `,
        q` data = as.data.frame(lapply(data, function (x) {replace(x, x == 999, NA)})) `,
);
#my $p = $R->run(q` print(data) `);
#print "$p\n";

#specifiying colours and names
$R->run(q` my.colours = c( "C2H6" = "#696537", "C3H8" = "#f9c600", "NC4H10" = "#76afca", "IC4H10" = "#dc3522", "NC5H12" = "#8c6238", "IC5H12" = "#9bb08f", "NC6H14" = "#8b1537", "NC7H16" = "#ba8b01", "NC8H18" = "#0352cb", "C2H4" = "#86b650", "C3H6" = "#6c254f", "BUT1ENE" = "#ee6738", "MEPROPENE" = "#58691b", "C5H8" = "#8ed6d5", "BENZENE" = "#f3aa7f", "TOLUENE" = "#c65d6c", "MXYL" = "#888a87", "OXYL" = "#0e5c28", "PXYL" = "#b569b3", "EBENZ" = "#2c9def" ) `,
        q` my.names = c("C2H6" = "Ethane", "C3H8" = "Propane", "NC4H10" = "Butane", "IC4H10" = "2-Methylpropane", "NC5H12" = "Pentane", "IC5H12" = "2-Methylbutane", "NC6H14" = "Hexane", "NC7H16" = "Heptane", "NC8H18" = "Octane", "C2H4" = "Ethene", "C3H6" = "Propene", "BUT1ENE" = "Butene", "MEPROPENE" = "2-Methylpropene", "C5H8" = "Isoprene", "BENZENE" = "Benzene", "TOLUENE" = "Toluene", "MXYL" = "m-Xylene", "OXYL" = "o-Xylene", "PXYL" = "p-Xylene", "EBENZ" = "Ethylbenzene") `,
        q` data$NMVOC = factor(data$NMVOC, levels = c("C2H6", "C3H8", "NC4H10", "IC4H10", "NC5H12", "IC5H12", "NC6H14", "NC7H16", "NC8H18", "C2H4", "C3H6", "BUT1ENE", "MEPROPENE", "C5H8", "BENZENE", "TOLUENE", "MXYL", "OXYL", "PXYL", "EBENZ")) `,
);

$R->run(q` plot = ggplot(data, aes(x = MCMv3.2, y = TOPP, colour = NMVOC)) `,
        q` plot = plot + geom_point(shape = 19, size = 55) `,
        q` plot = plot + facet_wrap( ~ Mechanism, ncol = 2) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + xlab("\nMCM v3.2 TOPP (molecules(Ox)/molecules(VOC))\n") `,
        q` plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC))\n") `,
        q` plot = plot + geom_abline(intercept = 0, slope = 1, size = 5) `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 180, face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 180, face = "bold")) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(panel.margin = unit(1, "cm")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(axis.ticks.length = unit(2, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + scale_colour_manual(values = my.colours, labels = my.names, guide = guide_legend(label.position = "bottom", label.theme = element_text(size = 145, angle = 45), label.hjust = 0.5, label.vjust = 0.9)) `,
);

$R->run(q` CairoPDF(file = "first_day_values.pdf", width = 152, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_TOPPs {
    my (@files) = @_;
    
    foreach my $file (@files) {
        my ($filename) = $file =~ m|comparison/(.*?)\.txt$|; 
        my ($index, $extra) = split /_/, $filename;
        open FILE, "<$file" or die $!;
        foreach (<FILE>) {
            next if ($_ =~ /^Working/ or $_ =~ /^CH4/); 
            $_ =~ s/\[|\]//g;
            my ($species, $TOPP) = split /=>/, $_;
            foreach my $item ($species, $TOPP) {
                $item =~ s/^\s+|\s+$//g;
            }
            $TOPP =~ s/\s+/,/g;
            my @TOPP = split /,/, $TOPP;
            $species .= "_$index"; 
            $TOPP{$species} = $TOPP[0];
        }
        close FILE;
    }
}
