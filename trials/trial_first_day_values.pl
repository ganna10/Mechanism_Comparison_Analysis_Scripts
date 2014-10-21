#! /usr/bin/perl
# compare first day TOPP values of all VOCs from all mechanisms to MCM v3.2 TOPP values in a separate facet including y = x line
# trying to plot everything on one plot, different colours for the NMVOCs and different shapes for the mechansisms
# Version 0: Jane Coates 14/02/2014

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/work/users/jco/MECCA/Mechanism_Comparison/TOPP_values_comparison";
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
$TOPP{'MOZART-4_NC5H12'} = $TOPP{'MOZART-4_BIGALK'};
delete $TOPP{'MOZART-4_BIGALK'};
$TOPP{'MOZART-4_BUT1ENE'} = $TOPP{'MOZART-4_BIGENE'};
delete $TOPP{'MOZART-4_BIGENE'};
$TOPP{'MOZART-4_C5H8'} = $TOPP{'MOZART-4_ISOP'};
delete $TOPP{'MOZART-4_ISOP'};

$TOPP{'RADM2_C2H6'} = $TOPP{'RADM2_ETH'};
delete $TOPP{'RADM2_ETH'};
$TOPP{'RADM2_C3H8'} = $TOPP{'RADM2_HC3'};
delete $TOPP{'RADM2_HC3'};
$TOPP{'RADM2_NC5H12'} = $TOPP{'RADM2_HC5'};
delete $TOPP{'RADM2_HC5'};
$TOPP{'RADM2_NC8H18'} = $TOPP{'RADM2_HC8'};
delete $TOPP{'RADM2_HC8'};
$TOPP{'RADM2_C2H4'} = $TOPP{'RADM2_OL2'};
delete $TOPP{'RADM2_OL2'};
$TOPP{'RADM2_C3H6'} = $TOPP{'RADM2_OLT'};
delete $TOPP{'RADM2_OLT'};
$TOPP{'RADM2_MEPROPENE'} = $TOPP{'RADM2_OLI'};
delete $TOPP{'RADM2_OLI'};
$TOPP{'RADM2_C5H8'} = $TOPP{'RADM2_ISO'};
delete $TOPP{'RADM2_ISO'};
$TOPP{'RADM2_TOLUENE'} = $TOPP{'RADM2_TOL'};
delete $TOPP{'RADM2_TOL'};
$TOPP{'RADM2_MXYL'} = $TOPP{'RADM2_XYL'};
delete $TOPP{'RADM2_XYL'};

$TOPP{'RACM_C2H6'} = $TOPP{'RACM_ETH'};
delete $TOPP{'RACM_ETH'};
$TOPP{'RACM_C3H8'} = $TOPP{'RACM_HC3'};
delete $TOPP{'RACM_HC3'};
$TOPP{'RACM_NC5H12'} = $TOPP{'RACM_HC5'};
delete $TOPP{'RACM_HC5'};
$TOPP{'RACM_NC8H18'} = $TOPP{'RACM_HC8'};
delete $TOPP{'RACM_HC8'};
$TOPP{'RACM_C2H4'} = $TOPP{'RACM_ETE'};
delete $TOPP{'RACM_ETE'};
$TOPP{'RACM_C3H6'} = $TOPP{'RACM_OLT'};
delete $TOPP{'RACM_OLT'};
$TOPP{'RACM_MEPROPENE'} = $TOPP{'RACM_OLI'};
delete $TOPP{'RACM_OLI'};
$TOPP{'RACM_C5H8'} = $TOPP{'RACM_ISO'};
delete $TOPP{'RACM_ISO'};
$TOPP{'RACM_TOLUENE'} = $TOPP{'RACM_TOL'};
delete $TOPP{'RACM_TOL'};
$TOPP{'RACM_MXYL'} = $TOPP{'RACM_XYL'};
delete $TOPP{'RACM_XYL'};

$TOPP{'RACM2_C2H6'} = $TOPP{'RACM2_ETH'};
delete $TOPP{'RACM2_ETH'};
$TOPP{'RACM2_C3H8'} = $TOPP{'RACM2_HC3'};
delete $TOPP{'RACM2_HC3'};
$TOPP{'RACM2_NC5H12'} = $TOPP{'RACM2_HC5'};
delete $TOPP{'RACM2_HC5'};
$TOPP{'RACM2_NC8H18'} = $TOPP{'RACM2_HC8'};
delete $TOPP{'RACM2_HC8'};
$TOPP{'RACM2_C2H4'} = $TOPP{'RACM2_ETE'};
delete $TOPP{'RACM2_ETE'};
$TOPP{'RACM2_C3H6'} = $TOPP{'RACM2_OLT'};
delete $TOPP{'RACM2_OLT'};
$TOPP{'RACM2_MEPROPENE'} = $TOPP{'RACM2_OLI'};
delete $TOPP{'RACM2_OLI'};
$TOPP{'RACM2_C5H8'} = $TOPP{'RACM2_ISO'};
delete $TOPP{'RACM2_ISO'};
$TOPP{'RACM2_BENZENE'} = $TOPP{'RACM2_BEN'};
delete $TOPP{'RACM2_BEN'};
$TOPP{'RACM2_TOLUENE'} = $TOPP{'RACM2_TOL'};
delete $TOPP{'RACM2_TOL'};
$TOPP{'RACM2_MXYL'} = $TOPP{'RACM2_XYM'};
delete $TOPP{'RACM2_XYM'};
$TOPP{'RACM2_OXYL'} = $TOPP{'RACM2_XYO'};
delete $TOPP{'RACM2_XYO'};
$TOPP{'RACM2_PXYL'} = $TOPP{'RACM2_XYP'};
delete $TOPP{'RACM2_XYP'};

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(plyr) `);
$R->run(q` library(reshape2) `);
$R->run(q` library(reshape) `);
$R->run(q` library(grid) `);
$R->run(q` library(gridExtra) `);
$R->run(q` library(gtable) `);
$R->run(q` library(Cairo) `);

$R->run(q` NMVOC = {} `,
        q` Mechanism = {} `,
        q` TOPP.value = {} `,
);

foreach my $species (sort keys %TOPP) {
    my ($mechanism, $name) = split /_/, $species;
    my $R_name = $R->set('name', $name);
    if ($mechanism =~ /MOZART/) {
        my $R_mech = $R->set('mechanism', 'MOZART');
    } elsif ($mechanism =~ /CBM/) {
        my $R_mech = $R->set('mechanism', 'CBM4');
    } else {
        my $R_mech = $R->set('mechanism', $mechanism);
    }
    my $R_data = $R->set('TOPP', $TOPP{$species});
    $R->run(q` NMVOC = cbind(NMVOC, name) `,
            q` Mechanism = cbind(Mechanism, mechanism) `,
            q` TOPP.value = cbind(TOPP.value, TOPP) `,
    );
}

#create dataframe after converting the matrices above to vectors
$R->run(q` NMVOC = c(NMVOC) `,
        q` Mechanism = c(Mechanism) `,
        q` TOPP.value = c(TOPP.value) `,
        q` data = data.frame(NMVOC, Mechanism, TOPP.value) `,
);

#manually specify order of NMVOC
$R->run(q` casted = cast(data, NMVOC ~ Mechanism, value = 'TOPP.value') `);
$R->run(q` melted = melt(casted, id.vars = c("NMVOC", "MCMv3.2"), variable.name = "TOPPs", variable.name = "Mechanism") `,
    #q` colnames(melted)[3:4] = c("Mechanism", "TOPP.value")`,
);
#$R->run(q` melted$NMVOC <- factor(casted$NMVOC, levels = c("C2H6", "C3H8", "NC4H10", "IC4H10", "NC5H12", "IC5H12", "NC6H14", "NC7H16", "NC8H18", "C2H4", "C3H6", "BUT1ENE", "MEPROPENE", "C5H8", "BENZENE", "TOLUENE", "MXYL", "OXYL", "PXYL", "EBENZ"), ordered = TRUE) `);

my $print = $R->run(q` print(melted) `);
print $print, "\n";
#specifiying colours, shapes and names
$R->run(q` my.colours = c( "C2H6" = "#696537", "C3H8" = "#f9c600", "NC4H10" = "#76afca", "IC4H10" = "#dc3522", "NC5H12" = "#8c6238", "IC5H12" = "#9bb08f", "NC6H14" = "#8b1537", "NC7H16" = "#ba8b01", "NC8H18" = "#0352cb", "C2H4" = "#86b650", "C3H6" = "#6c254f", "BUT1ENE" = "#ee6738", "MEPROPENE" = "#58691b", "C5H8" = "#8ed6d5", "BENZENE" = "#f3aa7f", "TOLUENE" = "#c65d6c", "MXYL" = "#888a87", "OXYL" = "#0e5c28", "PXYL" = "#b569b3", "EBENZ" = "#2c9def" ) `,
        q` my.names = c("C2H6" = "Ethane", "C3H8" = "Propane", "NC4H10" = "Butane", "IC4H10" = "2-Methylpropane", "NC5H12" = "Pentane", "IC5H12" = "2-Methylbutane", "NC6H14" = "Hexane", "NC7H16" = "Heptane", "NC8H18" = "Octane", "C2H4" = "Ethene", "C3H6" = "Propene", "BUT1ENE" = "Butene", "MEPROPENE" = "2-Methylpropene", "C5H8" = "Isoprene", "BENZENE" = "Benzene", "TOLUENE" = "Toluene", "MXYL" = "m-Xylene", "OXYL" = "o-Xylene", "PXYL" = "p-Xylene", "EBENZ" = "Ethylbenzene") `,
        q` my.shapes = c(0, 1, 2, 3, 4, 6, 7, 9) `,
);

$R->run(q` plot.lines = function () { list( geom_point(size = 20), xlab("\nMCM v3.2 TOPP (molecules(Ox)/molecules(VOC))"), scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), scale_y_continuous(limits = c(0,7), breaks = seq(0, 7, 1)), geom_abline(intercept = 0, slope = 1), theme_bw(), theme(axis.text.x = element_text(size = 110), axis.text.y = element_text(size = 110), axis.title.x = element_text(size = 120, face = "bold"), axis.title.y = element_text(size = 120, face = "bold")), scale_colour_manual(labels = my.names, values = my.colours, guide = guide_legend(label.position = "bottom", label.theme = element_text(size = 170, angle = 45), label.hjust = 0.5, label.vjust = 0.9)), theme(legend.position = "bottom", legend.key.size = unit(10, "cm"), legend.key = element_blank()), scale_shape_manual(values = my.shapes) ) } `); 

$R->run(q` plot = ggplot(data = melted, aes(x = MCMv3.2, y = TOPP.value, colour = NMVOC, shaped = Mechanism)) `,
        q` plot = plot + plot.lines() `, 
        q` plot = plot + ylab("\nTOPP (molecules(Ox)/molecules(VOC))\n") `,
        
        q` CairoPDF(file = "trial_first_day_values.pdf", width = 100, height = 100) `,

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
            $TOPP{"${index}_$species"} = $TOPP[0];
        }
        close FILE;
    }
 }
