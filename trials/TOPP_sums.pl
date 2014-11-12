#! /usr/bin/perl
# print out cumulative sums of TOPP values on one plot
# Version 0: Jane Coates 25/06/2014

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;

my $base = "/work/users/jco/MECCA/Mechanism_Comparison/TOPP_sum_comparison";
my $mcm3_1_file = "$base/MCMv3.1_TOPP_sums.txt";
my $mcm3_2_file = "$base/MCMv3.2_TOPP_sums.txt";
my $cri_file = "$base/CRIv2_TOPP_sums.txt";
my $mozart_file = "$base/MOZART-4_TOPP_sums.txt";
my $radm2_file = "$base/RADM2_TOPP_sums.txt";
my $racm_file = "$base/RACM_TOPP_sums.txt";
my $racm2_file = "$base/RACM2_TOPP_sums.txt";
my $cbm4_file = "$base/CBM-IV_TOPP_sums.txt";
my $cb05_file = "$base/CB05_TOPP_sums.txt";

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

get_plot(\%TOPP);

sub get_plot {
    my ($plot_data) = @_;
    my %plot_hash = %$plot_data;
    
    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(gtable) `);
    $R->run(q` library(Cairo) `);

    #time axis
    my $n_days = 7;
    my @time_axis = (1..$n_days);
    $R->set('times', [@time_axis]); 
    $R->set('rep.number', scalar(keys %plot_hash));

    $R->run(q` times = rep(times, rep.number) `,
            q` NMVOC = {} `,
            q` Mechanism = {} `,
            q` TOPP.value = {} `,
    );

    foreach my $species (sort keys %plot_hash) {
        my ($name, $mechanism) = split /_/, $species;
        my $R_name = $R->set('name', $name);
        my $R_mech = $R->set('mechanism', $mechanism);
        my $R_data = $R->set('TOPP', [@{$plot_hash{$species}}]);
        $R->run(q` NMVOC = cbind(NMVOC, rep(name, 7)) `,
                q` Mechanism = cbind(Mechanism, rep(mechanism, 7)) `,
                q` TOPP.value = cbind(TOPP.value, TOPP) `,
        );
    }

    #create dataframe after converting the matrices above to vectors
    $R->run(q` NMVOC = c(NMVOC) `,
            q` Mechanism = c(Mechanism) `,
            q` TOPP.value = c(TOPP.value) `,
            q` data = data.frame(times, NMVOC, Mechanism, TOPP.value) `,
    );

    #order facets by C-number
    $R->run(q` data$NMVOC = factor(data$NMVOC, levels = c('C2H6', 'C3H8', 'NC4H10', 'IC4H10', 'NC5H12', 'IC5H12', 'NC6H14', 'NC7H16', 'NC8H18', 'C2H4', 'C3H6', 'BUT1ENE', 'MEPROPENE', 'C5H8', 'BENZENE', 'TOLUENE', 'MXYL', 'PXYL', 'OXYL', 'EBENZ')) `);

    #change facet labels
    $R->run(q` facet_labels = list( 'C2H6' = 'Ethane',
                                    'C3H8' = 'Propane',
                                    'NC4H10' = 'Butane',
                                    'IC4H10' = '2-Methylpropane',
                                    'NC5H12' = 'Pentane',
                                    'IC5H12' = '2-Methylbutane',
                                    'NC6H14' = 'Hexane',
                                    'NC7H16' = 'Heptane',
                                    'NC8H18' = 'Octane',
                                    'C2H4' = 'Ethene',
                                    'C3H6' = 'Propene',
                                    'BUT1ENE' = 'Butene',
                                    'MEPROPENE' = '2-Methylpropene',
                                    'C5H8' = 'Isoprene',
                                    'BENZENE' = 'Benzene',
                                    'TOLUENE' = 'Toluene',
                                    'MXYL' = 'm-Xylene',
                                    'PXYL' = 'p-Xylene',
                                    'OXYL' = 'o-Xylene',
                                    'EBENZ' = 'Ethylbenzene') `
    );
    
    #function for labeller
    $R->run(q` facet_labeller = function(variable,value) {return(facet_labels[value])} `); 
    $R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);
    $R->run(q` plot.lines = function () { list( geom_line(size = 3), 
                                                geom_point(size = 5), 
                                                scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, 1)), expand_limits(x = 1, y = 0), 
                                                facet_grid(. ~ NMVOC, scales = "free", space = "free", labeller = facet_labeller), 
                                                scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 1)), 
                                                theme_bw(), 
                                                theme(axis.text.x = element_text(size = 30)), 
                                                theme(axis.text.y = element_text(size = 30)), 
                                                theme(legend.title = element_text(size = 50, face = "bold")), 
                                                theme(legend.text = theme_text(size = 40)), 
                                                theme(legend.key.size = unit(4, "cm")), 
                                                theme(strip.background = element_blank()), 
                                                theme(strip.text.x = element_text(size = 40, face = "bold", vjust = 1)), 
                                                theme(axis.title.x = element_blank()), 
                                                theme(axis.title.y = element_blank()), 
                                                scale_colour_manual(values = my.colours), 
                                                theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) ,
                                                theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) ,
                                                theme(axis.ticks.length = unit(.3, "cm")),
                                                theme(axis.ticks.margin = unit(.2, "cm")),
                                                theme(axis.ticks = element_line(size = 2)) ,
                                                theme(panel.margin=unit(0.2 , "lines")), 
                                                theme(legend.key = element_blank()) ) } `,
    );
    
    my $plot = $R->run( q` plot1 = ggplot(data = subset(data, NMVOC == "C2H6" | NMVOC == "C3H8" | NMVOC == "NC4H10" | NMVOC == "IC4H10" | NMVOC == "NC5H12"), aes(x = times, y = TOPP.value, colour = Mechanism)) `,
                        q` plot1 = plot1 + plot.lines() `, 

                        q` plot2 = ggplot(data = subset(data, NMVOC == "IC5H12" | NMVOC == "NC6H14" | NMVOC == "NC7H16" | NMVOC == "NC8H18" | NMVOC == "C2H4"), aes(x = times, y = TOPP.value, colour = Mechanism)) `,
                        q` plot2 = plot2 + plot.lines() `, 

                        q` plot3 = ggplot(data = subset(data, NMVOC == "C3H6" | NMVOC == "BUT1ENE" | NMVOC == "MEPROPENE" | NMVOC == "C5H8" | NMVOC == "BENZENE"), aes(x = times, y = TOPP.value, colour = Mechanism)) `,
                        q` plot3 = plot3 + plot.lines() `, 

                        q` plot4 = ggplot(data = subset(data, NMVOC == "TOLUENE" | NMVOC == "MXYL" | NMVOC == "OXYL" | NMVOC == "PXYL" | NMVOC == "EBENZ"), aes(x = times, y = TOPP.value, colour = Mechanism)) `,
                        q` plot4 = plot4 + plot.lines() `, 

                        q` y.label = textGrob("TOPP (molecules(Ox)/molecules(VOC))", rot = 90, gp = gpar(fontsize = 50, fontface = "bold"), vjust = 0.5) `,
                        q` legend = gtable_filter(ggplot_gtable(ggplot_build(plot1)), "guide-box") `,
                        q` CairoPDF(file = "TOPP_sums_all_species.pdf", width = 30, height = 20) `,
                        q` main.plot = grid.arrange(y.label, 
                                                    arrangeGrob(plot1 + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                                plot2 + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                                plot3 + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                                plot4 + theme(legend.position="none"), 
                                                                nrow = 4), 
                                                    sub = textGrob("Time (days)\n", gp = gpar(fontsize = 50, fontface = "bold"), vjust = 0.5), 
                                                    legend, 
                                                    nrow = 1, ncol = 3,
                                                    widths=unit.c(unit(4, "lines"), unit(1, "npc") - unit(4, "lines") - legend$width, legend$width)) `,
                        q` print(main.plot) `,
                        q` dev.off() `,
    );

    $R->stop();
}


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
            $TOPP{$species} = \@TOPP;
        }
        close FILE;
    }
}
