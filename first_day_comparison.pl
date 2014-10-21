#! /usr/bin/perl
# compare first day TOPP values of all VOCs from all mechanisms to MCM v3.2 TOPP values in a separate facet including y = x line
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
    my ($name, $mechanism) = split /_/, $species;
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

#reshape data frame and manually specify order of NMVOC
$R->run(q` casted = cast(data, NMVOC ~ Mechanism, value = 'TOPP.value') `);
$R->run(q` casted$NMVOC <- factor(casted$NMVOC, levels = c("C2H6", "C3H8", "NC4H10", "IC4H10", "NC5H12", "IC5H12", "NC6H14", "NC7H16", "NC8H18", "C2H4", "C3H6", "BUT1ENE", "MEPROPENE", "C5H8", "BENZENE", "TOLUENE", "MXYL", "OXYL", "PXYL", "EBENZ"), ordered = TRUE) `);

#specifiying colours and names
$R->run(q` my.colours = c( "C2H6" = "#696537", "C3H8" = "#f9c600", "NC4H10" = "#76afca", "IC4H10" = "#dc3522", "NC5H12" = "#8c6238", "IC5H12" = "#9bb08f", "NC6H14" = "#8b1537", "NC7H16" = "#ba8b01", "NC8H18" = "#0352cb", "C2H4" = "#86b650", "C3H6" = "#6c254f", "BUT1ENE" = "#ee6738", "MEPROPENE" = "#58691b", "C5H8" = "#8ed6d5", "BENZENE" = "#f3aa7f", "TOLUENE" = "#c65d6c", "MXYL" = "#888a87", "OXYL" = "#0e5c28", "PXYL" = "#b569b3", "EBENZ" = "#2c9def" ) `,
        q` my.names = c("C2H6" = "Ethane", "C3H8" = "Propane", "NC4H10" = "Butane", "IC4H10" = "2-Methylpropane", "NC5H12" = "Pentane", "IC5H12" = "2-Methylbutane", "NC6H14" = "Hexane", "NC7H16" = "Heptane", "NC8H18" = "Octane", "C2H4" = "Ethene", "C3H6" = "Propene", "BUT1ENE" = "Butene", "MEPROPENE" = "2-Methylpropene", "C5H8" = "Isoprene", "BENZENE" = "Benzene", "TOLUENE" = "Toluene", "MXYL" = "m-Xylene", "OXYL" = "o-Xylene", "PXYL" = "p-Xylene", "EBENZ" = "Ethylbenzene") `,
);

$R->run(q` plot.lines = function () { list( geom_point(shape = 19, size = 55), 
                                            xlab("\nMCM v3.2 TOPP (molecules(Ox)/molecules(VOC))"), 
                                            scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), 
                                            scale_y_continuous(limits = c(0,7), breaks = seq(0, 7, 1)), 
                                            geom_abline(intercept = 0, slope = 1, size = 3), 
                                            theme_bw(), 
                                            theme(axis.text.x = element_text(size = 150), 
                                            axis.text.y = element_text(size = 150), 
                                            axis.title.x = element_text(size = 180, face = "bold"), 
                                            axis.title.y = element_blank()), 
                                            theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) ,
                                            theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) ,
                                            theme(axis.ticks.length = unit(.9, "cm")),
                                            theme(axis.ticks.margin = unit(.8, "cm")),
                                            theme(axis.ticks = element_line(size = 4)) ,
                                            scale_colour_manual(labels = my.names, values = my.colours, guide = guide_legend(label.position = "bottom", label.theme = element_text(size = 170, angle = 45), label.hjust = 0.5, label.vjust = 0.9)), 
                                            theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "cm"), legend.key = element_blank()) ) } `,
); 

$R->run(q` mcm3.1.plot = ggplot(data = casted, aes(x = MCMv3.2, y = MCMv3.1, colour = NMVOC)) `,
        q` mcm3.1.plot = mcm3.1.plot + plot.lines() `, 
        
        q` cri.plot = ggplot(data = casted, aes(x = MCMv3.2, y = CRIv2, colour = NMVOC)) `,
        q` cri.plot = cri.plot + plot.lines() `, 
        
        q` mozart.plot = ggplot(data = casted, aes(x = MCMv3.2, y = MOZART, colour = NMVOC)) `,
        q` mozart.plot = mozart.plot + plot.lines() `, 
        
        q` radm2.plot = ggplot(data = casted, aes(x = MCMv3.2, y = RADM2, colour = NMVOC)) `,
        q` radm2.plot = radm2.plot + plot.lines() `, 
        
        q` racm.plot = ggplot(data = casted, aes(x = MCMv3.2, y = RACM, colour = NMVOC)) `,
        q` racm.plot = racm.plot + plot.lines() `, 
        
        q` racm2.plot = ggplot(data = casted, aes(x = MCMv3.2, y = RACM2, colour = NMVOC)) `,
        q` racm2.plot = racm2.plot + plot.lines() `, 
        
        q` cbm4.plot = ggplot(data = casted, aes(x = MCMv3.2, y = CBM4, colour = NMVOC)) `,
        q` cbm4.plot = cbm4.plot + plot.lines() `, 
        
        q` cb05.plot = ggplot(data = casted, aes(x = MCMv3.2, y = CB05, colour = NMVOC)) `,
        q` cb05.plot = cb05.plot + plot.lines() `, 
        
        q` legend = gtable_filter(ggplot_gtable(ggplot_build(mcm3.1.plot)), "guide-box") `,
        q` y.label = textGrob("TOPP (molecules(Ox)/molecules(VOC))\n", rot = 90, gp = gpar(fontsize = 220, fontface = "bold"), vjust = 0.7) `,
        q` CairoPDF(file = "first_day_values.pdf", width = 180, height = 200) `,
        q` main.plot = grid.arrange(y.label, 
                                    arrangeGrob(mcm3.1.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(a) MCM v3.1", size = 100, face = "bold") + scale_x_continuous(breaks = NULL), 
                                                cri.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(b) CRI v2", size = 100, face = "bold") + theme(axis.ticks = element_blank()), 
                                                radm2.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(c) RADM2", size = 100, face = "bold") + scale_x_continuous(breaks = NULL), 
                                                racm.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(d) RACM", size = 100, face = "bold") + theme(axis.ticks = element_blank()), 
                                                racm2.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(e) RACM2", size = 100, face = "bold") + scale_x_continuous(breaks = NULL), 
                                                mozart.plot + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(f) MOZART-4", size = 100, face = "bold") + theme(axis.ticks = element_blank()), 
                                                cbm4.plot + theme(legend.position="none") + annotate("text", x = 1, y = 6.5, label = "(g) CBM-IV", size = 100, face = "bold"), 
                                                cb05.plot + theme(legend.position="none", axis.text.y = element_blank()) + annotate("text", x = 1, y = 6.5, label = "(h) CB05", size = 100, face = "bold") + theme(axis.ticks = element_blank()), 
                                                ncol = 2), 
                                    blank = grid.rect(gp = gpar(col = "white")), legend, 
                                    ncol=2,
                                    widths=unit.c(unit(26, "lines"), unit(1, "npc") - unit(26, "lines")),
                                    nrow=2,
                                    heights = c(175, 25)) `,

        q` print(main.plot) `, 
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
