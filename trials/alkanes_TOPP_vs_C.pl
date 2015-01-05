#! /usr/bin/env perl
# Correlate first day TOPP and final cumulative TOPP with C number of parent VOC
# Version 0: Jane Coates 8/12/2014

use strict;
use diagnostics;
use Statistics::R;

my $dir = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper/TOPP_plots";
opendir DIR, $dir or die "Can't open $dir : $!";
my @files = grep { $_ =~ /_TOPP_values/ } readdir DIR;
closedir DIR;

my (%TOPP_first_day, %TOPP_cumulative);
foreach my $file (@files) {
    my @lines = split /\n/, read_file("$dir/$file");
    (my $mechanism = $file) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next unless ($line =~ /^NC|^IC|^C2H6|^C3H8|^HC|^BIGALK|^ETH/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $TOPP_first_day{$mechanism}{$VOC} = $TOPPs[0];
        my $sum = 0;
        $sum += $_ foreach (@TOPPs);
        $TOPP_cumulative{$mechanism}{$VOC} = $sum;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(grid) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(plyr) `,
        q` library(Cairo) `,
);

$R->run(q` first.day.data = data.frame(Mechanism = as.numeric(0), VOC = as.numeric(0), C.number = as.numeric(0), TOPP = as.numeric(0)) `,
        q` cumulative.data = data.frame(Mechanism = as.numeric(0), VOC = as.numeric(0), C.number = as.numeric(0), TOPP = as.numeric(0)) `,
);
foreach my $mechanism (sort keys %TOPP_first_day) {
    $R->set('mechanism', $mechanism);
    foreach my $VOC (sort keys %{$TOPP_first_day{$mechanism}}) {
        my $C_number = get_C_number($VOC);
        my $name = get_chemical_name($VOC);
        $R->set('c.number', $C_number);
        $R->set('voc', $name);
        $R->set('first.day', $TOPP_first_day{$mechanism}{$VOC});
        $R->set('cumulative', $TOPP_cumulative{$mechanism}{$VOC});
        $R->run(q` cumulative.data = rbind(cumulative.data, c(mechanism, voc, c.number, cumulative)) `);
        $R->run(q` first.day.data = rbind(first.day.data, c(mechanism, voc, c.number, first.day)) `);
    }
}

$R->run(q` my.colours = c( "Ethane" = "#696537", "Propane" = "#f9c600", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb", "Ethene" = "#86b650", "Propene" = "#6c254f", "Butene" = "#ee6738", "2-Methylpropene" = "#58691b", "Isoprene" = "#8ed6d5", "Benzene" = "#f3aa7f", "Toluene" = "#c65d6c", "m-Xylene" = "#888a87", "o-Xylene" = "#0e5c28", "p-Xylene" = "#b569b3", "Ethylbenzene" = "#2c9def" ) `);

$R->run(q` lm_eqn = function(df){   m = lm(TOPP ~ C.number, df);
                                    eq <- substitute(italic(TOPP) == a + b %.% italic(C.number)*","~~italic(r)^2~"="~r2, 
                                                        list(a = format(coef(m)[1], digits = 4), 
                                                                b = format(coef(m)[2], digits = 4), 
                                                                r2 = format(summary(m)$r.squared, digits = 5)));
                                    as.character(as.expression(eq)) } `);
                                                    
$R->run(q` plotting = function (data, filename) {   eq = ddply(data, .(Mechanism), lm_eqn);
                                                    plot = ggplot(data, aes(x = C.number, y = TOPP, colour = VOC)) ;
                                                    plot = plot + geom_point(shape = 19) ;
                                                    plot = plot + theme_bw();
                                                    plot = plot + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black");
                                                    plot = plot + geom_text(data = eq, aes(x = 2.5, y = 1.5, label = V1), size = 1.3, parse = TRUE);
                                                    plot = plot + facet_wrap( ~ Mechanism);
                                                    plot = plot + xlab("Carbon Number of VOC");
                                                    plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC))");
                                                    plot = plot + theme(legend.title = element_blank());
                                                    plot = plot + theme(axis.title = element_text(face = "bold"));
                                                    plot = plot + theme(strip.text = element_text(face = "bold"));
                                                    plot = plot + theme(strip.background = element_blank());
                                                    plot = plot + theme(panel.grid = element_blank());
                                                    plot = plot + theme(panel.border = element_rect(colour = "black"));
                                                    plot = plot + theme(legend.key = element_blank());
                                                    plot = plot + scale_colour_manual(values = my.colours);
                                                    CairoPDF(file = filename, width = 8, height = 8);
                                                    print(plot);
                                                    dev.off() } `);
    
$R->run(q` first.day.data = first.day.data[-1,] `,
        q` plotting(first.day.data, "Alkanes_first_day_TOPP_vs_C_number.pdf") `,
        q` cumulative.data = cumulative.data[-1,] `,
        #q` plotting(cumulative.data, "Alkanes_cumulative_TOPP_vs_C_number.pdf") `,
);
my $p = $R->run(q` print(first.day.data) `);
print $p, "\n";

$R->stop();

sub read_file {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    close $in;
    return $data;
}

sub get_C_number {
    my ($VOC) = @_;
    my $C;
    if ($VOC eq "C2H6" or $VOC eq "ETH" or $VOC eq "C2H4" or $VOC eq "OL2" or $VOC eq "ETE") {
        $C = 2;
    } elsif ($VOC eq "C3H8" or $VOC eq "HC3" or $VOC eq "C3H6" or $VOC eq "OLT") {
        $C = 3;
    } elsif ($VOC eq "NC4H10" or $VOC eq "IC4H10" or $VOC eq "BUT1ENE" or $VOC eq "MEPROPENE" or $VOC eq "BIGENE" or $VOC eq "OLI") {
        $C = 4;
    } elsif ($VOC eq "NC5H12" or $VOC eq "IC5H12" or $VOC eq "C5H8" or $VOC =~ /ISO/ or $VOC eq "BIGALK" or $VOC eq "HC5") {
        $C = 5;
    } elsif ($VOC eq "NC6H14" or $VOC eq "BENZENE" or $VOC eq "BEN") {
        $C = 6;
    } elsif ($VOC eq "NC7H16" or $VOC eq "TOLUENE" or $VOC eq "TOL") {
        $C = 7;
    } elsif ($VOC eq "NC8H18" or $VOC eq "HC8" or $VOC =~ /XY/ or $VOC eq "EBENZ") {
        $C = 8;
    } else {
        print "No C for $VOC\n";
    }
    return $C;
}

sub get_functionality {
    my ($VOC) = @_;
    my $group;
    if ($VOC eq "C2H6" or $VOC eq "C3H8" or $VOC eq "ETH" or $VOC eq "BIGALK" or $VOC =~ /^NC|^IC|^HC/) {
        $group = "Alkane";
    } elsif ($VOC eq "C2H4" or $VOC eq "C3H6" or $VOC =~ /^OL|^ISO/ or $VOC eq "MEPROPENE" or $VOC eq "BUT1ENE" or $VOC eq "BIGENE" or $VOC eq "ETE" or $VOC eq "C5H8") {
        $group = "Alkene";
    } elsif ($VOC =~ /^BEN|^TOL|XY/ or $VOC eq "EBENZ") {
        $group = "Aromatic";
    } else {
        print "No group for $VOC\n";
    }
    return $group;
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
