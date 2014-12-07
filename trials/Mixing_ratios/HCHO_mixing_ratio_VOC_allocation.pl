#! /usr/bin/env perl
# Allocation HCHO mixing ratio time series in each mechanism to parent VOCs
# Version 0: Jane Coates 18/11/2014
# Version 1: Jane Coates 5/12/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;

#my @runs = qw( MOZART_tagging );
#my @mechanisms = qw( MOZART-4 );
#my @species = qw( CH2O );
my @mechanisms = ( "(a) MCM v3.2", "(b) MCM v3.1", "(c) CRI v2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
my @species = qw( HCHO HCHO HCHO CH2O HCHO HCHO HCHO HCHO FORM );
my $index = 0;
my %plot_data;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $spc_file = "$base/${mechanism}_tagged/gas.spc";
    my $all_tagged_species = get_tagged_species($species[$index], $spc_file); 
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $all_tagged_species);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

my $times = $mecca->time;
$times = $times(1:$NTIME-2);
$times -= $times->at(0);
$times /= 86400 ;
$R->set('Time', [ map { $_ } $times->dog]);
$R->run(q` data = data.frame() `);

foreach my $run (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $VOC (sort keys %$ref) {
            my ($name) = get_name($VOC);
            $R->set('voc', $name);
            $R->set('mixing.ratio', [ map { $_ } $ref->{$VOC}->dog ]);
            $R->run(q` pre[voc] = mixing.ratio `);
        }
    }
    $R->set('mechanism', $run);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "VOC", value.name = "Mixing.Ratio") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c(  "CO" = "#2b9eb3" ,
                            "NO2" = "#b569b3" ,
                            "Methane" = "#1b695b" ,
                            "2-Methylbutane" = "#ae4901" ,
                            "Propane" = "#e7e85e" ,
                            "Toluene" = "#0e5c28" ,
                            "Butane" = "#f3aa7f" ,
                            "Ethene" = "#898989" ,
                            "O3" = "#1c3e3d" ,
                            "Pentane" = "#f9c500" ,
                            "Isoprene" = "#8c1531" ,
                            "Ethane" = "#86b650" ,
                            "m-Xylene" = "#ef6638" ,
                            "Propene" = "#0352cb" ,
                            "NO" = "#c9a415" ,
                            "Hexane" = "#9bb18d" ,
                            "2-Methylpropane" = "#a67c52" ,
                            "2-Methylpropene" = "#dc3522" ,
                            "Benzene" = "#f7c56c" ,
                            "o-Xylene" = "#4c9383" ,
                            "Ethylbenzene" = "#ba8b01" ,
                            "Others" = "#58691b" ) `,
                            q` data$VOC = factor(data$VOC, levels = c("Methane", "Pentane", "Ethene", "Propane", "Toluene", "Propene", "Isoprene", "Ethane", "2-Methylbutane", "Butane", "m-Xylene", "2-Methylpropane", "2-Methylpropene", "Others")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, fill = VOC)) `,
        q` plot = plot + geom_area(position = "stack") `,
        q` plot = plot + geom_line(position = "stack", colour = "black", show_guide = FALSE) `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("\nMixing Ratio (ppbv)\n") `,
        q` plot = plot + xlab("\nTime (days)\n") `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + theme(axis.title = element_text(size = 180, face = "bold")) `,
        q` plot = plot + theme(axis.text = element_text(size = 140)) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 140)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold", vjust = 1)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$VOC))) `,
);

$R->run(q` CairoPDF(file = "HCHO_VOC_allocation_mixing_ratio.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $tagged_species) = @_;
    
    my %mixing_ratios;
    foreach my $species (@$tagged_species) { 
        my $mixing_ratio = $mecca->tracer($species);
        next if ($mixing_ratio->sum == 0);
        my ($HCHO, $parent) = split /_/, $species;
        $mixing_ratios{$parent} += $mixing_ratio(1:$NTIME-2) * 1e9; #convert to ppbv
    }
    
    my $others_max = 20;
    foreach my $VOC (keys %mixing_ratios) {
        if ($mixing_ratios{$VOC}->sum < $others_max) {
            $mixing_ratios{"Others"} += $mixing_ratios{$VOC};
            delete $mixing_ratios{$VOC};
        }
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($mixing_ratios{$b}) <=> &$sort_function($mixing_ratios{$a}) } keys %mixing_ratios;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Others') ;
        push @final_sorted_data, { $_ => $mixing_ratios{$_} };
    } 
    push @final_sorted_data, { 'Others' => $mixing_ratios{'Others'} } if (defined $mixing_ratios{'Others'}); 

    return \@final_sorted_data;
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

sub read_file {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file for reading : $!";
    chomp(my @all = <$in>);
    close $in;
    return \@all;
} 

sub get_name {
    my ($parent) = @_;
    if ($parent eq "CH4") {
        $parent = "Methane";
    } elsif ($parent =~ /C2H6|ETH/) {
        $parent = "Ethane";
    } elsif ($parent =~ /C3H8|HC3/) {
        $parent = "Propane";
    } elsif ($parent eq "NC4H10") {
        $parent = "Butane";
    } elsif ($parent eq "IC4H10") {
        $parent = "2-Methylpropane";
    } elsif ($parent =~ /NC5H12|BIGALK|HC5/) {
        $parent = "Pentane";
    } elsif ($parent eq "IC5H12") {
        $parent = "2-Methylbutane";
    } elsif ($parent eq "NC6H14") {
        $parent = "Hexane";
    } elsif ($parent eq "NC7H16") {
        $parent = "Heptane";
    } elsif ($parent =~ /NC8H18|HC8/) {
        $parent = "Octane";
    } elsif ($parent =~ /C2H4|OL2|ETE/) {
        $parent = "Ethene";
    } elsif ($parent =~ /C3H6|OLT/) {
        $parent = "Propene";
    } elsif ($parent =~ /BUT1ENE|BIGENE/) {
        $parent = "Butene";
    } elsif ($parent =~ /MEPROPENE|OLI/) {
        $parent = "2-Methylpropene";
    } elsif ($parent =~ /C5H8|ISO/) {
        $parent = "Isoprene";
    } elsif ($parent =~ /^BEN/) {
        $parent = "Benzene";
    } elsif ($parent =~ /TOL/) {
        $parent = "Toluene";
    } elsif ($parent =~ /OXYL|XYO/) {
        $parent = "o-Xylene";
    } elsif ($parent =~ /PXYL|XYP/) {
        $parent = "p-Xylene";
    } elsif ($parent eq "EBENZ") {
        $parent = "Ethylbenzene";
    } elsif ($parent =~ /MXYL|XYL|XYM/) {
        $parent = "m-Xylene";
    } elsif ($parent =~ /Others/) {
        $parent = $parent;
    } else {
        print "No chemical name for $parent\n";
    }
    return $parent;
}
