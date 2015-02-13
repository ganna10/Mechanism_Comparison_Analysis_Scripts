#!/usr/bin/env perl
# Compare O3 concentration time series when with and without permutation reactions
# Version 0: Jane Coates 5/12/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @runs = qw( MOZART-4 RADM2 RACM RACM2 CBM-IV CB05);
my %plot_data;

foreach my $mechanism (@runs) {
    my $orig_base = "$base/${mechanism}_no_RO2";
    my $orig_boxmodel = "$orig_base/boxmodel";
    my $orig_mecca = MECCA->new($orig_boxmodel);
    my $orig_O3_mixing_ratio = $orig_mecca->tracer('O3');
    $plot_data{"Original_$mechanism"} = $orig_O3_mixing_ratio;

    my $perm_base = "$base/${mechanism}_perm";
    my $perm_boxmodel = "$perm_base/boxmodel";
    my $perm_mecca = MECCA->new($perm_boxmodel);
    my $perm_O3_mixing_ratio = $perm_mecca->tracer('O3');
    $plot_data{"Permutation_$mechanism"} = $perm_O3_mixing_ratio;
}

my $mecca = MECCA->new("/local/home/coates/MECCA/MOZART-4_perm/boxmodel");
my $times = $mecca->time;
my $NTIME = $times->nelem;
$times -= $times->at(0);
$times /= 86400;
my $time = pdl_to_array($times);

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(tidyr) `);
$R->run(q` library(Cairo) `);
$R->run(q` library(grid) `);

$R->set('time', [@$time]);
$R->run(q` mozart.data = data.frame(time) `);
$R->run(q` radm2.data = data.frame(time) `);
$R->run(q` racm.data = data.frame(time) `);
$R->run(q` racm2.data = data.frame(time) `);
$R->run(q` cbm4.data = data.frame(time) `);
$R->run(q` cb05.data = data.frame(time) `);

foreach my $run (sort keys %plot_data) {
    $plot_data{$run} = pdl_to_array($plot_data{$run}) ; 
    if ($run =~ /Original/) {
        $R->set('name', 'Original');
    } elsif ($run =~ /Permutation/) {
        $R->set('name', 'Modified');
    }
    $R->set('conc', [@{$plot_data{$run}}]); #convert to ppb
    if ($run =~ /MOZART/) {
        $R->run(q` mozart.data[name] = conc * 1e9 `);
    } elsif ($run =~ /RADM2/) {
        $R->run(q` radm2.data[name] = conc * 1e9 `);
    } elsif ($run =~ /RACM2/) {
        $R->run(q` racm2.data[name] = conc * 1e9 `);
    } elsif ($run =~ /RACM/) {
        $R->run(q` racm.data[name] = conc * 1e9 `);
    } elsif ($run =~ /CBM-IV/) {
        $R->run(q` cbm4.data[name] = conc * 1e9 `);
    } elsif ($run =~ /CB05/) {
        $R->run(q` cb05.data[name] = conc * 1e9 `);
    }
}

$R->run(q` arrange.data = function (data, mechanism) { data = gather(data, Run, Concentration, -time);
                                                       data["Mechanism"] = rep(mechanism, length(time));
                                                       return(data) } `,
);

$R->run(q` mozart.data = arrange.data(mozart.data, "MOZART-4") `,
        q` radm2.data = arrange.data(radm2.data, "RADM2") `,
        q` racm.data = arrange.data(racm.data, "RACM") `,
        q` racm2.data = arrange.data(racm2.data, "RACM2") `,
        q` cbm4.data = arrange.data(cbm4.data, "CBM-IV") `,
        q` cb05.data = arrange.data(cb05.data, "CB05") `,
        q` data = rbind(mozart.data, radm2.data, racm.data, racm2.data, cbm4.data, cb05.data) `,
        q` data$Mechanism = factor(data$Mechanism, levels = c("RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `,
);

#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` plot = ggplot(data = data, aes(x = time, y = Concentration, colour = Run)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 1)) `,
        #q` plot = plot + scale_y_continuous(limits = c(30, 350), breaks = seq(30, 350, 50)) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("O3 Mixing Ratio (ppbv)") `, 
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `, 
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(1.01, 1.01)) `,
        q` plot = plot + theme(legend.justification = c(1.01, 1.01)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_colour_manual(values = c("#6c254f", "#4c9383")) `,

        q` CairoPDF(file = "O3_mixing_ratio_comparison.pdf", width = 6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub pdl_to_array {
    my ($pdl) = @_;
    $pdl = $pdl(1:$NTIME-2);
    my @array = map { $_ } $pdl->dog;
    return \@array;
}
