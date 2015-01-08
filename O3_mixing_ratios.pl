#! /usr/bin/env perl
# Compare O3 mixing ratio time series
# Version 0: Jane Coates 18/12/2014
# Version 1: Jane Coates 8/1/2015 Adding OH to plot

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$NTIME-2);

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" ); 
my @species = qw(O3 OH);
my %mixing_ratio;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    foreach my $species (@species) {
        my $mixing_ratio = $mecca->tracer($species);
        $mixing_ratio{$mechanism}{$species} = $mixing_ratio(1:$NTIME-2) * 1e9;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame() `); 
foreach my $mechanism (sort keys %mixing_ratio) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `);
    foreach my $species (sort keys %{$mixing_ratio{$mechanism}}) {
        $R->set('species', $species);
        $R->set('mixing.ratio', [map { $_ } $mixing_ratio{$mechanism}{$species}->dog]);
        $R->run(q` pre[species] = mixing.ratio `);
    }
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Species, Mixing.Ratio, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#ef6638", "CRIv2" = "#b569b3", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Species, scales = "free" ) `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
);

$R->run(q` CairoPDF(file = "O3_mixing_ratios.pdf", width = 8.0, height = 5.6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
