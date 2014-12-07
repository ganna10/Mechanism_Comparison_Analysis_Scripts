#! /usr/bin/env perl
# Compare mixing ratio time series O3, OH and NO source
# Version 0: Jane Coates 5/12/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt ;
my $n_days = int $ntime / $n_per_day;

my @mechanisms = ("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my @species = qw(O3 OH NO-source);
my $index = 0;
my %mixing_ratio;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    foreach my $species (@species) {
        next if ($species eq "NO-source");
        my $mixing_ratio = $mecca->tracer($species);
        $mixing_ratio{$mechanisms[$index]}{$species} = $mixing_ratio(1:$ntime-2) * 1e9;
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

foreach my $mechanism (sort keys %mixing_ratio) {
    $R->set('mechanism', $mechanism);
    foreach my $species (sort keys %{$mixing_ratio{$mechanism}}) {
        $R->run(q` data = data.frame(Time) `);
        $R->set('species', $species);
        $R->set('mixing.ratio', [map { $_ } $mixing_ratio{$mechanism}{$species}->dog]);
        $R->run(q` data[mechanism] = mixing.ratio `);
    }
}
my $p = $R->run(q` print(data)`);
print $p, "\n";

#$R->run(q` data = gather(data, Mechanism, Mixing.Ratio, -Time) `);
#$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
#        q` plot = plot + geom_line() `,
#        q` plot = plot + scale_colour_manual(values = my.colours) `,
#        q` plot = plot + ggtitle("O3") `,
#        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
#        q` plot = plot + xlab("Time (Days)") `,
#        q` plot = plot + scale_y_continuous(limits = c(30, 350), breaks = seq(30, 350, 50)) `,
#        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
#        q` plot = plot + theme_bw() `,
#        q` plot = plot + theme(panel.grid = element_blank()) `,
#        q` plot = plot + theme(legend.title = element_blank()) `,
#        q` plot = plot + theme(legend.key = element_blank()) `,
#        q` plot = plot + theme(panel.border = element_rect(colour = "black", size = 0.5)) `,
#        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
#        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
#        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
#);
#
#$R->run(q` CairoPDF(file = "mixing_ratios.pdf") `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

$R->stop();
