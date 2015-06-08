#! /usr/bin/env perl
# Compare O3 mixing ratio time series
# Version 0: Jane Coates 18/12/2014
# Version 1: Jane Coates 8/1/2015 Adding OH to plot
# Version 2: Jane Coates 6/2/2015 Removing OH from plot, calculating difference on first day
# Version 3: Jane Coates 5/3/2015 calculating difference of first day O3 in RACM2

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
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" ); 
my %mixing_ratio;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $mixing_ratio = $mecca->tracer("O3");
    $mixing_ratio{$mechanism}{"O3"} = $mixing_ratio(1:$NTIME-2) * 1e9;
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
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0.5)) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("O3 Mixing Ratio (ppbv)") `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.position = c(1, 1.058)) `,
        q` plot = plot + theme(legend.justification = c(1, 1.058)) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
);

#$R->run(q` CairoPDF(file = "O3_mixing_ratios.pdf", width = 8.0, height = 5.6) `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

$R->stop();

#calculate difference on first day between RADM2 and RACM (largest and lowest O3 mixing ratios)
my $radm2_O3 = $mixing_ratio{"RADM2"}{"O3"};
my $day2 = 2 * $n_per_day;
my $radm2_O3_day2 = $radm2_O3(0:$day2);
my $radm2_day2_max = $radm2_O3_day2->max;
print "RADM2: $radm2_day2_max\n";

my $mcm_O3 = $mixing_ratio{"MCMv3.2"}{"O3"};
my $mcm_O3_day2 = $mcm_O3(0:$day2);
my $mcm_day2_max = $mcm_O3_day2->max;
print "MCM v3.2 $mcm_day2_max\n";
print "MCM - RADM2 ", $mcm_day2_max - $radm2_day2_max, "\n";

my $racm2_O3 = $mixing_ratio{"RACM2"}{"O3"};
my $racm2_O3_day2 = $racm2_O3(0:$day2);
my $racm2_day2_max = $racm2_O3_day2->max;
print "RACM2: $racm2_day2_max\n";
print "MCM - RACM2 ", $mcm_day2_max - $racm2_day2_max, "\n";
