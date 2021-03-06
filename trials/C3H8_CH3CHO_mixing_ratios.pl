#! /usr/bin/env perl
# Compare CH3CHO mixing ratio time series
# Version 0: Jane Coates 9/7/2015

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

my @mechanisms = ( "MCMv3.2", "RADM2", "RACM", "RACM2" );
my $index = 0;

my %mixing_ratio;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $species;
    if ($mechanism =~ /MCM/) {
        $species = "CH3CHO_C3H8";
    } elsif ($mechanism =~ /RAD/ or $mechanism eq "RACM") {
        $species = "ALD_HC3";
    } elsif ($mechanism eq "RACM2") {
        $species = "ACD_HC3";
    } 
    my $mixing_ratio = $mecca->tracer($species);
    $mixing_ratio{$mechanisms[$index]} = $mixing_ratio(1:$NTIME-2) * 1e9;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame(Time) `);

foreach my $mechanism (sort keys %mixing_ratio) {
    $R->set('mechanism', $mechanism);
    $R->set('mixing.ratio', [map { $_ } $mixing_ratio{$mechanism}->dog]);
    $R->run(q` data[mechanism] = mixing.ratio `);
}
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

$R->run(q` data = gather(data, Mechanism, Mixing.Ratio, -Time) `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0.5)) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
);

$R->run(q` CairoPDF(file = "C3H8_ALD_mixing_ratios.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
