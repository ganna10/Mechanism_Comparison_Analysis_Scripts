#! /usr/bin/env perl
# Compare H2O2 mixing ratio time series
# Version 0: Jane Coates 26/11/2014
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
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$NTIME-2);

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my $index = 0;
my %mixing_ratio;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $mixing_ratio = $mecca->tracer("H2O2");
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
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "H2O2_mixing_ratios.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
