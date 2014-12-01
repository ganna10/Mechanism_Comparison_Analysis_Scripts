#! /usr/bin/env perl
# Compare OH mixing ratio time series
# Version 0: Jane Coates 18/11/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$NTIME-2);

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "MCM v3.2", "MCM v3.1", "CRI v2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my $index = 0;

my %mixing_ratio;

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $mixing_ratio = $mecca->tracer("OH");
    $mixing_ratio{$mechanisms[$index]} = $mixing_ratio(1:$NTIME-2) * 1e9;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame(Time) `);

foreach my $mechanism (sort keys %mixing_ratio) {
    $R->set('mechanism', $mechanism);
    $R->set('mixing.ratio', [map { $_ } $mixing_ratio{$mechanism}->dog]);
    $R->run(q` data[mechanism] = mixing.ratio `);
}

$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Mechanism", value.name = "Mixing.Ratio") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
);

$R->run(q` CairoPDF(file = "OH_mixing_ratios.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
