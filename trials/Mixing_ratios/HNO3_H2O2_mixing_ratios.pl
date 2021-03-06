#! /usr/bin/env perl
# compare HNO3 concentrations in each mechanism
# Version 0: Jane Coates 12/10/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$NTIME-2);

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my $index = 0;

my %plot_data;
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $HNO3 = $mecca->tracer("HNO3");
    my $H2O2 = $mecca->tracer("H2O2");
    my $ratio = $HNO3 / $H2O2;
    $plot_data{$mechanisms[$index]} = $ratio(1:$NTIME-2) * 1e9;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('time', [map { $_ } $time->dog]);
$R->run(q` data = data.frame(time) `);
foreach my $run (sort keys %plot_data) {
    $R->set('mechanism', $run);
    $R->set('conc', [map { $_ } $plot_data{$run}->dog]);
    $R->run(q` data[mechanism] = conc `);
}

$R->run(q` data = melt(data, id.vars = c("time"), variable.name = "Mechanism", value.name = "Mixing.ratio") `);
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

$R->run(q` plot = ggplot(data, aes(x = time, y = Mixing.ratio, colour = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HNO3_to_H2O2_mixing_ratio.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
