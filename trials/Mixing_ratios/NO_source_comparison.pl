#! /usr/bin/env perl
# compare NO source between all mechanisms
# Version 0: Jane Coates 29/9/2014
# Version 1: Jane Coates 5/12/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$NTIME-2);

my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my $index = 0;

my %plot_data;
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);

    my $reaction = $kpp->producing_from("NO", "UNITY");
    my $reaction_number = $kpp->reaction_number(@$reaction);
    my $rate = $mecca->rate($reaction_number);
    $plot_data{$mechanisms[$index]} = $rate(1:$NTIME-2);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('time', [map { $_ } $time->dog]);
$R->run(q` data = data.frame(time) `);

foreach my $run (sort keys %plot_data) {
    $R->set('mechanism', $run);
    $R->set('rate', [map { $_ } $plot_data{$run}->dog]);
    $R->run(q` data[mechanism] = rate `);
}

$R->run(q` data = gather(data, Mechanism, Rate, -time) `);
$R->run(q` my.colours = c(  "CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `,);

$R->run(q` plot = ggplot(data, aes(x = time, y = Rate, colour = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "NO_source_all_mechanisms.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
