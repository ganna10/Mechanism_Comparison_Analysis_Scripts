#! /usr/bin/env perl
# Compare % of IC5H12 mixing ratios at beginning of 2nd day in MCM 3.2, CBM-IV and CB05
# Version 0: Jane Coates 3/3/2015

use strict;
use diagnostics;
use PDL;
use MECCA;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @mechanisms = qw( MCMv3.2 CBM-IV CB05 );
my %data;

my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $species;
    if ($mechanism =~ /MCM/) {
        $species = "IC5H12";
    } else {
        $species = "PAR_IC5H12";
    }
    my $mr = $mecca->tracer($species);
    print "$mechanism at day 1: ", $mr->at(1), "\n";
    print "$mechanism at day 2: ", $mr->at($n_per_day), "\n";
    $mr *= 5 if ($mechanism =~ /MCM/);
    $data{$mechanism} = $mr(1:$ntime-2) * 1e9;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` data = data.frame(Time) `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}->dog ]);
    $R->run(q` data[mechanism] = mixing.ratio `);
}
$R->run(q` data = gather(data, Mechanism, Mixing.Ratio, -Time) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
);

$R->run(q` CairoPDF(file = "IC5H12_mixing_ratios_MCM_CBs.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
