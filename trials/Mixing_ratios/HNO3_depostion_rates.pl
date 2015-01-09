#! /usr/bin/env perl
# HNO3 depostion rates in each mechanism
# Version 0: Jane Coates 9/1/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400 ;

my @mechanisms = qw( MCMv3.2 MCMv3.2 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my %data;
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn);
    my $consumers = $kpp->consuming("HNO3");
    my $reaction;
    foreach (@$consumers) {
        my $reaction_string = $kpp->reaction_string($_);
        next unless ($reaction_string =~ /=\sUNITY/);
        print "$reaction_string\n";
        $reaction = $_;
    }
    my $reaction_number = $kpp->reaction_number($reaction);
    my $rate = $mecca->rate($reaction_number);
    $data{$mechanism} = $rate;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
);

$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` data = data.frame(Time) `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->set('rate', [ map { $_ } $data{$mechanism}->dog ]);
    $R->run(q` data[mechanism] = rate `);
}
$R->run(q` data = gather(data, Mechanism, Rate, -Time) `);
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HNO3_deposition_rates.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
