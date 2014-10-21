#! /usr/bin/env perl
# compare NO source between all mechanisms
# Version 0: Jane Coates 29/9/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
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
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);

    my $reaction = $kpp->producing_from("NO", "UNITY");
    my $reaction_number = $kpp->reaction_number(@$reaction);
    my $rate = $mecca->rate($reaction_number);
    $plot_data{$mechanisms[$index]} = $rate(1:$NTIME-2);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('time', [map { $_ } $time->dog]);
$R->run(q` data = data.frame(time) `);

foreach my $run (sort keys %plot_data) {
    $R->set('mechanism', $run);
    $R->set('rate', [map { $_ } $plot_data{$run}->dog]);
    $R->run(q` data[mechanism] = rate `);
}

$R->run(q` data = melt(data, id.vars = "time", variable.name = "Mechanism", value.name = "Rate") `);
$R->run(q` my.colours = c(  "CB05" = "#0352cb",
                            "CBM-IV" = "#b569b3",
                            "CRIv2" = "#ef6638", 
                            "MCMv3.1" = "#000000",
                            "MCMv3.2" = "#dc3522",
                            "MOZART-4" = "#cc9900",
                            "RACM" = "#6c254f",
                            "RACM2" = "#4682b4",
                            "RADM2" = "#035c28") `,);
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` plot = ggplot(data, aes(x = time, y = Rate, colour = Mechanism)) `,
        q` plot = plot + geom_line(size = 3) `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab(expression(bold(paste("Reaction Rate (molecules ", cm^-3, " ", s^-1, ")")))) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 85, face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 85, face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 60)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 60)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(4, "cm")) `,
        q` plot = plot + theme(legend.text = element_text(size = 40)) `,
        q` plot = plot + theme(legend.justification = c(0.99, 0.99)) `,
        q` plot = plot + theme(legend.position = c(0.99, 0.99)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "NO_source_all_mechanisms.pdf", width = 57, height = 40) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
