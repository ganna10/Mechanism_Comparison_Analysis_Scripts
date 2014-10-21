#! /usr/bin/env perl
# Calculate total OPE of each mechanism's model run base on Kleinman:2002 paper (P(Ox)/P(HNO3))
# Version 0: Jane Coates 12/10/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use Statistics::R;
use PDL;
use PDL::NiceSlice;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
#my @runs = qw( RACM2_tagged CRI_tagging );
#my @mechanisms = qw( RACM2 CRIv2 );
my $index = 0;

my (%families, %weights, %plot_data);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    $plot_data{$mechanisms[$index]} = get_data($kpp, $mecca, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
);

my @days = ("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7");
$R->set('Time', [@days]);
$R->run(q` data = data.frame(Time) `);
foreach my $run (sort keys %plot_data) {
    $R->set('Mechanism', $run);
    $R->set('OPE', [map { $_ } $plot_data{$run}->dog]);
    $R->run(q` data[Mechanism] = OPE `);
}
$R->run(q` data = data[order(Time),] `,
        q` data = data[1:7,] `,
        q` data = melt(data, id.vars = c("Time"), variable.name = "Mechanism", value.name = "OPE") `,
);
my $p = $R->run(q` print(data) `);
print "$p\n";


$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = OPE)) `,
        q` plot = plot + geom_bar(stat = "identity") `, 
        q` plot = plot + facet_grid(Time ~ .) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05"))) `,
        q` plot = plot + ylab("\nOx Production Efficiency (molecules (Ox) / molecules (HNO3))\n") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text.y = element_text(size = 160, face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 160, face = "bold")) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
);

$R->run(q` CairoPDF(file = "Total_OPE.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism) = @_;
    my @loop = ("Ox_$mechanism", "HNO3");

    my ($producers, $producer_yields, %production_rates);
    foreach my $species (@loop) {
        if (exists $families{$species}) {
            $kpp->family({
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
        } else {
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
        }
        print "No producers found for $species\n" if (@$producers == 0);
        
        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number);
            next if ($rate->sum == 0);
            $production_rates{$species} += $rate(1:$NTIME-2);
        } 
    }

    my $dt = $mecca->dt->at(0); #model time step
    my $n_per_day = 43200 / $dt;
    my $n_days = int ($NTIME / $n_per_day);
    my $OPE = $production_rates{"Ox_$mechanism"} / $production_rates{"HNO3"} ;
    my $reshaped_OPE = $OPE->copy->reshape($n_per_day, $n_days);
    my $daytime_OPE = $reshaped_OPE->sumover;
    return $daytime_OPE;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 
