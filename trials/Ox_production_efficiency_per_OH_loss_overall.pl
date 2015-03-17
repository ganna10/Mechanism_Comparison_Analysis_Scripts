#! /usr/bin/env perl
# Calculate OxPE for each of the degradation product sizes of pentane and toluene degradation by normalising by total OH loss rate
# Version 0: Jane Coates 17/3/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
#my @mechanisms = qw( RADM2 );
my (%n_carbon, %families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn);
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O1D O NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };

    my $oh_loss_rate = 0; #calculate total OH loss rate for all reactions
    my $oh_consumers = $kpp->consuming("OH");
    my $oh_consumer_yields = $kpp->effect_on("OH", $oh_consumers);
    print "No consumers for OH\n"  if (@$oh_consumers == 0);
    for (0..$#$oh_consumers) {
        my $reaction = $oh_consumers->[0];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $oh_consumer_yields->[$_];
        next if ($rate->sum == 0);
        $oh_loss_rate = $rate(1:$NTIME-2);
    }
    $oh_loss_rate = $oh_loss_rate->reshape($N_PER_DAY, $N_DAYS);
    $oh_loss_rate = $oh_loss_rate->sumover;
    
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism, $oh_loss_rate);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    $R->set('mechanism', $mechanism);
    $R->set('oxpe', [ map { $_ } $data{$mechanism}->dog ]);
    $R->run(q` pre[mechanism] = oxpe `,
            q` pre = gather(pre, Mechanism, OxPE, -Time) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#ef6638", "CRIv2" = "#b569b3", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + geom_point() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "Ox_PE_overall_per_OH_loss.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $oh_loss_rate) = @_;
    
    my ($production_rates, $consumption_rates);
    my ($producers, $consumers, $producer_yields, $consumer_yields);
    if (exists $families{"Ox_$mechanism"}) {
        $kpp->family({
                name    => "Ox_$mechanism",
                members => $families{"Ox_$mechanism"},
                weights => $weights{"Ox_$mechanism"},
        });
        $producers = $kpp->producing("Ox_$mechanism");
        $producer_yields = $kpp->effect_on("Ox_$mechanism", $producers);
        $consumers = $kpp->consuming("Ox_$mechanism");
        $consumer_yields = $kpp->effect_on("Ox_$mechanism", $consumers);
    } else {
        print "No family for Ox_$mechanism\n";
    }
    print "No producers for Ox_$mechanism\n" if (@$producers == 0);
    print "No consumers for Ox_$mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_]; 
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0); 
        $production_rates += $rate(1:$NTIME-2);
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $consumption_rates += $rate(1:$NTIME-2);
    }

    $production_rates = $production_rates->reshape($N_PER_DAY, $N_DAYS);
    $production_rates = $production_rates->sumover;
    $consumption_rates = $consumption_rates->reshape($N_PER_DAY, $N_DAYS);
    $consumption_rates = $consumption_rates->sumover;
    $production_rates += $consumption_rates; #net production
    $production_rates /= -$oh_loss_rate; #normalise by total OH loss rate 
    return $production_rates(0:13:2);
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
