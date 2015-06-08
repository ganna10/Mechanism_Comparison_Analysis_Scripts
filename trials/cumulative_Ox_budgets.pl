#!/usr/bin/env perl 
# cumulative Ox production budget plots from all tagged mechanisms at end of model run
# Version 0: Jane Coates 8/6/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base_dir = "/local/home/coates/MECCA/";
my $mecca = MECCA->new("$base_dir/CB05_tagged/boxmodel"); 
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my (%data, %families, %weights);
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 ); 
#my @mechanisms = qw( RACM2 RADM2 );

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base_dir/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base_dir/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base_dir/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]; 
    $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
    $data{$mechanism} = get_data($mecca, $kpp);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
        q` library(dplyr) `,
);

$R->run(q` data = data.frame(Dummy = c(1)) `);
foreach my $mechanism (keys %data) {
    $R->set('mechanism', $mechanism);
    $R->set('ox', $data{$mechanism});
    $R->run(q` data[mechanism] = ox `);
}

$R->run(q` data = gather(data, Mechanism, Ox) `,
        q` data = data[-1,] `,
);
$R->run(q` write.table(data, file = "cumulative_Ox.csv", sep = ",", quote = FALSE, row.names = FALSE) `);
#my $p = $R->run(q` print(calc) `);
#print $p, "\n";

$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Ox)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ylab("Reaction Rate (molecules cm-3)") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 0.4)) `, 
);

$R->run(q` CairoPDF( file = "Cumulative_Ox_production_budgets.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp) = @_;
    my ($production, $producers, $producer_yields);
    my $species = "Ox";
    if (exists $families{$species}) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers);
    } else {
        print "No family for $species\n";
    }
    print "No producers found for $species\n" if (@$producers == 0);
    
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $production += $rate->sum;
    }

    return $production;
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
