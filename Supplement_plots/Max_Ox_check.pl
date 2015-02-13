#! /usr/bin/env perl
# plot total Ox and O3 production on day 1 over NO emission percentages
# Version 0: Jane Coates 7/2/2015 

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/MECCA/Max_Ox_Check";
my $mecca = MECCA->new("$base/CRIv2_tagged_1.0/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 43200 / $dt;
my $n_days = int $ntime / $n_per_day;

my @runs = qw( tagged_0.5 tagged_0.8 tagged_0.9 tagged_1.0 tagged_1.1 tagged_1.2 tagged_1.5 );
my @mechanism = qw(MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05) ;
my (%families, %weights, %data);

foreach my $mechanism (@mechanism) {
    foreach my $run (@runs) {
        my $boxmodel = "$base/${mechanism}_$run/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/${mechanism}_$run/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $ro2_file = "$base/${mechanism}_$run/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox"} = [ qw(O3 O1D O NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
        (my $label = $run) =~ s/tagged_//;
        $data{$mechanism}{$label} = get_data($mecca, $kpp);
    }
}

foreach my $mechanism (sort keys %data) {
    my @Ox_prod = ();
    my @O3_prod = ();
    foreach my $run (sort keys %{$data{$mechanism}}) {
        push @Ox_prod, $data{$mechanism}{$run}{'Ox'};
        push @O3_prod, $data{$mechanism}{$run}{'O3'};
    }
    $data{$mechanism}{"Ox"} = \@Ox_prod;
    $data{$mechanism}{"O3"} = \@O3_prod;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(scales) `,
);

$R->set('NO.Percent', [qw(0.5 0.8 0.9 1.0 1.1 1.2 1.5)]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(NO.Percent) `);
    $R->set('mechanism', $mechanism);
    $R->set('ox.prod', [@{$data{$mechanism}{'Ox'}}]);
    $R->set('o3.prod', [@{$data{$mechanism}{'O3'}}]);
    $R->run(q` pre$Ox = ox.prod `,
            q` pre$O3 = o3.prod `,
            q` pre$Mechanism = rep(mechanism, length(NO.Percent)) `,
            q` pre = gather(pre, Item, Production, -NO.Percent, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `);

$R->run(q` plot = ggplot(data, aes(x = NO.Percent, y = Production, colour = Item)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + xlab("Percent of NO Emissions") `,
        q` plot = plot + ylab("Net Production (molecules cm-3)") `,
        q` plot = plot + scale_x_continuous(labels = percent) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.7)) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + scale_colour_manual(values = c("#6c254f", "#4c9383"), limits = c("Ox", "O3")) `,
        q` plot = plot + theme(legend.position = c(1, 0.8)) `,
        q` plot = plot + theme(legend.justification = c(1, 0.8)) `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_NO_emissions.pdf", width = 8, height = 11) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp) = @_;
    my (%production, %consumption, $producers, $producer_yields, $consumers, $consumer_yields);
    
    foreach my $species (qw( O3 Ox )) {
        if (defined $families{$species}) {
            $kpp->family({
                    name    => $species,
                    members => $families{$species},
                    weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        } else { 
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        }
        print "No producers found for $species\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            $production{$species} += $rate;
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            $consumption{$species} += $rate;
        }
    }
    foreach my $species (keys %production) {
        $production{$species} += $consumption{$species}; 
        $production{$species} = $production{$species}->reshape($n_per_day, $n_days);
        $production{$species} = $production{$species}->sumover;
        $production{$species} = $production{$species}->at(0) * $dt;
    }
    return \%production; 
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open FILE, $file or die $!; 
    my @ro2;
    for (<FILE>) {
        push @ro2, split /\s+/, $_; 
    }
    close FILE;
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

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
}
