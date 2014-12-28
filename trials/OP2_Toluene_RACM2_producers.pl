#! /usr/bin/env perl
# OP2 production contributions during Toluene degradation 
# Version 0: Jane Coates 27/12/2014

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/RACM2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0); 
my $n_per_day = 43200 / $dt;
my $n_days = int ($NTIME / $n_per_day);
my (%families, %weights, %plot_data);

my $eqnfile = "$base/RACM2_tagged/gas.eqn";
my $kpp = KPP->new($eqnfile); 
my ($consumers, $producers, $consumer_yields, $producer_yields, %production_reaction_rates, %consumption_reaction_rates);
my $species = "OP2_TOL";
$consumers = $kpp->consuming($species);
$producers = $kpp->producing($species);
$consumer_yields = $kpp->effect_on($species, $consumers); 
$producer_yields = $kpp->effect_on($species, $producers);  
die "No producers found for $species\n" if (@$producers == 0);
die "No consumers found for $species\n" if (@$consumers == 0);

my $prod_others_max = 8e5;
for (0..$#$producers) { #get rates for all producing reactions
    my $reaction = $producers->[$_];
    my $reaction_number = $kpp->reaction_number($reaction);
    my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
    next if ($rate->sum == 0); # do not include reactions that do not occur 
    my $reaction_string = $kpp->reaction_string($reaction);
    $reaction_string =~ s/_(.*?)\b//g;
    my ($reactants, $products) = split / = /, $reaction_string;
    $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
}

for (0..$#$consumers) { #get rates for all consuming reactions
    my $reaction = $consumers->[$_];
    my $reaction_number = $kpp->reaction_number($reaction);
    my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
    next if ($rate->sum == 0); # do not include reactions that do not occur
    my $reaction_string = $kpp->reaction_string($reaction);
    $reaction_string =~ s/_(.*?)\b//g;
    my ($reactants, $products) = split / = /, $reaction_string;
    $consumption_reaction_rates{$reactants} += $rate(1:$NTIME-2);
} 

remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
foreach my $reaction (keys %production_reaction_rates) {
    if ($production_reaction_rates{$reaction}->sum < $prod_others_max) {
        $production_reaction_rates{"Production Others"} += $production_reaction_rates{$reaction};
        delete $production_reaction_rates{$reaction};
    }
}

foreach my $reaction (keys %production_reaction_rates) {
    my $reshape = $production_reaction_rates{$reaction}->copy->reshape($n_per_day, $n_days);
    my $integrate = $reshape->sumover;
    $integrate = $integrate(0:13:2);
    $production_reaction_rates{$reaction} = $integrate;
}

my $sort_function = sub { $_[0]->sum };
my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;

my @final_sorted_data;
foreach (@prod_sorted_data) { 
    next if ($_ eq 'Production Others') ;
    push @final_sorted_data, { $_ => $production_reaction_rates{$_} };
} 
push @final_sorted_data, { 'Production Others' => $production_reaction_rates{'Production Others'} } if (defined $production_reaction_rates{'Production Others'}); 
$plot_data{"RACM2"} = \@final_sorted_data;

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
);

$R->run(q` CairoPDF(file = "Toluene_OP2_production_contributors.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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
