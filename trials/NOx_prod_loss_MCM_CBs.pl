#! /usr/bin/env perl
# compare NOx family production and loss budgets in CBM-IV and CB05 to MCM v3.2
# Version 0: Jane Coates 29/1/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 CBM-IV CB05 );
#my @mechanisms = qw( CB05 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    $families{$mechanism} = [ qw( NO NO2 NO3 N2O5 HO2NO2 CH3O2NO2_CH4 CH3O2NO2_C2H6 CH3O2NO2_C3H8 CH3O2NO2_NC4H10 CH3O2NO2_IC4H10 CH3O2NO2_NC5H12 CH3O2NO2_IC5H12 CH3O2NO2_NC6H14 CH3O2NO2_NC7H16 CH3O2NO2_NC8H18 CH3O2NO2_C2H4 CH3O2NO2_C3H6 CH3O2NO2_MEPROPENE CH3O2NO2_BUT1ENE CH3O2NO2_C5H8 CH3O2NO2_BENZENE CH3O2NO2_TOLUENE CH3O2NO2_MXYL CH3O2NO2_PXYL CH3O2NO2_OXYL CH3O2NO2_EBENZ ) ];
    $weights{$mechanism} = { N2O5 => 2 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog ]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#$R->run(q` HO2NO2 = filter(data, Reaction == "HO2NO2" | Reaction == "HO2 + NO2") `);
#my $p = $R->run(q` print(HO2NO2) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction, group = Reaction)) `,
        q` plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") `,
        q` plot = plot + geom_bar(data = subset(data, Rate > 0), stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
);

$R->run(q` CairoPDF(file = "NOx_loss_production_budgets_by_reactions.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    my ($producers, $producer_yields, $consumers, $consumer_yields, %production_rates, %consumption_rates);
    if (exists $families{$mechanism}) {
        $kpp->family({
                name    => $mechanism,
                members => $families{$mechanism},
                weights => $weights{$mechanism},
        });
        $producers = $kpp->producing($mechanism);
        $producer_yields = $kpp->effect_on($mechanism, $producers);
        $consumers = $kpp->consuming($mechanism);
        $consumer_yields = $kpp->effect_on($mechanism, $consumers);
    } else {
        print "No NOx family found for $mechanism\n";
    }
    print "No producers found for NOx family in $mechanism\n" if (@$producers == 0);
    print "No consumers found for NOx family in $mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $production_rates{$reactants} += $rate(1:$NTIME-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $consumption_rates{$reactants} += $rate(1:$NTIME-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

    my $others_max = 5e8;
    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        if ($integrate->sum < $others_max) {
            $production_rates{"Production Others"} += $integrate(0:13:2);
            delete $production_rates{$reaction};
        } else {
            $production_rates{$reaction} = $integrate(0:13:2);
        }
    }

    foreach my $reaction (keys %consumption_rates) {
        my $reshape = $consumption_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        if ($integrate->sum > -$others_max) {
            $consumption_rates{"Consumption Others"} += $integrate(0:13:2);
            delete $consumption_rates{$reaction};
        } else {
            $consumption_rates{$reaction} = $integrate(0:13:2);
        }
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production_rates{$b} <=> $production_rates{$a}) } keys %production_rates;
    my @sorted_cons = reverse sort { &$sort_function($consumption_rates{$b} <=> $consumption_rates{$a}) } keys %consumption_rates;
    my @sorted_data;

    foreach (@sorted_cons) {
        next if ($_ =~ /Others/);
        push @sorted_data, { $_ => $consumption_rates{$_} };
    }
    push @sorted_data, { "Consumption Others" => $consumption_rates{"Consumption Others"} } if (defined $consumption_rates{"Consumption Others"});
    foreach (@sorted_prod) {
        next if ($_ =~ /Others/);
        push @sorted_data, { $_ => $production_rates{$_} };
    }
    push @sorted_data, { "Production Others" => $production_rates{"Production Others"} } if (defined $production_rates{"Production Others"});
    return \@sorted_data;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
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
