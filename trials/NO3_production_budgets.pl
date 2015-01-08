#! /usr/bin/env perl
# NO3 production budgets in each mechanism
# Version 0: Jane Coates 8/1/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $N_TIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $N_TIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 RADM2 );
my (%families, %weights, %data);
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    $families{"NO3_$mechanism"} = [ qw( O3 NO3 NO2 N2O5 ) ];
    $weights{"NO3_$mechanism"} = { N2O5 => 2 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
);

$R->set('Time', [("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
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
        q` plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") `,
        q` plot = plot + geom_bar(data = subset(data, Rate > 0), stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "NO3_production_budgets_reactions.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    my ($producers, $producer_yields, $consumers, $consumer_yields, %production, %consumption);
    my $species = "NO3_$mechanism";
    if (exists $families{$species}) {
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
        print "No family for $species\n";
    }
    print "No producers for $species\n" if (@$producers == 0);
    print "No consumers for $species\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $string = $kpp->reaction_string($reaction);
        $string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $string;
        $production{$reactants} += $rate(1:$N_TIME-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $string = $kpp->reaction_string($reaction);
        $string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $string;
        $consumption{$reactants} += $rate(1:$N_TIME-2);
    }
    remove_common_processes(\%production, \%consumption);

    my $others = 6e8;
    foreach my $item (keys %consumption) {
        if ($consumption{$item}->sum > -$others) {
            $consumption{"C Others"} += $consumption{$item};
            delete $consumption{$item};
        }
    }
    
    foreach my $item (keys %consumption) {
        my $reshape = $consumption{$item}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $consumption{$item} = $integrate;
    }

    foreach my $item (keys %production) {
        if ($production{$item}->sum < $others) {
            $production{"Others"} += $production{$item};
            delete $production{$item};
        }
    }
    
    foreach my $item (keys %production) {
        my $reshape = $production{$item}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production{$item} = $integrate;
    }
    
    my $sort_function = sub { $_[0]->sum } ;
    my @sorted_prod = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    my @sorted_cons = reverse sort { &$sort_function($consumption{$b}) <=> &$sort_function($consumption{$a}) } keys %consumption;
    my @sorted_data;
    foreach (@sorted_cons) {
        next if ($_ eq "C Others");
        push @sorted_data, { $_ => $consumption{$_} };
    }
    push @sorted_data, {"C Others" => $consumption{"C Others"} } if (defined $consumption{"C Others"});
    foreach (@sorted_prod) {
        next if ($_ eq 'Others');
        push @sorted_data, { $_ => $production{$_} }
    }
    push @sorted_data, { 'Others' => $production{'Others'} } if (defined $production{'Others'});
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
