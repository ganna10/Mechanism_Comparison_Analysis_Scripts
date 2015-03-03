#! /usr/bin/env perl
# XO2 production budget in IC5H12 degradation, as this has highest emissions of alkanes with C >= 4
# Version 0: Jane Coates 2/3/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @mechanisms = qw( CBM-IV CB05 );
my %data;

my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file); 
    my $producers = $kpp->producing("XO2_IC5H12");
    my $producer_yields = $kpp->effect_on("XO2_IC5H12", $producers);
    print "No producers found in $mechanism\n" if (@$producers == 0);
    
    my %production_rates;
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        print "$reaction_string : $producer_yields->[$_]\n";
        $production_rates{$reactants} = $rate(1:$NTIME-2);
    }

    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production_rates{$reaction} = $integrate(0:13:2);
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    my @final_sorted;
    push @final_sorted, { $_ => $production_rates{$_} } foreach (@sorted_prod);
    $data{$mechanism} = \@final_sorted;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `,
            q` pre$Mechanism = rep(mechanism, length(Time)) `,
    );
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog ]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->run(q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `, 
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism )`,
);

$R->run(q` CairoPDF(file = "XO2_budget.pdf") `,
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
