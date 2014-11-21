#! /usr/bin/env perl
# OH production during CH4 degradation
# Version 0: Jane Coates 19/11/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#my @runs = qw( MOZART_tagging );
#my @mechanisms = qw( MOZART-4 );
my @runs = qw( MCM_3.2_tagged MOZART_tagging  );
my @mechanisms = ( "MCMv3.2", "MOZART-4" );
my $index = 0;
my (%families, %weights, %plot_data);
$families{"HO2x"} = [ qw(HO2 HO2NO2) ];

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);

    my ($producers, $producer_yields, $consumers, $consumer_yields, %production, %consumption);
    my @loop = ("OH", "HO2x");
    foreach my $species (@loop) {
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
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        }
        print "No consumers found for OH in $run\n" if (@$consumers == 0);
        print "No producers found for OH in $run\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq "CH4");
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            my ($reactants, $products) = split / = /, $reaction_string;
            $production{$reactants} += $rate(1:$NTIME-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq "CH4");
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            my ($reactants, $products) = split / = /, $reaction_string;
            $consumption{$reactants} += $rate(1:$NTIME-2);
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others_max = 2e5;
    foreach (keys %production) {
        if ($production{$_}->sum < $others_max) {
            $production{"Others"} += $production{$_};
            delete $production{$_};
        }
    }

    foreach my $reaction (keys %production) {
        my $reshape = $production{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production{$reaction} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Others') ;
        push @final_sorted_data, { $_ => $production{$_} };
    } 
    push @final_sorted_data, { 'Others' => $production{'Others'} } if (defined $production{'Others'}); 

    $plot_data{$mechanisms[$index]} = \@final_sorted_data;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [ ("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7") ]);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "Reaction", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill =  Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
);

$R->run(q` CairoPDF(file = "CH4_OH_production_budgets.pdf") `,
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
