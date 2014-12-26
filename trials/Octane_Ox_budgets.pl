#! /usr/bin/env perl
# Compare Ox production and consumption budgets from NC8H18 degradation 
# Version 0: Jane Coates 25/12/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MOZART-4_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
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
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes( x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") `,
        q` plot = plot + geom_bar(data = subset(data, Rate > 0), stat = "identity") `,
        q` plot = plot + facet_wrap(~ Mechanism) `,
);

$R->run(q` CairoPDF(file = "NC8H18_Ox_budgets.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
    my @loop = ("Ox_$mechanism", "HO2x");
    my (%production, %consumption);

    foreach my $species (@loop) {
        my ($producers, $producer_yields, $consumers, $consumer_yields);
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
            print "No family found for $species\n";
        }
        print "No producers found for $species\n" if (@$producers == 0);
        print "No consumers found for $species\n" if (@$consumers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            next unless ($reaction =~ /NC8H18|BIGALK|HC8/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            $rate *= 0.01 if ($mechanism eq "MOZART-4"); 
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $products) = split / = /, $reaction_string;
            $production{$reactants} += $rate(1:$NTIME-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            next unless ($reaction =~ /NC8H18|BIGALK|HC8/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            $rate *= 0.01 if ($mechanism eq "MOZART-4"); 
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $products) = split / = /, $reaction_string;
            $consumption{$reactants} += $rate(1:$NTIME-2);
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others = 9e5;
    foreach my $reaction (keys %production) {
        my $reshape = $production{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$reaction} = $integrate;
    }
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others) {
            $production{"Production Others"} += $production{$reaction};
            delete $production{$reaction};
        }
    }

    foreach my $reaction (keys %consumption) {
        my $reshape = $consumption{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $consumption{$reaction} = $integrate;
    }
    foreach my $reaction (keys %consumption) {
        if ($consumption{$reaction}->sum > -$others) {
            $consumption{"Consumption Others"} += $consumption{$reaction};
            delete $consumption{$reaction};
        }
    }

    my $sort_function = sub { $_[0]->sum } ;
    my @sorted_prod = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    my @sorted_cons = reverse sort { &$sort_function($consumption{$b}) <=> &$sort_function($consumption{$a}) } keys %consumption;
    my @sorted_data;
    foreach (@sorted_cons) {
        next if ($_ eq "Consumption Others") ;
        push @sorted_data, { $_ => $consumption{$_} };
    }
    push @sorted_data, { 'Consumption Others' => $consumption{'Consumption Others'} } if (defined $consumption{'Consumption Others'});
    foreach (@sorted_prod) {
        next if ($_ eq 'Production Others');
        push @sorted_data, { $_ => $production{$_} }
    }
    push @sorted_data, { 'Production Others' => $production{'Production Others'} } if (defined $production{'Production Others'});
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
