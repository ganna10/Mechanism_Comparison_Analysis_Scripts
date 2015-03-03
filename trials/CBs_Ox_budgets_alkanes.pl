#! /usr/bin/env perl
# Plot reactions contributing to Ox production in CB05, CBM-IV for all alkanes
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
my (%families, %weights, %data);

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
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O1D O NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $parent (sort keys %{$data{$mechanism}}) {
        $R->set('alkane', $parent);
        $R->run(q` pre = data.frame(Time) `,
                q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre$Alkane = rep(alkane, length(Time)) `,
        );
        foreach my $ref (@{$data{$mechanism}{$parent}}) {
            foreach my $reaction (sort keys %$ref) {
                $R->set('reaction', $reaction);
                $R->set('rates', [ map { $_ } $ref->{$reaction}->dog ]);
                $R->run(q` pre[reaction] = rates `);
            }
        }
        $R->run(q` pre = gather(pre, Reaction, Rates, -Time, -Alkane, -Mechanism) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Alkane = factor(data$Alkane, levels = c("NC4H10", "IC4H10", "NC5H12", "IC5H12", "NC6H14", "NC7H16", "NC8H18")) `);
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("CBM-IV", "CB05")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rates, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid( Mechanism ~ Alkane ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.text.y = element_text(angle = 0)) `,
);

$R->run(q` CairoPDF(file = "CBMs_Alkanes_Ox_production_budgets.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
    my @loop = ("Ox_$mechanism", "HO2x");
    my (%production_rates, %consumption_rates);

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
            print "No family for $species in $mechanism\n";
        } 
        print "No producers of $species in $mechanism\n" if (@$producers == 0);
        print "No consumers of $species in $mechanism\n" if (@$consumers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            next unless ($reaction =~ /NC|IC/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            print "$reaction_string : $producer_yields->[$_]\n";
            my ($reactants, $products) = split / = /, $reaction_string;
#            $reactants = "OH + PAR" if ($reactants eq "ROR");
#            if ($reactants =~ /XO2/) {
#                my $operator = "XO2_" . $parent;
#                my $op_producers = $kpp->producing($operator);
#                my $op_producer_yields = $kpp->effect_on($operator, $op_producers); 
#                die "No producers found for $operator\n" if (@$op_producers == 0);
#
#                for (0..$#$op_producers) { #get rates for all producing reactions
#                    my $reaction = $op_producers->[$_];
#                    my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
#                    my $reaction_number = $kpp->reaction_number($reaction);
#                    my $rate = $op_producer_yields->[$_] * $mecca->rate($reaction_number); 
#                    next if ($rate->sum == 0); # do not include reactions that do not occur 
#                    my $reaction_string = $kpp->reaction_string($reaction);
#                    $reaction_string =~ s/_(.*?)\b//;
#                    my ($reactants, $products) = split / = /, $reaction_string;
#                    $reactants = "OH + PAR" if ($reactants eq "ROR");
#                    $production_rates{$parent}{$reactants} += $rate(1:$NTIME-2);
#                } 
#            } else {
                $production_rates{$parent}{$reactants} += $rate(1:$NTIME-2);
                #}
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            next unless ($reaction =~ /NC|IC|C2H6|C3H8/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $products) = split / = /, $reaction_string;
            #$reactants = "OH + PAR" if ($reactants eq "ROR");
            $consumption_rates{$parent}{$reactants} += $rate(1:$NTIME-2);
        }
    }

    my $others = 1e7;
    foreach my $parent (sort keys %production_rates) {
        remove_common_processes($production_rates{$parent}, $consumption_rates{$parent});
        foreach my $reaction (sort keys %{$production_rates{$parent}}) {
            my $reshape = $production_rates{$parent}{$reaction}->reshape($N_PER_DAY, $N_DAYS);
            my $integrate = $reshape->sumover;
            if ($integrate->sum < $others) {
                $production_rates{$parent}{"Others"} += $integrate(0:13:2);
                delete $production_rates{$parent}{$reaction};
            } else {
                $production_rates{$parent}{$reaction} = $integrate(0:13:2);
            }
        }
    } 

    my $sort_function = sub { $_[0]->sum };
    foreach my $parent (sort keys %production_rates) {
        my @sorted_prod = sort { &$sort_function($production_rates{$parent}{$b}) <=> &$sort_function($production_rates{$parent}{$a}) } keys %{$production_rates{$parent}};
        my @final_sorted;
        foreach (@sorted_prod) {
            next if ($_ =~ /Others/);
            push @final_sorted, { $_ => $production_rates{$parent}{$_} };
        }
        push @final_sorted, { "Others" => $production_rates{$parent}{"Others"} } if (defined $production_rates{$parent}{"Others"});
        $production_rates{$parent} = \@final_sorted;
    } 
    return \%production_rates;
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
