#! /usr/bin/env perl
# Plot VOCs contibuting to Total OxPE in all mechanisms
# Version 0: Jane Coates 27/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2",  "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my $index = 0;
my (%data, %families, %weights);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    ($data{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            next if ($reaction =~ /O3 \+ OH|OH \+ O3/);
            my $chemical = get_chemical_name($reaction);
            $R->set('reaction', $chemical);
            $R->set('OxPE', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` pre[reaction] = OxPE `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Reaction, OxPE, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(pre) `);
#print $p, "\n";
$R->run(q` my.colours = c(  "Others" = "#696537", 
                            "Methane" = "#6c254f",
                            "CO" = "#f9c500", 
                            "Ethane" = "#0352cb", 
                            "Propane" = "#0e5c28", 
                            "2-Methylpropane" = "#e7e85e", 
                            "Butane" = "#ef6638", 
                            "Pentane" = "#b569b3", 
                            "2-Methylbutane" = "#4c9383", 
                            "Hexane" = "#86b650", 
                            "Ethene" = "#cc6329", 
                            "Propene" = "#2b9eb3", 
                            "2-Methylpropene" = "#f7c56c",
                            "Isoprene" = "#0c3f78",
                            "Toluene" = "#8c1531",
                            "m-Xylene" = "#6db875",
                            "o-Xylene" = "#f3aa7f",
                            "p-Xylene" = "#be2448" ) `); 
$R->run(q` data$Reaction = factor(data$Reaction, levels = c("Methane", "CO", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Ethene", "Propene", "2-Methylpropene", "Isoprene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Others")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 10)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Reaction))) `,
);

$R->run(q` CairoPDF(file = "Ox_PE_VOCs.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_$mechanism", "HO2x");
    my (%production, $consumption, %OxPE);

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
        }
        
        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production{$species}{$parent} += $rate(1:$NTIME-2);
            } else { 
                my $reaction_string = $kpp->reaction_string($reaction);
                $production{$species}{$reaction_string} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) {
            last if ($species eq "HO2x");
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            $consumption += $rate(1:$NTIME-2);
        }
    } 
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production{'HO2x'}{$_} for (keys %{ $production{'HO2x'} });

    foreach my $reaction( keys %{ $production{'HO2x'} }) {
        if ($mechanism =~ /CRI/) {
            $production{"Ox_$mechanism"}{$reaction} += $production{"Ox_$mechanism"}{'HO2 + NO = OH + NO2'} * $production{'HO2x'}{$reaction} / $ho2x_total_production;
        } else {
            $production{"Ox_$mechanism"}{$reaction} += $production{"Ox_$mechanism"}{'HO2 + NO = NO2 + OH'} * $production{'HO2x'}{$reaction} / $ho2x_total_production;
        }
    }
    delete $production{"Ox_$mechanism"}{'HO2 + NO = NO2 + OH'} if (defined $production{"Ox_$mechanism"}{'HO2 + NO = NO2 + OH'});
    delete $production{"Ox_$mechanism"}{'HO2 + NO = OH + NO2'} if (defined $production{"Ox_$mechanism"}{'HO2 + NO = OH + NO2'});

    my $others = 1e8;
    foreach my $reaction (keys %{$production{"Ox_$mechanism"}}) {
        if ($production{"Ox_$mechanism"}{$reaction}->sum < $others) {
            $production{"Ox_$mechanism"}{"Others"} += $production{"Ox_$mechanism"}{$reaction};
            delete $production{"Ox_$mechanism"}{$reaction};
        }
    }

    foreach my $reaction (keys %{$production{"Ox_$mechanism"}}) {
        my $OxPE = $production{"Ox_$mechanism"}{$reaction} / -$consumption;
        my $reshape = $OxPE->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $OxPE{$reaction} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum }; 
    my @sorted = sort { &$sort_function($OxPE{$b}) <=> &$sort_function($OxPE{$a}) } keys %OxPE;
    
    my @final_sorted_data;
    foreach (@sorted) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Others');
        push @final_sorted_data, { $_ => $OxPE{$_} };
    }

    push @final_sorted_data, { 'Others' => $OxPE{'Others'} } if (defined $OxPE{'Others'}); 
    return \@final_sorted_data;
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

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'CO + OH = HO2' or $VOC eq 'OH + CO = HO2') {
        $chemical_species = 'CO';
    } elsif ($VOC eq 'CH4') {
        $chemical_species = 'Methane';
    } elsif ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane';
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene';
    } elsif ($VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene";
    } elsif ($VOC eq 'BIGENE') {
        $chemical_species = "Butene";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene";
    } elsif ($VOC eq 'Production Others') {
        $chemical_species = 'Production Others';
    } elsif ($VOC eq 'Consumption Others') {
        $chemical_species = 'Consumption Others';
    } elsif ($VOC eq 'Others') {
        $chemical_species = 'Others';
    } elsif ($VOC eq "OH + NO2 = HNO3") {
        $chemical_species = "NO2 + OH = HNO3";
    } else {
        print "No chemical species found for $VOC\n";
        $chemical_species = $VOC;
    }
    return $chemical_species;
}
