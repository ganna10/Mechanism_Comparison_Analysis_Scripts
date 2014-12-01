#! /usr/bin/env perl
# Plot Ox production and consumption budgets 
# Version 0: Jane Coates 26/11/2014

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
my (%plot_data, %families, %weights);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 O O1D NO2 NO3 HO2NO2 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")]);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$mechanism}}) {
        foreach my $species (sort keys %$ref) {
            my $chemical = get_chemical_name($species);
            $R->set('species', $chemical);
            $R->set('rate', [map { $_ } $ref->{$species}->dog]);
            $R->run(q` pre[species] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = gather(pre, Species, Rate, -Time) `,
            q` pre$Mechanism = rep(mechanism, length(pre$Time)) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c(  "Production Others" = "#696537", 
                            "Methane" = "#6c254f",
                            "CO" = "#f9c500", 
                            "Consumption Others" = "#0352cb", 
                            "O3 = UNITY" = "#0e5c28", 
                            "N2O5 = NA + NA" = "#e7e85e", 
                            "Butane" = "#ef6638", 
                            "Pentane" = "#b569b3", 
                            "2-Methylbutane" = "#4c9383", 
                            "NO2 + OH = HNO3" = "#86b650", 
                            "O1D = OH + OH" = "#cc6329", 
                            "Propane" = "#2b9eb3", 
                            "2-Methylpropene " = "#f7c56c",
                            "Isoprene " = "#0c3f78",
                            "Toluene" = "#8c1531",
                            "m-Xylene " = "#6db875",
                            "o-Xylene " = "#f3aa7f",
                            "p-Xylene " = "#be2448" ) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Species)) `,
        q` plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") `,
        q` plot = plot + geom_bar(data = subset(data, Rate > 0), stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "Ox_production_consumpion_budgets.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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
    } elsif ($VOC eq "OH + NO2 = HNO3") {
        $chemical_species = "NO2 + OH = HNO3";
    } else {
        print "No chemical species found for $VOC\n";
        $chemical_species = $VOC;
    }
    return $chemical_species;
}

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 )];
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
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production{$parent} += $rate(1:$NTIME-2);
            } else {
                my $reaction_string = $kpp->reaction_string($reaction);
                $production{$reaction_string} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $consumption{$parent} += $rate(1:$NTIME-2);
            } else {
                my $reaction_string = $kpp->reaction_string($reaction);
                $consumption{$reaction_string} += $rate(1:$NTIME-2);
            }
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others = 4e8;
    foreach my $item (keys %production) {
        if ($production{$item}->sum < $others) {
            $production{"Production Others"} += $production{$item};
            delete $production{$item};
        }
    }

    foreach my $item (keys %consumption) {
        if ($consumption{$item}->sum > -$others) {
            $consumption{"Consumption Others"} += $consumption{$item};
            delete $consumption{$item};
        }
    }

    foreach my $item (keys %production) {
        my $reshape = $production{$item}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production{$item} = $integrate;
    }

    foreach my $item (keys %consumption) {
        my $reshape = $consumption{$item}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $consumption{$item} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };

    my @sorted_prod = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    my @sorted_cons = reverse sort { &$sort_function($consumption{$b}) <=> &$sort_function($consumption{$a}) } keys %consumption;

    my @final_sorted_data;
    foreach (@sorted_cons) { #sum up rates of reactions, starting with reaction with lowest sum, consumption others added separately 
        next if ($_ eq 'Consumption Others');
        push @final_sorted_data, { $_ => $consumption{$_} };
    }

    push @final_sorted_data, { 'Consumption Others' => $consumption{'Consumption Others'} } if (defined $consumption{'Consumption Others'}); 
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @final_sorted_data, { $_ => $production{$_} };
    }

    push @final_sorted_data, { 'Production Others' => $production{'Production Others'} } if (defined $production{'Production Others'}); 
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
