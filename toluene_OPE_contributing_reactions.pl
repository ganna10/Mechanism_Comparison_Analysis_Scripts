#! /usr/bin/env perl
# Calculate toluene OPE of each mechanism's model run base on Kleinman:2002 paper (P(Ox)/P(HNO3)) and express as contributing reactions
# Version 0: Jane Coates 15/10/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use Statistics::R;
use PDL;
use PDL::NiceSlice;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( CRI_tagging );
#my @mechanisms = qw( CRI );
my $index = 0;

my (%families, %weights, %plot_data);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my $mech_species = get_model_name('Toluene', $run);
    $plot_data{$mechanisms[$index]} = get_data($kpp, $mecca, $mechanisms[$index], $mech_species);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(grid) `,
        q` library(scales) `,
        q` library(gridExtra) `,
);

my @days = ("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7");
$R->set('Time', [@days]);
$R->run(q` plots = list() `);

$R->run(q` scientific_10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `);
$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "C6H5O2 + NO" = "#603912", "NO + RN10O2" = "#603912", "NO + TOLO2" = "#603912", "NO + TOLP" = "#603912",
                            "CH3O2 + NO" = "#6db875", "MO2 + NO" = "#6db875", "MEO2 + NO" = "#6db875",
                            "CH3CO3 + NO" = "#ae4901", "ACO3 + NO" = "#ae4901", "C2O3 + NO" = "#ae4901",
                            "CRES + OH" = "#898989", "CSL + OH" = "#898989",
                            "OH + TOL" = "black",
                            "MCATEC1O2 + NO" = "#8fd5d3",
                            "MECOACETO2 + NO" = "#bb8a01",
                            "NO + TLBIPERO2" = "#9bb08f",
                            "NO + TLFUO2" = "#5b671f",
                            "ACCOMECO3 + NO" = "#0e5c28",
                            "C2H5CO3 + NO" = "#f9c500",
                            "HOCH2CH2O2 + NO" = "#e7e85e",
                            "C2H5O2 + NO" = "#4b9483", "ETHP + NO" = "#4b9483",
                            "NO + RA16O2" = "#f8c56c",
                            "NO + TCO3" = "#6d6537",
                            "BAL1 + NO" = "#ef6638",
                            "BALP + NO" = "#cc6638",
                            "MCTP + NO" = "#4c9383",
                            "NO + PER1" = "#623812",
                            "NO + XO2" = "#0352cb", 
                            "KETP + NO" = "#c9a415",
                            "OH + ONIT" = "#6c254f", 
                            "HC3P + NO" = "#8ed6d2",
                            "NO + RCO3" = "#f3aa7f" ) `,
);

$R->run(q` plotting = function (data) { plot = ggplot(data, aes(x = Time, y = OPE, fill = reaction)) ;
                                        plot = plot + geom_bar(stat = "identity") ; 
                                        plot = plot + scale_y_continuous(limits = c(0, 1.5e-9), breaks = seq(0, 1.5e-9, 5e-10), labels = scientific_10) ;
                                        plot = plot + theme_bw() ;
                                        plot = plot + theme(axis.title.x = element_blank()) ;
                                        plot = plot + theme(axis.title.y = element_blank()) ;
                                        plot = plot + theme(strip.text.x = element_text(size = 160, face = "bold")) ;
                                        plot = plot + theme(strip.text.y = element_text(size = 160, face = "bold", angle = 0)) ;
                                        plot = plot + theme(strip.background = element_blank()) ;
                                        plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                        plot = plot + theme(axis.text.x = element_text(size = 140)) ;
                                        plot = plot + theme(panel.grid.major = element_blank()) ;
                                        plot = plot + theme(panel.grid.minor = element_blank()) ;
                                        plot = plot + theme(legend.key = element_blank()) ;
                                        plot = plot + theme(legend.key.size = unit(8, "cm")) ;
                                        plot = plot + theme(legend.title = element_blank()) ;
                                        plot = plot + theme(legend.text = element_text(size = 140)) ;
                                        plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                        plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                        plot = plot + scale_fill_manual(values = my.colours) ;
                                        return (plot) } `,
);

foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}} ) {
        foreach my $reaction (keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('OPE', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = OPE `);
        }
    }
    $R->run(q` data = data[order(Time),] `,
            q` data = data[1:7,] `,
            q` data = melt(data, id.vars = c("Time"), variable.name = "reaction", value.name = "OPE") `,
            q` reaction.levels = rev(levels(factor(data$reaction))) `,
            q` data$reaction = ordered(data$reaction, levels = reaction.levels) `,
            q` plot = plotting(data) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(length(plots)) `);
#print "$p\n"; 

$R->run(q` CairoPDF(file = "toluene_OPE_comparison.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[2]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[6]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[7]], 
                                                    plots[[8]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[9]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       sub = textGrob("\n", gp = gpar(fontsize = 30)) ,
                                       left = textGrob(expression(bold(paste("\nNormalised Ox Production Efficiency (molecules (Ox) / (molecules (VOC) molecules (HNO3))) ", s^-1))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC) = @_;
    $families{"Ox_${mechanism}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_${mechanism}", "HO2x_${mechanism}", "HNO3");

    my ($producers, $producer_yields, %production_rates);
    foreach my $species (@loop) {
        if (exists $families{$species}) {
            $kpp->family({
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
        } else {
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
        }
        print "No producers found for $species\n" if (@$producers == 0);
        
        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            if ($species eq "HNO3") {
                my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number);
                next if ($rate->sum == 0);
                $production_rates{$species} += $rate(1:$NTIME-2);
            } elsif ($species =~ /Ox/) {
                my ($number, $parent) = split /_/, $reaction;
                next unless (defined $parent and $parent eq $VOC);
                my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number);
                next if ($rate->sum == 0);
                my $reaction_string = $kpp->reaction_string($reaction);
                $reaction_string =~ s/_(.*?)\b//g;
                my ($reactants, $products) = split / = /, $reaction_string;
                if ($reactants =~ /XO2/ and $mechanism =~ /RADM2|RACM|CB/) {
                    my $operator = "XO2_" . $VOC;
                    my $op_producers = $kpp->producing($operator);
                    my $op_producer_yields = $kpp->effect_on($operator, $op_producers); 
                    die "No producers found for $operator\n" if (@$op_producers == 0);

                    for (0..$#$op_producers) { #get rates for all producing reactions
                        my $reaction = $op_producers->[$_];
                        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
                        my $reaction_number = $kpp->reaction_number($reaction);
                        my $rate = $op_producer_yields->[$_] * $mecca->rate($reaction_number); 
                        next if ($rate->sum == 0); # do not include reactions that do not occur 
                        my $reaction_string = $kpp->reaction_string($reaction);
                        $reaction_string =~ s/_(.*?)\b//;
                        my ($reactants, $products) = split / = /, $reaction_string;
                        $production_rates{$reactants} += $rate(1:$NTIME-2);
                    } 
                } else {
                    $production_rates{$reactants} += $rate(1:$NTIME-2);
                }
            }
        } 
    }

    my $prod_others = 5e6;
    foreach my $reaction (sort keys %production_rates) {
        next if ($reaction eq "HNO3");
        if ($production_rates{$reaction}->sum < $prod_others) {
            $production_rates{"Production Others"} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
    }

    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($mechanism =~ /CBM-IV/ or $mechanism =~ /CB05/) {
        my $name = "PAR_NC5H12";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
    } else {
        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (Product) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    foreach my $reaction (sort keys %production_rates) {
        next if ($reaction eq "HNO3");
        $production_rates{$reaction} /= $parent_emissions;
    }

    my $n_per_day = 43200 / $dt;
    my $n_days = int ($NTIME / $n_per_day);
    my %daily_OPEs;
    foreach my $reaction (sort keys %production_rates) {
        next if ($reaction eq "HNO3");
        print "$reaction\n";
        my $OPE = $production_rates{$reaction} / $production_rates{"HNO3"} ;
        my $reshaped_OPE = $OPE->copy->reshape($n_per_day, $n_days);
        my $daytime_OPE = $reshaped_OPE->sumover;
        $daily_OPEs{$reaction} = $daytime_OPE;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($daily_OPEs{$b}) <=> &$sort_function($daily_OPEs{$a}) } keys %daily_OPEs;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $daily_OPEs{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $daily_OPEs{'Production Others'} } if (defined $daily_OPEs{'Production Others'}); 

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

sub get_model_name {
    my ($VOC, $mechanism) = @_;
    my $species;
    if ($VOC eq "Pentane") {
        if ($mechanism  =~ /MCM|CRI|CB/) {
            $species = "NC5H12";
        } elsif ($mechanism =~ /MOZART/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RADM|RACM/) {
            $species = "HC5";
        } else {
            print "No mechanism for $VOC\n";
        }
    } elsif ($VOC eq "Toluene") {
        if ($mechanism =~ /MCM|CRI|MOZART|CB/) {
            $species = "TOLUENE";
        } elsif ($mechanism =~ /RADM|RACM/) {
            $species = "TOL";
        } else {
            print "No mechanism for $VOC\n";
        }
    } else {
        print "No species found for $VOC\n";
    }
}
