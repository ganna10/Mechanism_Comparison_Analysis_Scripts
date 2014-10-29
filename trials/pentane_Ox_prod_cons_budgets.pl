#! /usr/bin/env perl
# pentane Ox production and consumption budget for each mechanism
# Version 0: Jane Coates 24/10/2014

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
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog;
my @time_blocks;
foreach my $time (@time_axis) {
    if ($time <= 12) {
        push @time_blocks, "Day 1";
    } elsif ($time > 12 and $time <= 24) {
        push @time_blocks, "Night 1";
    } elsif ($time > 24 and $time <= 36) {
        push @time_blocks, "Day 2";
    } elsif ($time > 36 and $time <= 48) {
        push @time_blocks, "Night 2";
    } elsif ($time > 48 and $time <= 60) {
        push @time_blocks, "Day 3",
    } elsif ($time > 60 and $time <= 72) {
        push @time_blocks, "Night 3";
    } elsif ($time > 72 and $time <= 84) {
        push @time_blocks, "Day 4";
    } elsif ($time > 84 and $time <= 96) {
        push @time_blocks, "Night 4";
    } elsif ($time > 96 and $time <= 108) {
        push @time_blocks, "Day 5";
    } elsif ($time > 108 and $time <= 120) {
        push @time_blocks, "Night 5";
    } elsif ($time > 120 and $time <= 132) {
        push @time_blocks, "Day 6";
    } elsif ($time > 132 and $time <= 144) {
        push @time_blocks, "Night 6";
    } elsif ($time > 144 and $time <= 156) {
        push @time_blocks, "Day 7";
    } else {
        push @time_blocks, "Night 7";
    }
}

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( CBM4_tagging CB05_tagging );
#my @mechanisms = ( "(h) CBM-IV", "(i) CB05" );
my $index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $run);
        ($plot_data{$mechanisms[$index]}, $legend{$mechanisms[$index]}) = get_data($kpp, $mecca, $mechanisms[$index], $mech_species);
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(plyr) `,
        q` library(grid) `,
        q` library(gridExtra) `,
        q` library(dplyr) `,
        q` library(scales) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "CH3O2 + NO" = "#e7e85e", "MEO2 + NO" = "#e7e85e", "MO2 + NO" = "#e7e85e",
                            "CH3CO3 + NO" = "#8b1537", "C2O3 + NO" = "#8b1537", "ACO3 + NO" = "#8b1537",
                            "CRES + OH" = "#cc6329", "CSL + OH" = "#cc6329",
                            "OH + TOL" = "#dc3522",
                            "NO + TLBIPERO2" = "#9bb08f",
                            "MCATEC1O2 + NO" = "#8c6238",
                            "C6H5O2 + NO" = "#4c9383",
                            "NO + TLFUO2" = "#76afca",
                            "MECOACETO2 + NO" = "#f9c600",
                            "MCATEC1O + O3" = "#0352cb", 
                            "C6H5O + O3" = "#86b650", 
                            "NPHEN1O + O3" = "#6c254f", 
                            "Consumption Others" = "#ee6738",
                            "ETHP + NO" = "#f9c600", 
                            "KETP + NO" = "#76afca", 
                            "HC3P + NO" = "#8c6238", 
                            "OH + ONIT" = "#9bb08f", 
                            "ADDC + O3" = "#0352cb", 
                            "CSL + NO3" = "#86b650", "CRES + NO3" = "#86b650",
                            "ADDT + O3" = "#6c254f", 
                            "NO2 + PHO" = "#8ed6d5",
                            "CRO + NO2" = "#898989") `,
);

#plot = plot + scale_y_continuous(limits = c(-25, 15), breaks = seq(-25, 15, 5)) ;
#plot = plot + scale_fill_manual(limits = legend, values = my.colours) ;
$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = time, y = rate, fill = reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, rate > 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, rate < 0), stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 140, angle = 45, vjust = 0.5)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                                            plot = plot + theme(plot.title = element_text(size = 200, face = "bold")) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.text = element_text(size = 140)) ;
                                                            plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                                            return(plot) } `,
); 

$R->set('time', [@time_blocks]);
$R->run(q` plots = list() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [@{$ref->{$reaction}}]);
            $R->run(q` data[reaction] = rate * 1e4 `);
        }
    }
    $R->set('legend', [@{$legend{$run}}]);
    $R->set('mechanism', $run);
    $R->run(q` data = ddply(data, .(time), colwise(sum)) `,
            q` data = data[1:7,] `,
            q` data = melt(data, id.vars = c("time"), variable.name = "reaction", value.name = "rate") `,
            q` reaction.levels = levels(factor(data$reaction)) `,
            q` data$reaction = ordered(data$reaction, levels = reaction.levels) `,
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` CairoPDF(file = "pentane_Ox_intermediates.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(arrangeGrob(plots[[1]] , 
                                                plots[[2]] ,
                                                plots[[3]] ,
                                                plots[[4]] ,
                                                plots[[5]] ,
                                                plots[[6]] ,
                                                plots[[7]] ,
                                                plots[[8]] ,
                                                plots[[9]] ,
                                                nrow = 3), 
                                   nrow = 1, ncol = 1,
                                   left = textGrob(expression(bold(paste("Molecules (intermediate) ", s^-1, "/ Molecules (VOC) x ", 10^4))), rot = 90, gp = gpar(fontsize = 85), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC) = @_;
    $families{"HO2x_${mechanism}"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_${mechanism}", "HO2x_${mechanism}");

    my ($producers, $producer_yields, %production_reaction_rates, $consumers, $consumer_yields, %consumption_reaction_rates);
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
            print "No family found for $species\n";
        }
        print "No producers found for $species\n" if (@$producers == 0);
        print "No consumers found for $species\n" if (@$consumers == 0);
        
        my $prod_max = 5e8;
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($number, $VOC) = split /_/, $reaction;
            next unless (defined $VOC and $VOC eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string;
            if ($reactants =~ /XO2/ and $species =~ /RACM|CB|RADM/) {
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
                    $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
                } 
            } else {
                $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($number, $VOC) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $VOC and $VOC eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string; 
            $consumption_reaction_rates{$reactants} += $rate(1:$NTIME-2); 
        }
    
        remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
        foreach my $reaction (sort keys %production_reaction_rates) {
            if ($production_reaction_rates{$reaction}->sum < $prod_max) {
                $production_reaction_rates{"Production Others"} += $production_reaction_rates{$reaction};
                delete $production_reaction_rates{$reaction};
            }
        }

        foreach my $reaction (sort keys %consumption_reaction_rates) {
            if ($consumption_reaction_rates{$reaction}->sum > -$prod_max) {
                $consumption_reaction_rates{"Consumption Others"} += $consumption_reaction_rates{$reaction};
                delete $consumption_reaction_rates{$reaction};
            }
        }
    }

    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($mechanism =~ /CB/) {
        my $name = "PAR_NC5H12";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
    } else {
        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $production_reaction_rates{$_} /= $parent_emissions foreach (sort keys %production_reaction_rates);
    $consumption_reaction_rates{$_} /= $parent_emissions foreach (sort keys %consumption_reaction_rates);
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    my @sorted_cons = reverse sort { &$sort_function($consumption_reaction_rates{$b}) <=> &$sort_function($consumption_reaction_rates{$a}) } keys %consumption_reaction_rates;

    my @final_sorted_data;
    foreach (@sorted_cons) { #sum up rates of reactions, starting with reaction with lowest sum, consumption others added separately 
        next if ($_ eq 'Consumption Others');
        push @final_sorted_data, { $_ => $consumption_reaction_rates{$_} };
    }

    push @final_sorted_data, { 'Consumption Others' => $consumption_reaction_rates{'Consumption Others'} } if (defined $consumption_reaction_rates{'Consumption Others'}); 
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @final_sorted_data, { $_ => $production_reaction_rates{$_} };
    }

    push @final_sorted_data, { 'Production Others' => $production_reaction_rates{'Production Others'} } if (defined $production_reaction_rates{'Production Others'}); 

    my (@plot_data, @legend_pos, @legend_neg, @legend);
    foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            if ($ref->{$item}->sum > 0) {
                push @legend_pos, $item;
            } else {
                push @legend_neg, $item;
            }
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    } 
    push @legend, reverse @legend_pos;
    push @legend, @legend_neg;
    return (\@plot_data, \@legend);
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
