#! /usr/bin/env perl
# plot reactions contribution to Ox production from toluene degradation in MCM v3.2, RACM2, CBM4 and CB05 mechanisms
# Version 0: Jane Coates 25/9/2014
# Version 1: Jane Coates 21/10/2014 including assignment of XO2 producing reactions to budget
# Version 2: Jane Coates 10/11/2014 including HO2x production and updates to script
# Version 3: Jane Coates 26/12/2014 including all mechanisms in analysis and re-factoring code

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#MCM emissions used for TOPP calculation of de-lumped VOC
my $kpp = KPP->new("$base/MCMv3.2_tagged/gas.eqn");
my $emission_reaction = $kpp->producing_from("TOLUENE", "UNITY");
next if (@$emission_reaction == 0);
my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
my $emission_rate = $mecca->rate($reaction_number);
$emission_rate = $emission_rate(1:$NTIME-2);
my $mcm_emission = $emission_rate->sum * $dt;

my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my (%families, %weights, %plot_data, %legend);
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    ($plot_data{$mechanism}, $legend{$mechanism}) = get_data($mecca, $kpp, "Ox_$mechanism");
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "CH3O2 + NO" = "#f9c500", "MEO2 + NO" = "#f9c500", "MO2 + NO" = "#f9c500",
                            "CH3CO3 + NO" = "#8b1537", "C2O3 + NO" = "#8b1537", "ACO3 + NO" = "#8b1537",
                            "CRES + OH" = "#cc6329", "CSL + OH" = "#cc6329",
                            "OH + TOL" = "#1c3e3d",
                            "NO + TLBIPERO2" = "#b569b3", "NO + RA16O2" = "#b569b3", "NO + TOLO2" = "#b569b3", "NO + TOLP" = "#b569b3",
                            "OH + TOLUENE" = "#8ed6dc",
                            "NO + TCO3" = "#6c254f",
                            "TO2" = "#f3aa7f",
                            "TLBIPERO" = "#4c9383",
                            "HCHO + OH" = "#ae4901", "CH2O + OH" = "#ae4901",
                            "GLYOXAL + OH" = "#a67c52",
                            "CH3O" = "#0e5c28",
                            "CO + OH" = "#77aecc", 
                            "Consumption Others" = "#ee6738",
                            "KETP + NO" = "#76afca", 
                            "HC3P + NO" = "#8c6238", 
                            "OH + ONIT" = "#9bb08f", 
                            "ADDC + O3" = "#0352cb", 
                            "CSL + NO3" = "#86b650", "CRES + NO3" = "#86b650",
                            "ADDT + O3" = "#6c254f", 
                            "NO2 + PHO" = "#8ed6d5",
                            "C2H5O2 + NO" = "#c9a415", "ETHP + NO" = "#c9a415",
                            "HOCH2CH2O2 + NO" = "#e7e85e",
                            "NO + RN10O2" = "#1c3e3d",
                            "BIGALD + hv" = "#000000", 
                            "CRO + NO2" = "#898989") `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = time, y = rate, fill = reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, rate > 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, rate < 0), stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(axis.title = element_blank()) ;
                                                            plot = plot + scale_y_continuous(limits = c(-35, 35), breaks = seq(-35, 35, 5), expand = c(0, 0.0)) ;
                                                            plot = plot + scale_x_discrete(expand = c(0, 0)) ;
                                                            plot = plot + theme(axis.text.x = element_text(face = "bold", size = 20, angle = 45, hjust = 0.8, vjust = 0.7)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 18)) ;
                                                            plot = plot + theme(plot.title = element_text(size = 22, face = "bold")) ;
                                                            plot = plot + theme(panel.grid = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(panel.border = element_rect(colour = "black")) ;
                                                            plot = plot + theme(legend.justification = c(1.0, 0.0)) ;
                                                            plot = plot + theme(legend.position = c(1.0, 0.0)) ;
                                                            plot = plot + scale_fill_manual(values = my.colours, limits = legend) ;
                                                            return(plot) } `,
); 

$R->set('time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
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
    $R->run(q` data = gather(data, reaction, rate, -time ) `,
            q` reaction.levels = levels(factor(data$reaction)) `,
            q` data$reaction = ordered(data$reaction, levels = reaction.levels) `,
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` CairoPDF(file = "TOL_Ox_intermediates.pdf", width = 15.6, height = 22.0) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[9]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                                    plots[[7]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[8]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[6]] + theme(axis.title.y = element_blank()),
                                                    plots[[2]] + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[1]] + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob("\nMolecules (intermediate) s-1 / Molecules (VOC)", gp = gpar(fontface = "bold", fontsize = 26), rot = 90, vjust = 0.5) ) `, 
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $Ox) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("HO2x", $Ox);
    
    my (%production_reaction_rates, %consumption_reaction_rates, $producers, $producer_yields, $consumers, $consumer_yields);
    foreach my $species (@loop) {
        if (exists $families{$species}) { #get family reaction numbers and yields
            $kpp->family({ 
                    name    => $species,
                    members => $families{$species},
                    weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);  
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);  
        } else { #get reaction numbers and yields
            print "No family found for $species\n";
        }

        die "No producers found for $species\n" if (@$producers == 0);
        die "No consumers found for $species\n" if (@$consumers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($number, $VOC) = split /_/, $reaction;
            next unless (defined $VOC and $VOC =~ /TOL/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string;
            if ($reactants =~ /XO2/ and $species =~ /RA|CB/) {
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
            next unless (defined $VOC and $VOC =~ /TOL/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string; 
            $consumption_reaction_rates{$reactants} += $rate(1:$NTIME-2); 
        }
    }

    my $prod_max = 2e7;
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

    my $emission_rate;
    if ($Ox =~ /RA|MOZ/) {
        $emission_rate = $mcm_emission;
    } else {
        my $parent;
        if ($Ox =~ /CB/) {
            $parent = "TOL_TOLUENE";
        } else {
            $parent = "TOLUENE";
        }
        my $emission_reaction = $kpp->producing_from($parent, "UNITY");
        my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
        $emission_rate = $mecca->rate($reaction_number); 
        $emission_rate = $emission_rate(1:$NTIME-2);
        $emission_rate = $emission_rate->sum * $dt; 
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $production_reaction_rates{$_} /= $emission_rate foreach (sort keys %production_reaction_rates);
    $consumption_reaction_rates{$_} /= $emission_rate foreach (sort keys %consumption_reaction_rates);

    foreach my $reaction (keys %production_reaction_rates) {
        if ($Ox =~ /MOZ/) {
            $production_reaction_rates{$reaction} *= 0.478;
        } elsif ($Ox =~ /RADM2/ or $Ox =~ /RACM\b/) {
            $production_reaction_rates{$reaction} *= 0.667;
        } elsif ($Ox =~ /RACM2/) {
            $production_reaction_rates{$reaction} *= 0.868;
        }
        my $reshape = $production_reaction_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_reaction_rates{$reaction} = $integrate;
    }

    foreach my $reaction (keys %consumption_reaction_rates) {
        my $reshape = $consumption_reaction_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $consumption_reaction_rates{$reaction} = $integrate;
    }
    
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
    open FILE, $file or die $!; 
    my @ro2;
    for (<FILE>) {
        push @ro2, split /\s+/, $_; 
    }
    close FILE;
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
