#! /usr/bin/env perl
# plot reactions contribution to Ox production from toluene degradation in MCM v3.2, CBM4 and CB05 mechanisms
# Version 0: Jane Coates 25/9/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
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

my $NMVOC = "TOLUENE";
my @runs = qw( MCM_3.2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) CBM-IV", "(c) CB05" );
#my @runs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my $index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3};
    ($plot_data{$mechanisms[$index]}, $legend{$mechanisms[$index]}) = get_data($NMVOC, $mecca, $kpp, "Ox_$mechanisms[$index]");
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "CH3O2 + NO = CH3O + NO2" = "#8c1531", "MEO2 + NO = FORM + HO2 + NO2" = "#8c1531",
                            "CH3CO3 + NO = CH3O2 + NO2" = "#1c3e3d", "C2O3 + NO = HCHO + HO2 + NO2 + XO2" = "#1c3e3d", "C2O3 + NO = MEO2 + NO2" = "#1c3e3d",
                            "NO + XO2 = NO2" = "#ef6638",
                            "NO + TLBIPERO2 = NO2 + TLBIPERO" = "#8fd5d3",
                            "MCATEC1O2 + NO = MCATEC1O + NO2" = "#c9a415",
                            "C6H5O2 + NO = C6H5O + NO2" = "#b569b3",
                            "NO + TLFUO2 = NO2 + TLFUO" = "#2c9dad",
                            "MECOACETO2 + NO = MECOACETO + NO2" = "#0b4074") `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = time, y = rate, fill = reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 70)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 60)) ;
                                                            plot = plot + theme(plot.title = element_text(size = 90, face = "bold")) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.text = element_text(size = 60)) ;
                                                            plot = plot + theme(legend.key.size = unit(6, "cm")) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 7e7), breaks = seq(0, 7e7, 1e7)) ;
                                                            plot = plot + scale_fill_manual(limits = legend, values = my.colours) ;
                                                            return(plot) } `,
); 

$R->set('time', [@time_blocks]);
$R->run(q` plots = list() `);
foreach my $run (sort keys %plot_data) {
    $R->set('mechanism', $run);
    $R->run(q` pre = data.frame(time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [@{$ref->{$reaction}}]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('legend', [@{$legend{$run}}]);
    $R->run(q` pre = ddply(pre, .(time), colwise(sum)) `,
            q` pre = pre[1:7,] `,
            q` pre = melt(pre, id.vars = c("time"), variable.name = "reaction", value.name = "rate") `,
            q` pre$Mechanism = rep(mechanism, length(pre$time)) `,
            q` reaction.levels = levels(factor(pre$reaction)) `,
            q` pre$reaction = ordered(pre$reaction, levels = reaction.levels) `,
            q` plot = plotting(pre, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` CairoPDF(file = "CBMs_MCM_TOL_Ox_breakdown.pdf", width = 100, height = 71) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] , 
                                                    plots[[2]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[3]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 80), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($parent, $mecca, $kpp, $species) = @_;
    
    my (%production_reaction_rates, %consumption_reaction_rates, $producers, $producer_yields, $consumers, $consumer_yields);
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
    
    my $prod_max = 5e6;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($number, $VOC) = split /_/, $reaction;
        next unless (defined $VOC and $VOC eq $parent);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_$VOC\b//g;
            if ($rate->sum < $prod_max) {
                $production_reaction_rates{'Production Others'} += $rate(1:$NTIME-2);
            } else {
                $production_reaction_rates{$reaction_string} += $rate(1:$NTIME-2); 
            }
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my ($number, $VOC) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $VOC and $VOC eq $parent);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_$VOC\b//g;
        $consumption_reaction_rates{$reaction_string} += $rate(1:$NTIME-2); 
    }

    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;

    my @final_sorted_data;
    foreach (@sorted_data) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
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
