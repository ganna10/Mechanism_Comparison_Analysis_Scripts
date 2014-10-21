#!/usr/bin/perl -w
# calculate Ox production from tagged and non-tagged model runs
# Version 0: Jane Coates 27/09/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_no_tagging/boxmodel"); 
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

my @runs = qw( MCM_3.2_no_tagging MCM_3.2_tagged );
my @mechanisms = qw( no_tagging tagged );
my $index = 0;

my (%families, %weights, %plot_data);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{$mechanisms[$index]} = ( [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]); 
    $weights{$mechanisms[$index]} = ( { NO3 => 2, N2O5 => 3});
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(plyr) `);
$R->run(q` library(reshape2) `);
$R->run(q` library(grid) `);
$R->run(q` library(gridExtra) `);
$R->run(q` library(RColorBrewer) `);
$R->run(q` library(scales) `);
$R->run(q` library(Cairo) `);

$R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
        q` my.colours = c( "Production Others" = "#696537", 
                           "C2H5O = CH3CHO + HO2" = "#f9c600", 
                           "C2H5O2 + NO = C2H5O + NO2" = "#76afca", 
                           "HCHO + hv = CO + HO2 + HO2" = "#dc3522", 
                           "CH3CO3 + NO = CH3O2 + NO2" = "#8c6238", 
                           "HCHO + OH = CO + HO2" = "#9bb08f", 
                           "CH3O = HCHO + HO2" = "#8b1537", 
                           "CH3O2 + NO = CH3O + NO2" = "#e7e85e", 
                           "CO + OH = HO2" = "#2c9def", 
                           "O3 + OH = HO2" = "#603912", 
                           "IC3H7O2 + NO = IC3H7O + NO2" = "#b569b3",
                           "NC7H16" = "#f9c600", 
                           "EBENZ" = "#76afca", 
                           "BENZENE" = "#dc3522", 
                           "OXYL" = "#8c6238", 
                           "PXYL" = "#9bb08f", 
                           "NC6H14" = "#8b1537", 
                           "IC4H10" = "#e7e85e", 
                           "C3H6" = "#0352cb", 
                           "C2H6" = "#86b650", 
                           "MXYL" = "#6c254f", 
                           "C5H8" = "#ee6738", 
                           "C2H4" = "#58691b", 
                           "NC5H12" = "#8ed6d5", 
                           "C3H8" = "#f3aa7f", 
                           "TOLUENE" = "#c65d6c", 
                           "NC4H10" = "#888a87", 
                           "IC5H12" = "#0e5c28", 
                           "CH4" = "#b569b3") `,
        q` my.names = c("NC7H16" = "Heptane", 
                        "EBENZ" = "Ethylbenzene", 
                        "BENZENE" = "Benzene", 
                        "OXYL" = "o-Xylene", 
                        "PXYL" = "p-Xylene", 
                        "NC6H14" = "Hexane", 
                        "IC4H10" = "2-Methylpropane", 
                        "C3H6" = "Propene", 
                        "C2H6" = "Ethane", 
                        "MXYL" = "m-Xylene", 
                        "C5H8" = "Isoprene", 
                        "C2H4" = "Ethene", 
                        "NC5H12" = "Pentane", 
                        "C3H8" = "Propane", 
                        "TOLUENE" = "Toluene", 
                        "NC4H10" = "Butane", 
                        "IC5H12" = "2-Methylbutane", 
                        "CH4" = "Methane", 
                        "CO + OH = HO2" = "CO",
                        "Production Others" = "Production Others", 
                        "C2H5O = CH3CHO + HO2" = "C2H5O = CH3CHO + HO2", 
                        "C2H5O2 + NO = C2H5O + NO2" = "C2H5O2 + NO = C2H5O + NO2", 
                        "HCHO + hv = CO + HO2 + HO2" = "HCHO + hv = CO + HO2 + HO2", 
                        "CH3CO3 + NO = CH3O2 + NO2" = "CH3CO3 + NO = CH3O2 + NO2", 
                        "HCHO + OH = CO + HO2" = "HCHO + OH = CO + HO2", 
                        "CH3O = HCHO + HO2" = "CH3O = HCHO + HO2", 
                        "CH3O2 + NO = CH3O + NO2" = "CH3O2 + NO = CH3O + NO2") `,
);

$R->run(q` plotting = function (data, title) {   plot = ggplot(data, aes(x = time, y = rate, fill = VOC)) ;
                                                 plot = plot + geom_bar(stat = "identity") ;
                                                 plot = plot + theme_bw() ;
                                                 plot = plot + ggtitle(title) ;
                                                 plot = plot + theme(plot.title = element_text(size = 90, face = "bold")) ;
                                                 plot = plot + theme(axis.text.x = element_text(size = 70, angle = 45, vjust = 0.5)) ;
                                                 plot = plot + theme(axis.text.y = element_text(size = 60)) ;
                                                 plot = plot + theme(axis.title.x = element_blank()) ;
                                                 plot = plot + theme(axis.title.y = element_blank()) ;
                                                 plot = plot + theme(legend.key.size = unit(4, "cm")) ;
                                                 plot = plot + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) ;
                                                 plot = plot + theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) ;
                                                 plot = plot + theme(legend.text = element_text(size = 40)) ;
                                                 plot = plot + theme(legend.title = element_blank()) ;
                                                 plot = plot + theme(legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99)) ;
                                                 plot = plot + scale_fill_manual( limits = rev(VOC.levels), labels = my.names, values = my.colours) ;
                                                 plot = plot + scale_y_continuous(limits=c(0, 1.5e9), breaks=seq(0, 1.5e9, 2e8), label = scientific_10);
                                                 plot = plot + theme(legend.key = element_blank()) ;
                                                 return(plot) } `,
);

$R->set('time', [@time_blocks]);
$R->run(q` plots = list() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` data[name] = rate`); 
        }
    }
    $R->run(q` data = ddply(data, .(time), colwise(sum)) `,
            q` data = melt(data, id.vars = c("time"), variable.name = "VOC", value.name = "rate") `,
            q` VOC.levels = levels(factor(data$VOC)) `,
            q` data$VOC = ordered(data$VOC, levels = VOC.levels) `,
            q` plot = plotting(data, "(a) No Tagging") `,
            q` plots = c(plots, list(plot)) `,
    );
    #my $p = $R->run(q` print(data) `);
    #print $p, "\n";
} 

$R->run(q` CairoPDF(file = "MCMv3_2_Ox_no_tagging_tagging_budget.pdf", width = 57, height = 40) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] ,
                                                    plots[[2]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop(); 

sub get_data {
    my ($mecca, $kpp, $Ox) = @_;
    $families{'HO2x'} = ( [ qw( HO2 HO2NO2 ) ] );
    my @species = ($Ox, "HO2x");

    my (%production_reaction_rates, %species_production_rates, %consumption_reaction_rates, %species_consumption_rates); 
    foreach my $species (@species) { #get all production and consumption rates
        my ($producers, $producer_yields, $consumers, $consumer_yields);
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
            print "No family for $species found\n";
        }

        die "No producers found for $species\n" if (@$producers == 0);
        die "No consumers found for $species\n" if (@$consumers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string = $kpp->reaction_string($reaction);
            if (defined $parent) { # for tagged reactions
                #$string =~ s/_$parent//g; #removing tag from reaction strings
                #$species_production_rates{$species}{$parent}{$string} += $rate;
                $string = $parent; # in order to merge all production rates from all parent species reactions into 1 pdl
            }
            $production_reaction_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string = $kpp->reaction_string($reaction);
            if (defined $parent) { # for tagged reactions
                #$string =~ s/_$parent//g; #removing tag from reaction strings
                #$species_consumption_rates{$species}{$parent}{$string} += $rate;
                $string = $parent; # in order to merge all consumption rates from all parent species reactions into 1 pdl
            }
            $consumption_reaction_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    } 
    remove_common_processes($production_reaction_rates{'HO2x'}, $consumption_reaction_rates{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production_reaction_rates{'HO2x'}{$_} for (keys %{ $production_reaction_rates{'HO2x'} });

    foreach my $reaction( keys %{ $production_reaction_rates{'HO2x'} }) {
        $production_reaction_rates{$Ox}{$reaction} += $production_reaction_rates{$Ox}{'HO2 + NO = NO2 + OH'} * $production_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption_reaction_rates{$Ox}{$reaction} += $consumption_reaction_rates{$Ox}{'HO2 + O3 = OH'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption_reaction_rates{$Ox}{$reaction} += $consumption_reaction_rates{$Ox}{'HO2 + NO3 = NO2 + OH'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
    }
    delete $production_reaction_rates{$Ox}{'HO2 + NO = NO2 + OH'};
    delete $consumption_reaction_rates{$Ox}{'HO2 + O3 = OH'};
    delete $consumption_reaction_rates{$Ox}{'HO2 + NO3 = NO2 + OH'}; 
    remove_common_processes($production_reaction_rates{$Ox}, $consumption_reaction_rates{$Ox});

    my $prod_others_max = 8e7;
    my $sort_function = sub { $_[0]->sum };
    foreach my $item (keys %{$production_reaction_rates{$Ox}}) {
        if ($production_reaction_rates{$Ox}{$item}->sum < $prod_others_max) { #get production others
            $production_reaction_rates{$Ox}{'Production Others'} += $production_reaction_rates{$Ox}{$item};
            delete $production_reaction_rates{$Ox}{$item};
        }
    }

    #$prod_hash{$_}->where($prod_hash{$_} < 0) .= 0 foreach (keys %prod_hash); 
    my @sorted_prod = sort { &$sort_function($production_reaction_rates{$Ox}{$b}) <=> &$sort_function($production_reaction_rates{$Ox}{$a}) } keys %{$production_reaction_rates{$Ox}};

    my (@sorted_plot_data);
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production_reaction_rates{$Ox}{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production_reaction_rates{$Ox}{'Production Others'} } if (defined $production_reaction_rates{$Ox}{'Production Others'}); #add Production Others to the beginning 

    my @plot_data;
    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    } 
    return \@plot_data;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, "<:encoding(utf-8)", $file or die $!; 
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
