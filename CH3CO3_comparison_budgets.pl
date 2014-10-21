#! /usr/bin/env perl
# Plot CH3CO3 production and comparison budgets for MCM3.2 and MOZART mechanisms
# Version 0: Jane Coates 28/9/2014

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
foreach my $time (@time_axis) {#map to day and night
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

my @runs = qw( MCM_3.2_tagged MOZART_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MOZART-4" );
my @base_name = qw( CH3CO3 CH3CO3 );
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
#my @base_name = qw( ACO3 );
my $array_index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $spcfile = "$base/$run/gas.spc";
    my $all_tagged_species = get_tagged_species($base_name[$array_index], $spcfile); 
    $families{$mechanisms[$array_index]} = [ @$all_tagged_species ];
    ($plot_data{$mechanisms[$array_index]}, $legend{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index]);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(RColorBrewer) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` arrange.data = function (data) { data = ddply(data, .(Time), colwise(sum)) ;
                                            data = data[1:7,] ;
                                            data = melt(data, id.vars = c("Time"), variable.name = "Reaction", value.name = "Rate") ;
                                            return(data) } `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "Consumption Others" = "#1c3e3d",
                            "MEKO2 + NO = CH3CHO\n+ CH3CO3 + NO2" = "#0e5c28",
                            "CH3CO3 + NO = CH3O2 + NO2" = "#4c9383", "CH3CO3 + NO = CH3O2\n+ CO2 + NO2" = "#4c9383", 
                            "PAN = CH3CO3 + NO2" = "#8d1435", 
                            "CH3CO3 + NO2 = PAN" = "#6d6537", 
                            "CH3CHO + OH = CH3CO3" = "#e7e85e", "CH3CHO + OH = CH3CO3\n+ H2O" = "#e7e85e" ) `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, Rate >= 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") ;
                                                            plot = plot + scale_y_continuous(limits = c(-5e8, 5e8), breaks = seq(-5e8, 5e8, 1e8)) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(plot.title = element_text(size = 90, face = "bold")) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 70, angle = 45, vjust = 0.5)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 60)) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.text = element_text(size = 40)) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.key.size = unit(4, "cm")) ;
                                                            plot = plot + scale_fill_manual(limits = legend, values = my.colours) ;
                                                            return(plot) } `,
);

$R->set('Time', [@time_blocks]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [@{$ref->{$reaction}}]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('legend', [@{$legend{$run}}]);
    $R->set('mechanism', $run);
    $R->run(q` data = arrange.data(data) `,
            q` reaction.levels = levels(factor(data$Reaction)) `,
            q` data$Reaction = ordered(data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `, #add plot to list 
    );
}

$R->run(q` CairoPDF(file = "CH3CO3_budget_comparison.pdf", width = 57, height = 40) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] ,
                                                    plots[[2]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 85), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $species) = @_;

    my ($consumers, $producers, $consumer_yields, $producer_yields, %production_reaction_rates, %consumption_reaction_rates);
    if (exists $families{$species}) { 
        $kpp->family({ 
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $consumers = $kpp->consuming($species);
        $producers = $kpp->producing($species);
        $consumer_yields = $kpp->effect_on($species, $consumers); 
        $producer_yields = $kpp->effect_on($species, $producers);  
    } else {
        print "No family found for $species\n";
    }

    #check that species reactions are found
    die "No producers found for $species\n" if (@$producers == 0);
    die "No consumers found for $species\n" if (@$consumers == 0);

    my $prod_others_max = 4e7;
    my $max_string_width = 28;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        if ($reaction_string =~ /MEKO2 \+ NO/) {
            $reaction_string = "MEKO2 + NO = CH3CHO\n+ CH3CO3 + NO2";
        } elsif ($reaction_string =~ /CH3CO3 \+ H2O/) {
            $reaction_string = "CH3CHO + OH = CH3CO3\n+ H2O";
        }
        if ($rate->sum <= $prod_others_max) {
            $production_reaction_rates{"Production Others"} += $rate;
        } else {
            $production_reaction_rates{$reaction_string} += $rate;
        }
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        if ($reaction_string =~ /CH3CO3 \+ NO = CH3O2 \+ CO2/) {
            $reaction_string = "CH3CO3 + NO = CH3O2\n+ CO2 + NO2";
        }
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 

    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    
    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    my @cons_sorted_data = reverse sort { &$sort_function($consumption_reaction_rates{$b}) <=> &$sort_function($consumption_reaction_rates{$a}) } keys %consumption_reaction_rates;

    my @final_sorted_data;
    foreach (@cons_sorted_data) { 
        next if ($_ eq 'Consumption Others') ;
        push @final_sorted_data, { $_ => $consumption_reaction_rates{$_} };
    } 
    push @final_sorted_data, {'Consumption Others' => $consumption_reaction_rates{'Consumption Others'}} if (defined $consumption_reaction_rates{'Consumption Others'}) ;

    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Production Others') ;
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
            my @new_array = splice @rate_array, 1, 504;
            push @plot_data, { $item => \@new_array };
        }
    } 
    push @legend, reverse @legend_pos;
    push @legend, @legend_neg;
    return (\@plot_data, \@legend);
}

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

sub get_tagged_species {
    my ($species, $spc_file) = @_;
    my $all_species = read_file($spc_file);
    my @tagged_species;
    foreach my $line (@$all_species) {
        next unless ($line =~ /^$species/);
        $line =~ s/\s=.*$//;
        push @tagged_species, $line;
    }
    return \@tagged_species;
}

sub read_file {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file for reading : $!";
    chomp(my @all = <$in>);
    close $in;
    return \@all;
} 
