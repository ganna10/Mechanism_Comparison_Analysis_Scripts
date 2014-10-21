#!/usr/bin/perl
# radical budget plot attributed back to the parent VOCs - only with net effects of organic functional groups, facetted by functional group
# Version 0 : Jane Coates 26/06/2014 
# Version 1 : Jane Coates 29/08/2014 including inorganic radicals in radical family, re-factoring code 

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

#Create x-axis for plot in hours
my $run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mecca = MECCA->new($run); 
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

my (%families, %weights, %plot_data, %legend);
my $base = "/work/users/jco/MECCA";
#my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
#my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my @runs = qw( RADM2_tagged ) ;
my @mechanisms = qw( RADM2 );
my $array_index = 0;

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    my $pan_file = "$base/$run/All_pans.txt";
    my @pans = get_species($pan_file);
    #$families{$mechanisms[$array_index]} = [ @radicals ];
    $families{$mechanisms[$array_index]} = [ @radicals, "PAN_OL2", "PAN_ETH", "PAN_HC3", "PAN_HC5", "PAN_HC8", "PAN_ISO", "PAN_OLI", "PAN_OLT", "PAN_TOL", "PAN_XYL" ];
#    $weights{$mechanisms[$array_index]} = { "H2O2" => 2 };
    ($plot_data{$mechanisms[$array_index]}, $legend{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index]);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(plyr) `,
        q` library(grid) `,
        q` library(gridExtra) `,
        q` library(RColorBrewer) `,
);

$R->run(q` arrange.data = function (data) { data = ddply(data, .(time), colwise(sum)) ;
                                            data["PAN Extra Consumption"] = c(-3359245.9, -7939716.2, -6719866.4, -4616477.4, -2890404.2, -1778101.0, -1121882.7, -5589058.1, -5587978.2, -4039766.6, -2563523.8, -1553588.9, -951734.1, -606696.6) ;
                                            data = melt(data, id.vars = c("time"), variable.name = "Reaction", value.name = "Rate") ;
                                            return(data) } `,
);
#
#                                                            plot = plot + theme(plot.title = element_text(size = 50, face = "bold")) ;
#                                                            plot = plot + theme(axis.text.x = element_text(size = 30)) ;
#                                                            plot = plot + theme(axis.text.y = element_text(size = 30)) ;
#                                                            plot = plot + theme(axis.title.y = element_text(size = 55, face = "bold")) ;
#                                                            plot = plot + theme(legend.title = element_text(size = 50, face = "bold")) ;
#                                                            plot = plot + theme(legend.text = element_text(size = 25)) ;
#                                                            plot = plot + theme(legend.key = element_blank()) ;
#                                                            plot = plot + theme(legend.key.size = unit(2, "cm")) ;
#                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
#                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
$R->run(q` plotting = function (data, legend, mechanism) {  cols = colorRampPalette(brewer.pal(8,"Dark2"))(nlevels(data$Reaction)) ; 
                                                            plot = ggplot(data, aes(x = time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, Rate >= 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + xlab("\n") ;
                                                            plot = plot + ylab(expression(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))) ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + scale_fill_manual(values = cols) ;
                                                            return(plot) } `,
);
                                                            #plot = plot + scale_y_continuous(limits = c(-8e8, 8e8), breaks = seq(-8e8, 8e8, 2e8)) ;

$R->set('time', [@time_blocks]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` data = data.frame(time) `);
    foreach my $hash_ref (@{$plot_data{$mechanism}}) {
        foreach my $reaction (keys %$hash_ref) {
            $R->set('Reaction', $reaction);
            $R->set('Rate', [@{$hash_ref->{$reaction}}]);
            $R->run(q` data[Reaction] = Rate `);
        }
    }
    $R->set('legend', [@{$legend{$mechanism}}]);
    $R->set('mechanism', $mechanism);
    $R->run(q` data = arrange.data(data) `,
            q` reaction.levels = levels(factor(data$Reaction)) `,
            q` data$Reaction = ordered(data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `, #add plot to list 
    );
}

$R->run(q` data = ddply(data, .(time), colwise(sum)) `,
        q` data["PAN Extra Consumption"] = c(-3359245.9, -7939716.2, -6719866.4, -4616477.4, -2890404.2, -1778101.0, -1121882.7, -5589058.1, -5587978.2, -4039766.6, -2563523.8, -1553588.9, -951734.1, -606696.6) `,
);


#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(
    #q` CairoPDF(file = "new_radicals_budget_assigned_without_PAN.pdf", width = 10, height = 7) `,
        q` CairoPDF(file = "new_radicals_budget_assigned_with_PAN.pdf", width = 10, height = 7) `,
        q` args.list = c(plots, list(nrow = 3)) `,
        q` multiplot = do.call(grid.arrange, plots) `,
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

    my $prod_others_max = 2e7;
    my $max_string_width = 35;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        $reaction_string =~ s/\n$//;
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
        $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        $reaction_string =~ s/\n$//;
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 
    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    ##net budget analysis
    my ($total_production, $total_consumption);
    $total_production += $production_reaction_rates{$_} foreach (keys %production_reaction_rates);
    $total_consumption += $consumption_reaction_rates{$_} foreach (keys %consumption_reaction_rates);
    my $net = $total_production + $total_consumption;
    $net = $net(1:$NTIME-2);
    my @net_rates = map {$_} $net->dog;
    my $R = Statistics::R->new();
    $R->run(q` library(plyr) `);
    $R->set('time', [@time_blocks]);
    $R->set('net', [@net_rates]);
    $R->run(q` data = data.frame(time, net) `,
            q` data = ddply(data, .(time), colwise(sum)) `,
    );
    my $p = $R->run(q` print(data) `);
    print $p, "\n";
    $R->stop;
    
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

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    push @ro2, split /\s+/, $_ for (<$in>) ;
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            push @no2_reservoirs, $products->[0] if (@$products == 1) ;
        }   
    }   
    return @no2_reservoirs;
} 

sub get_species { #get radicals from file
    my ($file) = @_; 

    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my %hash;
    my @lines = <$in>;
    close $in;
    my @separate = map { split /\s/, $_ } @lines;
    $hash{$_} += 1 foreach (@separate);
    return keys %hash;
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
