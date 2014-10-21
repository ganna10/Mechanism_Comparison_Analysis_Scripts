#!/usr/bin/perl
# calculate differences between rates when including PAN in radical family
# Version 0 : Jane Coates 15/9/2014

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

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{"withPAN"} = [ @radicals, "PAN_OL2", "PAN_ETH", "PAN_HC3", "PAN_HC5", "PAN_HC8", "PAN_ISO", "PAN_OLI", "PAN_OLT", "PAN_TOL", "PAN_XYL" ];
    $families{"withoutPAN"} = [ @radicals ];
#    $weights{"withPAN"} = { "H2O2" => 2 };
    ($plot_data{"withPAN"}, $legend{"withPAN"}) = get_data($kpp, $mecca, "withPAN");
    ($plot_data{"withoutPAN"}, $legend{"withoutPAN"}) = get_data($kpp, $mecca, "withoutPAN");
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
                                            data = melt(data, id.vars = c("time"), variable.name = "Reaction", value.name = "Rate") ;
                                            return(data) } `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  cols = colorRampPalette(brewer.pal(8,"Dark2"))(nlevels(data$Reaction)) ; 
                                                            plot = ggplot(data, aes(x = time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, Rate >= 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + xlab("\n") ;
                                                            plot = plot + ylab(expression(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))) ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(plot.title = element_text(size = 50, face = "bold")) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 30)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 30)) ;
                                                            plot = plot + theme(axis.title.y = element_text(size = 55, face = "bold")) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.title = element_text(size = 50, face = "bold")) ;
                                                            plot = plot + theme(legend.text = element_text(size = 25)) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.key.size = unit(2, "cm")) ;
                                                            plot = plot + scale_fill_manual(values = cols) ;
                                                            return(plot) } `,
);
                                                            #plot = plot + scale_y_continuous(limits = c(-8e8, 8e8), breaks = seq(-8e8, 8e8, 2e8)) ;

$R->set('time', [@time_blocks]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` plot.data = data.frame(time) `);
    foreach my $hash_ref (@{$plot_data{$mechanism}}) {
        foreach my $reaction (keys %$hash_ref) {
            $R->set('Reaction', $reaction);
            $R->set('Rate', [@{$hash_ref->{$reaction}}]);
            $R->run(q` plot.data[Reaction] = Rate `);
        }
    }
    $R->set('legend', [@{$legend{$mechanism}}]);
    $R->set('mechanism', $mechanism);
    $R->run(q` plot.data = ddply(plot.data, .(time), colwise(sum)) `);
    if ($mechanism eq 'withPAN') {
        $R->run(q` plot.data["ACO3 + NO2 = PAN"] = c( -2.83e+08, -4.83e+08, -3.76e+08, -2.48e+08, -1.53e+08, -9.36e+07, -5.91e+07, -4.60e+08, -4.54e+08, -3.28e+08, -2.08e+08, -1.27e+08, -7.80e+07, -4.98e+07) `,
                q` plot.data["PAN = ACO3 + NO2"] = c( 2.15e+08, 4.69e+08, 3.78e+08, 2.54e+08, 1.58e+08, 9.70e+07, 6.11e+07, 4.59e+08, 4.58e+08, 3.31e+08, 2.10e+08, 1.27e+08, 7.81e+07, 4.98e+07) `,
        );
    }
    $R->run(q` plot.data = melt(plot.data, id.vars = c("time"), variable.name = "Reaction", value.name = "Rate") `);
    $R->run(#q` plot.data = arrange.data(plot.data) `,
            q` reaction.levels = levels(factor(plot.data$Reaction)) `,
            q` plot.data$Reaction = ordered(plot.data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(plot.data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `, #add plot to list 
    );
}

#my $p = $R->run(q` print(head(plot.data)) `);
#print "$p\n";

$R->run(q` CairoPDF(file = "difference_plots.pdf", width = 100, height = 100) `,
        q` args.list = c(plots, list(ncol = 2)) `,
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

    my $prod_others_max = 1e0;
    my $max_string_width = 35;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        #$reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        #$reaction_string =~ s/\n$//;
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
        #$reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        #$reaction_string =~ s/\n$//;
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 
    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    my ($total_p, $total_c);
    $total_p += $production_reaction_rates{$_} foreach (keys %production_reaction_rates);
    $total_p = $total_p(1:$NTIME-2);
    print "$species: Total prod = ", $total_p->sum, "\n";
    $total_c += $consumption_reaction_rates{$_} foreach (keys %consumption_reaction_rates);
    $total_c = $total_c(1:$NTIME-2);
    print "$species: Total cons = ", $total_c->sum, "\n";
    
    my $R = Statistics::R->new();
    $R->run(q` library(reshape2) `,
            q` library(plyr) `,
    );
    $R->set('time', [@time_blocks]);
    $R->run(q` prod = data.frame(time) `);
    $R->run(q` cons = data.frame(time) `);
    foreach my $reaction (keys %production_reaction_rates) {
        my @rates = map { $_ } $production_reaction_rates{$reaction}->dog;
        @rates = splice @rates, 1, 504;
        $R->set('Reaction', $reaction);
        $R->set('Rate', [@rates]);
        $R->run(q` prod[Reaction] = Rate `);
    }

    foreach my $reaction (keys %consumption_reaction_rates) {
        my @rates = map { $_ } $consumption_reaction_rates{$reaction}->dog;
        @rates = splice @rates, 1, 504;
        $R->set('Reaction', $reaction);
        $R->set('Rate', [@rates]);
        $R->run(q` cons[Reaction] = Rate `);
    }

    $R->run(q` prod = ddply(prod, .(time), colwise(sum)) `,
            q` cons = ddply(cons, .(time), colwise(sum)) `,
    );

    $R->run(q` options(scipen = 0) `, #force scientific notation
            q` options(digits = 3) `,
    );
#    my $p = $R->run(q` print(cat(rowSums(prod[,-1]), sep="\n")) `);
#    if ($species eq "withoutPAN") {
#        my $p = $R->run(q` print(cons["ACO3 + NO2 = PAN"]) `);
#        print "$species \n$p\n";
#    #    my $p1 = $R->run(q` print(cat(rowSums(cons[,-1]), sep="\n")) `);
#        my $p1 = $R->run(q` print(prod["PAN = ACO3 + NO2"]) `);
#        print "$species \n$p1\n";
#    }

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
