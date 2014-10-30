#!/usr/bin/perl
# radical budget plot attributed back to the parent VOCs - only with net effects of organic functional groups, facetted by functional group
# Version 0 : Jane Coates 26/06/2014 
# Version 1 : Jane Coates 29/08/2014 including inorganic radicals in radical family, re-factoring code 
# Version 2 : Jane Coates 29/9/2014 including only production and day-time reactions

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
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
my $array_index = 0;

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{$mechanisms[$array_index]} = [ @radicals ];
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
);

$R->run(q` arrange.data = function (data) { data = ddply(data, .(time), colwise(sum)) ;
                                            data = data[1:7,] ;
                                            data = melt(data, id.vars = c("time"), variable.name = "Reaction", value.name = "Rate") ;
                                            return(data) } `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "O1D = OH + OH" = "#f7c56c",
                            "HONO + hv = NO + OH" = "#6c254f", "HONO + hv = OH + NO" = "#6c254f",
                            "HO2NO2 = HO2 + NO2" = "#f3aa7f",
                            "HO2 + NO2 = HO2NO2" = "#0e5c28",
                            "NO2 + OH = HNO3" = "#6db875", "OH + NO2 = HNO3" = "#6db875",
                            "NO + OH = HONO" = "#f9c500", "OH + NO = HONO" = "#f9c500",
                            "Consumption Others" = "#1c3e3d",
                            "HO2 + HO2 = H2O2" = "#a67c52",
                            "HO2 + OH = UNITY" = "#86b650", "OH + HO2 = UNITY" = "#86b650",
                            "CH3O2 + NO2 = CH3O2NO2" = "#77aecc", 
                            "CH3O2NO2 = CH3O2 + NO2" = "#8c1531",
                            "HCHO + hv = CO + 2 HO2" = "#dc3522", "CH2O + hv = CO + 2 HO2" = "#dc3522", "FORM + hv = CO + 2 HO2" = "#dc3522", "HCHO + hv = CO + 2.0 HO2" = "#dc3522", 
                            "PAN = CH3CO3 + NO2" = "#0c3f78", "PAN = C2O3 + NO2" = "#0c3f78", "PAN = ACO3 + NO2" = "#0c3f78",
                            "CH3CO3 + NO2 = PAN" = "#ef6638", "C2O3 + NO2 = PAN" = "#ef6638", "ACO3 + NO2 = PAN" = "#ef6638",
                            "PPN = C2H5CO3 + NO2" = "#c9a415", "PPN = NO2 + RCO3" = "#c9a415", "PANX = CXO3 + NO2" = "#c9a415",
                            "C2H5CO3 + NO2 = PPN" = "#58691b", "NO2 + RCO3 = PPN" = "#58691b", "CXO3 + NO2 = PANX" = "#58691b",
                            "TPAN = NO2 + TCO3" = "#4c9383", 
                            "NO2 + TCO3 = TPAN" = "#cc6638",
                            "CH4 + OH = HCHO + HO2 + XO2" = "#898989",
                            "NO + XO2 = NO2" = "#6d6537",
                            "OH + PAR = 0.06 ALD2 + 0.05 ALDX\n+ 0.11 HO2 + 0.11 PAROP + 0.76 ROR\n+ 0.87 XO2 + 0.13 XO2N" = "#e7e85e", "OH + PAR = 0.11 ALD2 + 0.11 HO2\n+ 0.11 PAROP + 0.76 ROR\n+ 0.87 XO2 + 0.13 XO2N" = "#e7e85e" ) `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 2.1e9), breaks = seq(0, 2.1e9, 5e8)) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(plot.title = element_text(size = 200, face = "bold")) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 140, angle = 45, vjust = 0.5)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.text = element_text(size = 140)) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + scale_fill_manual(limits = legend, values = my.colours) ;
                                                            return(plot) } `,
);

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

#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` CairoPDF(file = "radicals_budget_assigned.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(   arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
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
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
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

    my $prod_others_max = 7.2e7;
    my $max_string_width = 27;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        $reaction_string =~ s/\n$//;
        $reaction_string = string_mapping($reaction_string);
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
        $reaction_string = string_mapping($reaction_string);
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 

    sub string_mapping {
        my ($string) = @_;

        if ($string =~ /OH \+ PAR = 0\.11 ALD2/) {
            $string = "OH + PAR = 0.11 ALD2 + 0.11 HO2\n+ 0.11 PAROP + 0.76 ROR\n+ 0.87 XO2 + 0.13 XO2N";
        } elsif ($string =~ /OH \+ PAR = 0\.060 ALD2/) {
            $string = "OH + PAR = 0.06 ALD2 + 0.05 ALDX\n+ 0.11 HO2 + 0.11 PAROP + 0.76 ROR\n+ 0.87 XO2 + 0.13 XO2N";
        }
        return $string;
    }

    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    
    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    
    my @final_sorted_data;
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
