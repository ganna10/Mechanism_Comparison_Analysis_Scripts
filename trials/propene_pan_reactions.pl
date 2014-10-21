#!/usr/bin/perl
# PAN family budget reactions due to propene as stacked bar plots
# Version 0 : Jane Coates 31/7/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);
my @parents = ( 'C3H6', 'OLT');

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_panfile = "/work/users/jco/Mechanisms/MCM/MCM-3.2/All_pans.txt";
my @mcm_3_2_pans = get_pans($mcm_3_2_panfile); 
my @mcm_3_2_pan_producers = pan_producers($mcm_3_2_kpp, @mcm_3_2_pans);
my $ntime = $mcm_3_2_mecca->time->nelem; #number of time points 
$families{'PAN_mcm_3_2'} = [ @mcm_3_2_pans , @mcm_3_2_pan_producers ]; 
($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'}) = get_rates($parents[0], 'PAN_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_mecca = MECCA->new($radm2_run); 
my $radm2_eqnfile = "/work/users/jco/MECCA/RADM2_tagged/gas.eqn";
my $radm2_kpp = KPP->new($radm2_eqnfile); 
my $radm2_panfile = "/work/users/jco/Mechanisms/RADM2/All_pans.txt";
my @radm2_pans = get_pans($radm2_panfile); 
my @radm2_pan_producers = pan_producers($radm2_kpp, @radm2_pans);
$families{'PAN_radm2'} = [ @radm2_pans , @radm2_pan_producers ]; 
($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'}) = get_rates($parents[1], 'PAN_radm2', $radm2_kpp, $radm2_mecca);
my ($radm2_sorted_plot_data, $radm2_legend) = sort_data_for_plot($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA->new($racm_run); 
my $racm_eqnfile = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqnfile); 
my $racm_panfile = "/work/users/jco/Mechanisms/RACM/All_pans.txt";
my @racm_pans = get_pans($racm_panfile); 
my @racm_pan_producers = pan_producers($racm_kpp, @racm_pans);
$families{'PAN_racm'} = [ @racm_pans , @racm_pan_producers ]; 
($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'}) = get_rates($parents[1], 'PAN_racm', $racm_kpp, $racm_mecca);
my ($racm_sorted_plot_data, $racm_legend) = sort_data_for_plot($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'});
my $racm_plot_title = "(e) RACM";

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA->new($racm2_run); 
my $racm2_eqnfile = "/work/users/jco/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqnfile); 
my $racm2_panfile = "/work/users/jco/Mechanisms/RACM2/All_pans.txt";
my @racm2_pans = get_pans($racm2_panfile); 
my @racm2_pan_producers = pan_producers($racm2_kpp, @racm2_pans);
$families{'PAN_racm2'} = [ @racm2_pans , @racm2_pan_producers ]; 
($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'}) = get_rates($parents[1], 'PAN_racm2', $racm2_kpp, $racm2_mecca);
my ($racm2_sorted_plot_data, $racm2_legend) = sort_data_for_plot($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'});
my $racm2_plot_title = "(f) RACM2";

#Create x-axis for plot in hours
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog; 
#map to day and night
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

my ($budget_plot) = budget_plot({
        times           => \@time_blocks,
        mcm3_2_data     => $mcm_3_2_sorted_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title,
        mcm3_2_legend   => $mcm_3_2_legend,
        radm2_data      => $radm2_sorted_plot_data,
        radm2_title     => $radm2_plot_title,
        radm2_legend    => $radm2_legend,
        racm_data       => $racm_sorted_plot_data,
        racm_title      => $racm_plot_title,
        racm_legend     => $racm_legend,
        racm2_data      => $racm2_sorted_plot_data,
        racm2_title     => $racm2_plot_title,
        racm2_legend    => $racm2_legend,
});

sub get_pans {
    my ($pan_file) = @_;

    open my $in, '<:encoding(utf-8)', $pan_file or die "Can't open $pan_file : $!";
    my @lines = <$in>;
    close $in;
    my %hash;
    my @separate = map { split /\s/, $_ } @lines;
    $hash{$_} += 1 foreach (@separate);
    return keys %hash;
}

sub pan_producers {
    my ($kpp, @pans) = @_;
    
    my @pan_producers;
    foreach my $pan (@pans) {
        my $producers = $kpp->producing($pan);
        next if (@$producers == 0);
        foreach my $r_nr (@$producers) {
            my $reaction_string = $kpp->reaction_string($r_nr);
            next unless ($reaction_string =~ /\bNO2\b/);
            my ($reactants, $products) = split / = /, $reaction_string;
            $reactants =~ s/NO2//g;
            $reactants =~ s/\s\+\s//g;
            $reactants =~ s/^\s+|\s+$//g;
            push @pan_producers, $reactants;
        }
    }
    return @pan_producers;
}

sub get_rates {
    my ($VOC, $species, $kpp, $mecca) = @_;
    print "$species\n";

    my ($consumers, $producers, $consumer_yields, $producer_yields, %species_production_rates, %species_consumption_rates);
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
    }

    #check that species reactions are found
    die "No producers found for $species\n" if (@$producers == 0);
    die "No consumers found for $species\n" if (@$consumers == 0);

    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq $VOC);
        my $string = $kpp->reaction_string($reaction);
        my ($reactants, $products) = split ' = ', $string;
        $reactants =~ s/^\s+|\s+$//g;
        $reactants =~ s/_${parent}//g;
        $species_production_rates{$reactants} += $rate(1:$ntime-2); 
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq $VOC);
        my $string = $kpp->reaction_string($reaction);
        my ($reactants, $products) = split ' = ', $string;
        $reactants =~ s/^\s+|\s+$//g;
        $reactants =~ s/_${parent}//g;
        $species_consumption_rates{$reactants} += $rate(1:$ntime-2);
    }

    #get parent species emissions for each mechanism
    my $dt = $mecca->dt->at(0); #model time step
    my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
    my $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $species_production_rates{$_} /= $parent_emissions foreach (sort keys %species_production_rates);
    $species_consumption_rates{$_} /= $parent_emissions foreach (sort keys %species_consumption_rates);
    return (\%species_production_rates, \%species_consumption_rates);
}

sub sort_data_for_plot { #create hash with production of the reactions
    my ($production_rates, $consumption_rates) = @_;
    my %production_rates = %$production_rates;
    my %consumption_rates = %$consumption_rates;
    my (@production_others, @consumption_others, @sorted_plot_data); 
    my $prod_others_max = 1e-5;
    my $cons_others_max = -$prod_others_max;

    foreach my $item (keys %consumption_rates) {#sort consumption
        if ($consumption_rates{$item}->sum > $cons_others_max) { #get consumption others
            $consumption_rates{'Consumption Others'} += $consumption_rates{$item};
            delete $consumption_rates{$item};
        }
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_cons = reverse sort { &$sort_function($consumption_rates{$b}) <=> &$sort_function($consumption_rates{$a}) } keys %consumption_rates;

    foreach (@sorted_cons) { #sum up rates of reactions, starting with reaction with lowest sum, consumption others added separately 
        next if ($_ eq 'Consumption Others');
        push @sorted_plot_data, { $_ => $consumption_rates{$_} };
    }

    push @sorted_plot_data, { 'Consumption Others' => $consumption_rates{'Consumption Others'} } if (defined $consumption_rates{'Consumption Others'}); #add Consumption Others to the beginning 
    
    foreach my $item (keys %production_rates) {#sort production
        if ($production_rates{$item}->sum < $prod_others_max) { #get production others
            $production_rates{'Production Others'} += $production_rates{$item};
            delete $production_rates{$item};
        }
    }

    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;

    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); #add Production Others to the beginning 

    my (@plot_data, @legend_pos, @legend_neg, @legend);
    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            if ($ref->{$item}->sum > 0) {
                push @legend_pos, $item;
            } else {
                push @legend_neg, $item;
            }
            my $rate_list = join ":", $ref->{$item}->dog;
            my @rate_array = split /:/, $rate_list;
            push @plot_data, { $item => \@rate_array };
        }
    } 

    push @legend, reverse @legend_pos;
    push @legend, @legend_neg;

    return (\@plot_data, \@legend);
}

sub budget_plot { #create dataframe and then create plot
    my ($args) = @_;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(Cairo) `);
    $R->run(q` library(gtable) `);

    $R->set('time', [@{$args->{times}}]);

    #general plot R function
    $R->run(q` mech.plot = function(data, plot.title, plot.legend) {data = ddply(data, .(time), colwise(sum)) ;
                                                                    subset = data[1:7,] ;
                                                                    plot.data = melt( data = subset, id = names(subset)[1], measured = names(subset)[-1] ); 
                                                                    colnames(plot.data) = c("time", "Reaction", "rate"); 
                                                                    Reaction.levels = (levels(factor(plot.data$Reaction))); 
                                                                    plot.data$Reaction = ordered(plot.data$Reaction, levels = Reaction.levels); 
                                                                    plot.data = ddply( plot.data, .(Reaction)); 
                                                                    plot = ggplot(data = plot.data, aes(x = time, y = rate)); 
                                                                    plot = plot + geom_bar(data = subset(plot.data, rate > 0), aes(fill = Reaction), stat = "identity", width = 0.5) ;
                                                                    plot = plot + geom_bar(data = subset(plot.data, rate < 0), aes(fill = Reaction), stat = "identity", width = 0.5) ;
                                                                    plot = plot + ggtitle(plot.title); 
                                                                    plot = plot + theme_bw() ; 
                                                                    plot = plot + theme(legend.key.size = unit(10.0, "cm")) ; 
                                                                    plot = plot + theme(axis.text.x = element_text(size = 140)) ; 
                                                                    plot = plot + scale_fill_discrete(breaks = plot.legend);
                                                                    plot = plot + theme(axis.text.y = element_text(size = 150)) ; 
                                                                    plot = plot + theme(axis.title.x = element_blank()) ; 
                                                                    plot = plot + scale_y_continuous(limits = c(-50, 80), breaks = seq(-50, 80, 20)) ;
                                                                    plot = plot + theme(legend.text = element_text(size = 200)) ; 
                                                                    plot = plot + theme(legend.key = element_blank()) ; 
                                                                    plot = plot + theme(axis.title.y = element_blank()) ; 
                                                                    plot = plot + theme(plot.title = element_text(size = 200, face = "bold", vjust = 1)) ; 
                                                                    plot = plot + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) ;
                                                                    plot = plot + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) ;
                                                                    plot = plot + theme(legend.title = element_text(size = 250, face = "bold")); 
                                                                    return(plot) } `);
 
     #MCM v3.2
     $R->run(q` mcm3.2.data = data.frame(time)`);
     $R->set('mcm3.2.plot.title', $args->{mcm3_2_title});
     $R->set('mcm3.2.legend', $args->{mcm3_2_legend});
     foreach my $ref (@{$args->{mcm3_2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` mcm3.2.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, mcm3.2.plot.title, mcm3.2.legend) `);
 
     #RADM2
     $R->run(q` radm2.data = data.frame(time)`);
     $R->set('radm2.plot.title', $args->{radm2_title});
     $R->set('radm2.legend', $args->{radm2_legend});
     foreach my $ref (@{$args->{radm2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` radm2.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` radm2.plot = mech.plot(radm2.data, radm2.plot.title, radm2.legend) `);
 
     #RACM
     $R->run(q` racm.data = data.frame(time)`);
     $R->set('racm.plot.title', $args->{racm_title});
     $R->set('racm.legend', $args->{racm_legend});
     foreach my $ref (@{$args->{racm_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` racm.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` racm.plot = mech.plot(racm.data, racm.plot.title, racm.legend) `);
 
     #RACM2
     $R->run(q` racm2.data = data.frame(time)`);
     $R->set('racm2.plot.title', $args->{racm2_title});
     $R->set('racm2.legend', $args->{racm2_legend});
     foreach my $ref (@{$args->{racm2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` racm2.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` racm2.plot = mech.plot(racm2.data, racm2.plot.title, racm2.legend) `);
 
     $R->run(q` CairoPDF(file = "propene_pan_budget_reactions.pdf", width = 210, height = 200) `, 
            q` y.label = textGrob(expression(bold(paste("\nMolecules (Product) ", s^-1, "/Molecules (NMVOC) x ", 10^5))), rot = 90, gp = gpar(fontsize = 210), vjust = 0.6)`,
            q` main.plot = grid.arrange(y.label, 
                                        arrangeGrob(mcm3.2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    radm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    racm.plot ,
                                                    racm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    nrow = 2), 
                                       nrow = 1, ncol = 2,
                                       sub = textGrob("\n", gp = gpar(fontsize = 210, fontface = "bold"), vjust = 0.2), 
                                       widths=unit.c(unit(20, "lines"), unit(1, "npc") - unit(20, "lines"))) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );

    $R->stop(); 
}
