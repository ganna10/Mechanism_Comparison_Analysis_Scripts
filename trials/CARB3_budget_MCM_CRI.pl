#!/usr/bin/perl -w
# overlay plot of glyoxal production due to toluene degradation
# Version 0: Jane Coates 13/08/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my (%production_rates, %consumption_rates);
my $prod_others_max = 5.0e-5;

my $cri_species = "CARB3";
my $MCM_species = "GLYOX";

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $ntime = $mcm_3_2_mecca->time->nelem; #number of time points 
($production_rates{$MCM_species}, $consumption_rates{$MCM_species}) = get_rates($MCM_species, $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($prod_others_max, $production_rates{$MCM_species}, $consumption_rates{$MCM_species});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run); 
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile); 
($production_rates{$cri_species}, $consumption_rates{$cri_species}) = get_rates($cri_species, $cri_kpp, $cri_mecca);
my ($cri_sorted_plot_data, $cri_legend) = sort_data_for_plot($prod_others_max, $production_rates{$cri_species}, $consumption_rates{$cri_species});
my $cri_plot_title = "(b) CRI v2";

#Create x-axis for plot in hours -> time axis is the same for all model runs
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 3600;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list; 
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

my ($plot) = plot({#create plot
        times           => \@time_blocks,
        mcm3_2_data     => $mcm_3_2_sorted_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title,
        mcm3_2_legend   => $mcm_3_2_legend,
        cri_data        => $cri_sorted_plot_data,
        cri_title       => $cri_plot_title,
        cri_legend      => $cri_legend,
});

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

sub get_rates { #get production rates normalised by emissions of the reactions of the intermediate species of a VOC
    my ($species, $kpp, $mecca) = @_;

    my $look_up_species = $species . "_TOLUENE";
    my ($producers, $producer_yields, %species_production_rates, $consumers, $consumer_yields, %species_consumption_rates);
    $producers = $kpp->producing($look_up_species);
    $producer_yields = $kpp->effect_on($look_up_species, $producers);  
    $consumers = $kpp->consuming($look_up_species);
    $consumer_yields = $kpp->effect_on($look_up_species, $consumers);  
    die "No producers found for $look_up_species\n" if (@$producers == 0);
    die "No consumers found for $look_up_species\n" if (@$consumers == 0);
    
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq "TOLUENE");
		my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $string = $kpp->reaction_string($reaction);
        my ($reactants) = $string =~ /^(.*?)\s=/;
        $reactants =~ s/_$parent//g; #remove tag from string
        $species_production_rates{$reactants} += $rate(1:$ntime-2); 
    } 

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq "TOLUENE");
		my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $string = $kpp->reaction_string($reaction);
        my ($reactants) = $string =~ /^(.*?)\s=/;
        $reactants =~ s/_$parent//g; #remove tag from string
        $species_consumption_rates{$reactants} += $rate(1:$ntime-2); 
    } 

    #get parent species emissions for each mechanism
    my $dt = $mecca->dt->at(0); #model time step
    my $parent_source = $mecca->balance('TOLUENE'); #in molecules (VOC)/cm3/s
    my $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $species_production_rates{$_} /= $parent_emissions foreach (sort keys %species_production_rates);
    $species_consumption_rates{$_} /= $parent_emissions foreach (sort keys %species_consumption_rates);
    return (\%species_production_rates, \%species_consumption_rates);
}

sub sort_data_for_plot { #create hash with production of the reactions
    my ($prod_others_max, $production_rates, $consumption_rates) = @_;
    my %production_rates = %$production_rates;
    my %consumption_rates = %$consumption_rates;
    my (@production_others, @consumption_others, @sorted_plot_data); 
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

sub plot { #create dataframe and then create plot
    my ($args) = @_;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(Cairo) `);

    $R->set('time', [@{$args->{times}}]);
    
    #general plot R function
    $R->run(q` mech.plot = function(data, plot.title) {   data = ddply(data, .(time), colwise(sum)) ;
                                                                  plot.data = melt(data = data, id = names(data)[1], measured = names(data)[-1] ); 
                                                                  colnames(plot.data) = c("time", "reaction", "rate"); 
                                                                  reaction.levels = (levels(factor(plot.data$reaction))); 
                                                                  plot.data$reaction = ordered(plot.data$reaction, levels = reaction.levels); 
                                                                  plot = ggplot(data = plot.data, aes(x = time, y = rate, fill = reaction)); 
                                                                  plot = plot + geom_bar(data = subset(plot.data, rate > 0), stat = "identity", width = 0.6) ;
                                                                  plot = plot + geom_bar(data = subset(plot.data, rate < 0), stat = "identity", width = 0.6) ;
                                                                  plot = plot + ggtitle(plot.title); 
                                                                  plot = plot + scale_y_continuous(limits = c(-2.5, 3), breaks = seq(-2.5, 3, 0.5)) ;
                                                                  plot = plot + theme_bw() ; 
                                                                  plot = plot + theme(legend.key.size = unit(6, "cm")) ; 
                                                                  plot = plot + theme(axis.text.x = element_text(size = 90)) ; 
                                                                  plot = plot + theme(axis.text.y = element_text(size = 100)) ; 
                                                                  plot = plot + theme(legend.text = element_text(size = 110)) ; 
                                                                  plot = plot + theme(legend.key = element_blank()) ; 
                                                                  plot = plot + theme(axis.title.y = element_blank()) ; 
                                                                  plot = plot + theme(axis.title.x = element_blank()) ; 
                                                                  plot = plot + theme(plot.title = element_text(size = 130, face = "bold", vjust = 0)) ; 
                                                                  plot = plot + theme(legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95)) ; 
                                                                  plot = plot + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) ;
                                                                  plot = plot + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) ;
                                                                  return(plot) } `);
 
     #MCM v3.2
     $R->run(q` mcm3.2.data = data.frame(time)`);
     $R->set('mcm3.2.file.name', $args->{mcm3_2_file}); 
     $R->set('mcm3.2.plot.title', $args->{mcm3_2_title});
     $R->set('mcm3.2.legend', $args->{mcm3_2_legend});
     foreach my $ref (@{$args->{mcm3_2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` mcm3.2.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, mcm3.2.plot.title) `);
 
     #cri
     $R->run(q` cri.data = data.frame(time)`);
     $R->set('cri.file.name', $args->{cri_file}); 
     $R->set('cri.plot.title', $args->{cri_title});
     $R->set('cri.legend', $args->{cri_legend});
     foreach my $ref (@{$args->{cri_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` cri.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` cri.plot = mech.plot(cri.data, cri.plot.title) `); 
 
     #colours
     $R->run(q` mcm.colours = c("Production Others" = "#696537", "GLYOX + OH" = "#8d1435", "GLYOX + hv" = "#76afca", "HCOCOHPAN + OH" = "#dc3522", "OH + TLOBIPEROH" = "#8c6238", "TLBIPERO" = "#86b650", "Consumption Others" = "#ee6738") `,
             q` cri.colours = c("Production Others" = "#696537", "CARB3 + OH" = "#8d1435", "CARB3 + hv" = "#76afca", "NO + RA16O2" = "#0e5c28", "OH + RA16NO3" = "#e7e85e", "Consumption Others" = "#ee6738") `
     );

     #multiplot
     $R->run(q` CairoPDF(file = "Glyoxal_MCM_CRI_Ox_intermediates.pdf", width = 130, height = 70) `,
             q` y.label = textGrob(expression(bold(paste("Molecules (intermediate) ", s^-1, "/Molecules (VOC) x ", 10^4))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5)`,
             q` main.plot = grid.arrange(y.label, 
                                         arrangeGrob(mcm3.2.plot + scale_fill_manual( name = "Reaction", values = mcm.colours, limits = mcm3.2.legend),
                                                     cri.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(name = "Reaction", limits = cri.legend, values = cri.colours) ,
                                                     nrow = 1), 
                                         nrow = 1, 
                                         sub = textGrob("\n", gp = gpar(fontsize = 140)),
                                         widths = unit.c(unit(14, "lines"), unit(1, "npc") - unit(14, "lines"))) `,
             q` print(main.plot) `,
             q` dev.off() `,
     );
     
     $R->stop(); 
}
