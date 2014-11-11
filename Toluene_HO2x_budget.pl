#!/usr/bin/perl -w
# overlay plot of TOLUENE organic VOC intermediates reaction rates normalised by parent VOC total emissions, showing which intermediates are responsible for HO2x production, for CRI and MCM 3.2 - stacked bar plot
# Version 0: Jane Coates 24/06/2014
# Version 1: Jane Coates 22/09/2014 removed y-axis ticks in CRI plot, changed fill order
# Version 2: Jane Coates 26/09/2014 including RACM2

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);
my $prod_others_max = 1.5e-4;

my $species = "TOLUENE";
my $racm2_species = "TOL";

#MCMv3.2 data
my $mcm_3_2_run = "/local/home/coates/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/local/home/coates/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
$families{'HO2x_mcm_3_2'} = [ qw(HO2 HO2NO2 ) ];
$weights{'HO2x_mcm_3_2'} = { }; 
my $ntime = $mcm_3_2_mecca->time->nelem; #number of time points 
($production_rates{'HO2x_mcm_3_2'}, $consumption_rates{'HO2x_mcm_3_2'}) = get_rates($species, 'HO2x_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($prod_others_max, $production_rates{'HO2x_mcm_3_2'}, $consumption_rates{'HO2x_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#cri data
my $cri_run = "/local/home/coates/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run); 
my $cri_eqnfile = "/local/home/coates/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile); 
$families{'HO2x_cri'} = [ qw( HO2 HO2NO2 ) ];
$weights{'HO2x_cri'} = { }; 
($production_rates{'HO2x_cri'}, $consumption_rates{'HO2x_cri'}) = get_rates($species, 'HO2x_cri', $cri_kpp, $cri_mecca);
my ($cri_sorted_plot_data, $cri_legend) = sort_data_for_plot($prod_others_max, $production_rates{'HO2x_cri'}, $consumption_rates{'HO2x_cri'});
my $cri_plot_title = "(b) CRI v2";

#racm2 data
my $racm2_run = "/local/home/coates/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA->new($racm2_run); 
my $racm2_eqnfile = "/local/home/coates/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqnfile); 
$families{'HO2x_racm2'} = [ qw( HO2 HO2NO2 ) ];
$weights{'HO2x_racm2'} = { }; 
($production_rates{'HO2x_racm2'}, $consumption_rates{'HO2x_racm2'}) = get_rates($racm2_species, 'HO2x_racm2', $racm2_kpp, $racm2_mecca);
my ($racm2_sorted_plot_data, $racm2_legend) = sort_data_for_plot($prod_others_max, $production_rates{'HO2x_racm2'}, $consumption_rates{'HO2x_racm2'});
my $racm2_plot_title = "(c) RACM2";

#Create x-axis for plot in hours -> time axis is the same for all model runs
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

my ($plot) = plot({#create plot
        times           => \@time_blocks,
        mcm3_2_data     => $mcm_3_2_sorted_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title,
        mcm3_2_legend   => $mcm_3_2_legend,
        cri_data        => $cri_sorted_plot_data,
        cri_title       => $cri_plot_title,
        cri_legend      => $cri_legend,
        racm2_data      => $racm2_sorted_plot_data,
        racm2_title     => $racm2_plot_title,
        racm2_legend    => $racm2_legend,
});

sub get_rates { #get production rates normalised by emissions of the reactions of the intermediate species of a VOC
    my ($parent_species, $species, $kpp, $mecca) = @_;

    my $look_up_species = $species . "_" . $parent_species;
    my ($producers, $producer_yields, %species_production_rates, $consumers, $consumer_yields, %species_consumption_rates);
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
        die "No producers found for $species\n" if (@$producers == 0);
        die "No consumers found for $species\n" if (@$consumers == 0);
    } else {
        print "No family found!\n";
    }
    
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq $parent_species);
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
        next unless (defined $parent and $parent eq $parent_species);
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
    my $parent_source = $mecca->balance($parent_species); #in molecules (VOC)/cm3/s
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
            my @rate_array = map { $_ } $ref->{$item}->dog;
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
    $R->run(q` mech.plot = function(data, plot.title) {plot = ggplot(data = data, aes(x = time, y = rate, fill = Reaction)) ; 
                                                       plot = plot + geom_bar(data = subset(data, rate > 0), stat = "identity", width = 0.6) ;
                                                      plot = plot + geom_bar(data = subset(data, rate < 0), stat = "identity", width = 0.6) ;
                                                      plot = plot + ggtitle(plot.title); 
                                                      plot = plot + scale_y_continuous(limits = c(-1, 18), breaks = seq(-1, 18, 2)) ;
                                                      plot = plot + theme_bw() ; 
                                                      plot = plot + theme(legend.key.size = unit(4, "cm")) ; 
                                                      plot = plot + theme(axis.text.x = element_text(size = 70, angle = 45, vjust = 0.5)) ; 
                                                      plot = plot + theme(axis.text.y = element_text(size = 60)) ; 
                                                      plot = plot + theme(legend.text = element_text(size = 40)) ; 
                                                      plot = plot + theme(legend.title = element_blank()) ; 
                                                      plot = plot + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) ;
                                                      plot = plot + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) ;
                                                      plot = plot + theme(axis.title.y = element_blank()) ; 
                                                      plot = plot + theme(axis.title.x = element_blank()) ; 
                                                      plot = plot + theme(axis.ticks.length = unit(0.5, "cm")) ;
                                                      plot = plot + theme(axis.ticks.margin = unit(0.3, "cm")) ;
                                                      plot = plot + theme(plot.title = element_text(size = 90, face = "bold", vjust = 0)) ; 
                                                      plot = plot + theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.99)) ; 
                                                      return(plot) } `);
 
     #MCM v3.2
     $R->run(q` mcm3.2.data = data.frame(time)`);
     $R->set('mcm3.2.file.name', $args->{mcm3_2_file}); 
     $R->set('mcm3.2.plot.title', $args->{mcm3_2_title});
     $R->set('mcm3.2.legend', $args->{mcm3_2_legend});
     foreach my $ref (@{$args->{mcm3_2_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` mcm3.2.data[name] = rate*10000`); 
         }
     } 

     $R->run(q` mcm3.2.data = ddply(mcm3.2.data, .(time), colwise(sum))  `,
            q` mcm3.2.data = mcm3.2.data[1:7,] `,
            q` mcm3.2.data = melt(data = mcm3.2.data, id = names(mcm3.2.data)[1], measured = names(mcm3.2.data)[-1] ) `, 
            q` colnames(mcm3.2.data) = c("time", "Reaction", "rate") `, 
            q` mcm3.2.data$Reaction = factor(mcm3.2.data$Reaction, levels = c("Consumption Others", "CO + OH", "HCHO + OH", "HCHO + hv", "GLYOX + OH", "CH3O", "TLBIPERO", "OH + TOLUENE", "MGLYOX + hv", "Production Others")) `,
            q` mcm3.2.data = ddply( mcm3.2.data, .(Reaction)) `, 
    );
    $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, mcm3.2.plot.title) `);
 
     #CRI
     $R->run(q` cri.data = data.frame(time)`);
     $R->set('cri.file.name', $args->{cri_file}); 
     $R->set('cri.plot.title', $args->{cri_title});
     $R->set('cri.legend', $args->{cri_legend});
     foreach my $ref (@{$args->{cri_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` cri.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` cri.data = ddply(cri.data, .(time), colwise(sum))  `,
            q` cri.data = cri.data[1:7,] `,
            q` cri.data = melt(data = cri.data, id = names(cri.data)[1], measured = names(cri.data)[-1] ) `, 
            q` colnames(cri.data) = c("time", "Reaction", "rate") `, 
            q` cri.data$Reaction = factor(cri.data$Reaction, levels = c("Consumption Others", "CO + OH", "HCHO + OH", "HCHO + hv", "CARB3 + OH", "CH3O2 + NO", "C2H5O2 + NO", "HOCH2CH2O2 + NO", "NO + RN10O2", "CARB3 + hv", "Production Others")) `, 
            q` cri.data = ddply( cri.data, .(Reaction)) `, 
    );
     $R->run(q` cri.plot = mech.plot(cri.data, cri.plot.title) `); 
 
     #racm2
     $R->run(q` racm2.data = data.frame(time)`);
     $R->set('racm2.file.name', $args->{racm2_file}); 
     $R->set('racm2.plot.title', $args->{racm2_title});
     $R->set('racm2.legend', $args->{racm2_legend});
     foreach my $ref (@{$args->{racm2_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` racm2.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` racm2.data = ddply(racm2.data, .(time), colwise(sum))  `,
            q` racm2.data = racm2.data[1:7,] `,
            q` racm2.data = melt(data = racm2.data, id = names(racm2.data)[1], measured = names(racm2.data)[-1] ) `, 
            q` colnames(racm2.data) = c("time", "Reaction", "rate") `, 
            q` racm2.data$Reaction = factor(racm2.data$Reaction, levels = c("Consumption Others", "CO + OH", "HCHO + OH", "HCHO + hv", "MO2 + NO", "ETHP + NO", "EPX + O3", "TR2", "Production Others")) `, 
            q` racm2.data = ddply( racm2.data, .(Reaction)) `, 
    );
     $R->run(q` racm2.plot = mech.plot(racm2.data, racm2.plot.title) `); 
 
     #colours
     $R->run(q` mcm.colours = c("Production Others" = "#696537", "CO + OH" = "#f9c600", "HCHO + hv" = "#76afca", "HCHO + OH" = "#dc3522", "GLYOX + NO" = "#8c6238", "MGLYOX + hv" = "#9bb08f", "CH3CO3 + NO" = "#8b1537", "CH3O2 + NO" = "#e7e85e", "OH + TOLUENE" = "#0352cb", "TLBIPERO" = "#6c254f", "CH3O" = "#86b650", "GLYOX + OH" = "#8d1435", "Consumption Others" = "#ee6738") `,
             q` cri.colours = c("Production Others" = "#696537", "CO + OH" = "#f9c600", "HCHO + hv" = "#76afca", "HCHO + OH" = "#dc3522", "CH3CO3 + NO" = "#8b1537", "HOCH2CH2O2 + NO" = "#8c6238", "CH3O2 + NO" = "#e7e85e", "C2H5O2 + NO" = "#9bb08f", "NO + RA16O2" = "#0e5c28", "CARB3 + hv" = "#76afca", "CARB3 + OH" = "#8d1435", "NO + RN10O2" = "#b569b3", "Consumption Others" = "#ee6738") `,
             q` racm2.colours = c("Production Others" = "#696537", "CO + OH" = "#f9c600", "HCHO + hv" = "#76afca", "HCHO + OH" = "#dc3522", "EPX + O3" = "#0c3f78", "MO2 + NO" = "#e7e85e", "ETHP + NO" = "#9bb08f", "TR2" = "#cc6329", "Consumption Others" = "#ee6738") `,
     );

     #multiplot
     $R->run(q` CairoPDF(file = "TOL_HO2x_intermediates.pdf", width = 57, height = 40) `,
             q` y.label = textGrob(expression(bold(paste("Molecules (intermediate) ", s^-1, "/Molecules (VOC) x ", 10^4))), rot = 90, gp = gpar(fontsize = 85), vjust = 0.5)`,
             q` main.plot = grid.arrange(y.label, 
                                         arrangeGrob(mcm3.2.plot + scale_fill_manual( name = "Reaction", values = mcm.colours, limits = rev(levels(mcm3.2.data$Reaction))) , 
                                                     cri.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual( name = "Reaction", limits = rev(levels(cri.data$Reaction)), values = cri.colours) , 
                                                     racm2.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(values = racm2.colours, limits = rev(levels(racm2.data$Reaction))),
                                                     nrow = 1), 
                                         nrow = 1, 
                                         sub = textGrob("\n", gp = gpar(fontsize = 140)),
                                         widths = unit.c(unit(14, "lines"), unit(1, "npc") - unit(14, "lines"))) `,
             q` print(main.plot) `,
             q` dev.off() `,
     );
     
     $R->stop(); 
}
