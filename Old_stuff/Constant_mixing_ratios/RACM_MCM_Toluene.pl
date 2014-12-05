#!/usr/bin/perl -w
# overlay plot of TOLUENE organic VOC intermediates reaction rates normalised by parent VOC total emissions, showing which intermediates are responsible for Ox production, for RACM and MCM 3.2 - stacked bar plot
# Version 0: Jane Coates 06/05/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);
my $prod_others_max = 1.1e-4;

#prepare HO2 family data
$families{'HO2x_racm'} = [ qw( HO2 HO2NO2 ) ];
$families{'HO2x_mcm'} = [ qw( HO2 HO2NO2 ) ];

my $RACM_species = "TOL";
my $MCM_species = "TOLUENE";

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_ro2file = "/work/users/jco/MECCA/MCM_3.2_tagged/RO2_species.txt";
my @mcm_3_2_no2_reservoirs = get_no2_reservoirs($mcm_3_2_kpp, $mcm_3_2_ro2file);
$families{'Ox_mcm_3_2'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @mcm_3_2_no2_reservoirs ];
$weights{'Ox_mcm_3_2'} = { NO3 => 2, N2O5 => 3};

my $ntime = $mcm_3_2_mecca->time->nelem; #number of time points

($production_rates{'Ox_mcm_3_2'}, $consumption_rates{'Ox_mcm_3_2'}) = get_rates($MCM_species, 'Ox_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($prod_others_max, $production_rates{'Ox_mcm_3_2'}, $consumption_rates{'Ox_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) Toluene - Organic Ox Budget using MCM v3.2";

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA->new($racm_run); 
my $racm_eqnfile = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqnfile); 
my $racm_ro2file = "/work/users/jco/MECCA/RACM_tagging/RO2_species.txt";
my @racm_no2_reservoirs = get_no2_reservoirs($racm_kpp, $racm_ro2file);
$families{'Ox_racm'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @racm_no2_reservoirs ];
$weights{'Ox_racm'} = { NO3 => 2, N2O5 => 3};

($production_rates{'Ox_racm'}, $consumption_rates{'Ox_racm'}) = get_rates($RACM_species, 'Ox_racm', $racm_kpp, $racm_mecca);
my ($racm_sorted_plot_data, $racm_legend) = sort_data_for_plot($prod_others_max, $production_rates{'Ox_racm'}, $consumption_rates{'Ox_racm'});
my $racm_plot_title = "(b) Toluene - Organic Ox Budget using RACM";

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
        racm_data       => $racm_sorted_plot_data,
        racm_title      => $racm_plot_title,
        racm_legend     => $racm_legend,
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
    $R->run(q` mech.plot = function(data, legend.title, plot.title, legend) { data = ddply(data, .(time), colwise(sum)) ;
                                                                              plot.data = melt(data = data, id = names(data)[1], measured = names(data)[-1] ); 
                                                                              colnames(plot.data) = c("time", "reaction", "rate"); 
                                                                              reaction.levels = (levels(factor(plot.data$reaction))); 
                                                                              plot.data$reaction = ordered(plot.data$reaction, levels = reaction.levels); 
                                                                              plot.data = ddply( plot.data, .(reaction)); 
                                                                              cols = colorRampPalette(brewer.pal(9, "Accent"))(nlevels(plot.data$reaction)); 
                                                                              plot = ggplot(data = plot.data, aes(x = time, y = rate, fill = reaction)); 
                                                                              plot = plot + guides(fill = guide_legend(title = legend.title)); 
                                                                              plot = plot + geom_bar(data = subset(plot.data, rate > 0), stat = "identity", width = 0.6) ;
                                                                              plot = plot + geom_bar(data = subset(plot.data, rate < 0), stat = "identity", width = 0.6) ;
                                                                              plot = plot + ggtitle(paste("\n", plot.title, "\n")); 
                                                                              plot = plot + scale_y_continuous(limits = c(-25, 15), breaks = seq(-25, 15, 5)) ;
                                                                              plot = plot + theme_bw() ; 
                                                                              plot = plot + theme(legend.key.size = unit(6, "cm")) ; 
                                                                              plot = plot + theme(axis.text.x = element_text(size = 90)) ; 
                                                                              plot = plot + theme(axis.text.y = element_text(size = 100)) ; 
                                                                              plot = plot + theme(legend.text = element_text(size = 110)) ; 
                                                                              plot = plot + theme(legend.title = element_text(size = 120, face = "bold")) ; 
                                                                              plot = plot + theme(legend.key = element_blank()) ; 
                                                                              plot = plot + theme(axis.title.y = element_blank()) ; 
                                                                              plot = plot + theme(plot.title = element_text(size = 130, face = "bold", vjust = 0)) ; 
                                                                              plot = plot + theme(legend.justification = c(0.95, 0.05), legend.position = c(0.95, 0.05)) ; 
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
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` mcm3.2.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, "MCM v3.2", mcm3.2.plot.title, mcm3.2.legend) `);
 
     #RACM
     $R->run(q` racm.data = data.frame(time)`);
     $R->set('racm.file.name', $args->{racm_file}); 
     $R->set('racm.plot.title', $args->{racm_title});
     $R->set('racm.legend', $args->{racm_legend});
     foreach my $ref (@{$args->{racm_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` racm.data[name] = rate*10000`); 
         }
     } 
     $R->run(q` racm.plot = mech.plot(racm.data, "RACM", racm.plot.title, racm.legend) `); 
 
     #colours
     $R->run(q` mcm.colours = c("Production Others" = "#696537", "MECOACETO2 + NO" = "#f9c600", "NO + TLFUO2" = "#76afca", "C6H5O2 + NO" = "#dc3522", "MCATEC1O2 + NO" = "#8c6238", "NO + TLBIPERO2" = "#9bb08f", "CH3CO3 + NO" = "#8b1537", "CH3O2 + NO" = "#e7e85e", "MCATEC1O + O3" = "#0352cb", "C6H5O + O3" = "#86b650", "NPHEN1O + O3" = "#6c254f", "Consumption Others" = "#ee6738") `,
             q` racm.colours = c("Production Others" = "#696537", "ETHP + NO" = "#f9c600", "KETP + NO" = "#76afca", "NO + XO2" = "#dc3522", "ACO3 + NO" = "#8b1537", "HC3P + NO" = "#8c6238", "MO2 + NO" = "#e7e85e", "OH + ONIT" = "#9bb08f", "ADDC + O3" = "#0352cb", "CSL + NO3" = "#86b650", "ADDT + O3" = "#6c254f", "NO2 + PHO" = "#8ed6d5", "Consumption Others" = "#ee6738") `
     );

     #multiplot
     $R->run(q` CairoPDF(file = "TOL_MCM_RACM_Ox_intermediates.pdf", width = 130, height = 70) `,
             q` y.label = textGrob(expression(bold(paste("Molecules (intermediate) ", s^-1, "/Molecules (VOC) x ", 10^4))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5)`,
             q` main.plot = grid.arrange(y.label, 
                                         arrangeGrob(mcm3.2.plot + scale_fill_manual( name = "reaction", limits = mcm3.2.legend, values = mcm.colours) , 
                                                     racm.plot + theme(axis.text.y = element_blank()) + scale_fill_manual( name = "reaction", limits = racm.legend, values = racm.colours) , 
                                                     nrow = 1), 
                                         nrow = 1, 
                                         sub = textGrob("\n", gp = gpar(fontsize = 140)),
                                         widths = unit.c(unit(14, "lines"), unit(1, "npc") - unit(14, "lines"))) `,
             q` print(main.plot) `,
             q` dev.off() `,
     );
     
     $R->stop(); 
}
