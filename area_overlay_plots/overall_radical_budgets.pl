#!/usr/bin/perl
# Version 0 : Jane Coates 17/02/2014 radical budget plot attributed back to the parent VOCs

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);

#MCMv3.1 data
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA->new($mcm_3_1_run); 
my $mcm_3_1_eqnfile = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.eqn";
my $mcm_3_1_kpp = KPP->new($mcm_3_1_eqnfile); 
$families{'radicals_mcm_3_1'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points 
($production_rates{'radicals_mcm_3_1'}, $consumption_rates{'radicals_mcm_3_1'}) = get_rates('radicals_mcm_3_1', $mcm_3_1_kpp, $mcm_3_1_mecca);
my ($mcm_3_1_sorted_plot_data, $mcm_3_1_legend) = sort_data_for_plot($production_rates{'radicals_mcm_3_1'}, $consumption_rates{'radicals_mcm_3_1'});
my $mcm_3_1_plot_title = "(b) MCM v3.1";

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
$families{'radicals_mcm_3_2'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_mcm_3_2'}, $consumption_rates{'radicals_mcm_3_2'}) = get_rates('radicals_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($production_rates{'radicals_mcm_3_2'}, $consumption_rates{'radicals_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run); 
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile); 
$families{'radicals_cri'} = [ 'OH', 'HO2', 'HO2NO2' ];
($production_rates{'radicals_cri'}, $consumption_rates{'radicals_cri'}) = get_rates('radicals_cri', $cri_kpp, $cri_mecca); 
my ($cri_sorted_plot_data, $cri_legend) = sort_data_for_plot($production_rates{'radicals_cri'}, $consumption_rates{'radicals_cri'});
my $cri_plot_title = "(c) CRI v2";

#mozart data
my $mozart_run = "/work/users/jco/MECCA/MOZART_tagging/boxmodel";
my $mozart_mecca = MECCA->new($mozart_run); 
my $mozart_eqnfile = "/work/users/jco/MECCA/MOZART_tagging/gas.eqn";
my $mozart_kpp = KPP->new($mozart_eqnfile); 
$families{'radicals_mozart'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_mozart'}, $consumption_rates{'radicals_mozart'}) = get_rates('radicals_mozart', $mozart_kpp, $mozart_mecca);
my ($mozart_sorted_plot_data, $mozart_legend) = sort_data_for_plot($production_rates{'radicals_mozart'}, $consumption_rates{'radicals_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_mecca = MECCA->new($radm2_run); 
my $radm2_eqnfile = "/work/users/jco/MECCA/RADM2_tagged/gas.eqn";
my $radm2_kpp = KPP->new($radm2_eqnfile); 
$families{'radicals_radm2'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_radm2'}, $consumption_rates{'radicals_radm2'}) = get_rates('radicals_radm2', $radm2_kpp, $radm2_mecca);
my ($radm2_sorted_plot_data, $radm2_legend) = sort_data_for_plot($production_rates{'radicals_radm2'}, $consumption_rates{'radicals_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA->new($racm_run); 
my $racm_eqnfile = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqnfile); 
$families{'radicals_racm'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_racm'}, $consumption_rates{'radicals_racm'}) = get_rates('radicals_racm', $racm_kpp, $racm_mecca);
my ($racm_sorted_plot_data, $racm_legend) = sort_data_for_plot($production_rates{'radicals_racm'}, $consumption_rates{'radicals_racm'});
my $racm_plot_title = "(e) RACM";

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA->new($racm2_run); 
my $racm2_eqnfile = "/work/users/jco/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqnfile); 
$families{'radicals_racm2'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_racm2'}, $consumption_rates{'radicals_racm2'}) = get_rates('radicals_racm2', $racm2_kpp, $racm2_mecca);
my ($racm2_sorted_plot_data, $racm2_legend) = sort_data_for_plot($production_rates{'radicals_racm2'}, $consumption_rates{'radicals_racm2'});
my $racm2_plot_title = "(f) RACM2";

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA->new($cbm4_run);
my $cbm4_eqnfile = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqnfile);
$families{'radicals_cbm4'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_cbm4'}, $consumption_rates{'radicals_cbm4'}) = get_rates('radicals_cbm4', $cbm4_kpp, $cbm4_mecca);
my ($cbm4_sorted_plot_data, $cbm4_legend) = sort_data_for_plot($production_rates{'radicals_cbm4'}, $consumption_rates{'radicals_cbm4'});
my $cbm4_plot_title = "(h) CBM-IV";

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA->new($cb05_run);
my $cb05_eqnfile = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqnfile);
$families{'radicals_cb05'} = [ 'OH', 'HO2', 'HO2NO2' ]; 
($production_rates{'radicals_cb05'}, $consumption_rates{'radicals_cb05'}) = get_rates('radicals_cb05', $cb05_kpp, $cb05_mecca);
my ($cb05_sorted_plot_data, $cb05_legend) = sort_data_for_plot($production_rates{'radicals_cb05'}, $consumption_rates{'radicals_cb05'});
my $cb05_plot_title = "(i) CB05";

#Create x-axis for plot in hours
my $times = $mcm_3_1_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list;

my ($budget_plot) = budget_plot({
        y_max           => 2.1e7,
        y_min           => -1.7e7,
        breaks          => 5e6,
        times           => \@time_axis,
        mcm3_1_data     => $mcm_3_1_sorted_plot_data,
        mcm3_1_title    => $mcm_3_1_plot_title,
        mcm3_1_legend   => $mcm_3_1_legend,
        mcm3_2_data     => $mcm_3_2_sorted_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title,
        mcm3_2_legend   => $mcm_3_2_legend,
        cri_data        => $cri_sorted_plot_data,
        cri_title       => $cri_plot_title,
        cri_legend      => $cri_legend,
        mozart_data     => $mozart_sorted_plot_data,
        mozart_title    => $mozart_plot_title,
        mozart_legend   => $mozart_legend,
        radm2_data      => $radm2_sorted_plot_data,
        radm2_title     => $radm2_plot_title,
        radm2_legend    => $radm2_legend,
        racm_data       => $racm_sorted_plot_data,
        racm_title      => $racm_plot_title,
        racm_legend     => $racm_legend,
        racm2_data      => $racm2_sorted_plot_data,
        racm2_title     => $racm2_plot_title,
        racm2_legend    => $racm2_legend,
        cbm4_data       => $cbm4_sorted_plot_data,
        cbm4_title      => $cbm4_plot_title,
        cbm4_legend     => $cbm4_legend,
        cb05_data       => $cb05_sorted_plot_data,
        cb05_title      => $cb05_plot_title,
        cb05_legend     => $cb05_legend,
});

sub get_rates {
    my ($species, $kpp, $mecca) = @_;

    my ($consumers, $producers, $consumer_yields, $producer_yields, %species_production_rates, %production_reaction_rates, %species_consumption_rates, %consumption_reaction_rates);
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
        my $string;
        if (defined $parent) { # for tagged reactions
            $string = $kpp->reaction_string($reaction);
            $string =~ s/_$parent//g; #removing tag from reaction strings
            $species_production_rates{$species}{$parent}{$string} += $rate;
            $string = $parent; # in order to merge all production rates from all parent species reactions into 1 pdl
        } else { # for non-tagged reactions
            $string = $kpp->reaction_string($reaction);
        }
        $production_reaction_rates{$string} += $rate(1:$ntime-2); #attribute rates to each parent tagged species and the non-tagged reactions
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        my $string;
        if (defined $parent) {
            $string = $kpp->reaction_string($reaction);
            $string =~ s/_$parent//g; 
            $species_consumption_rates{$species}{$parent}{$string} += $rate;
            $string = $parent; 
        } else {
            $string = $kpp->reaction_string($reaction);
        }
        $consumption_reaction_rates{$string} += $rate(1:$ntime-2);
    }

    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    get_net_rate(\%production_reaction_rates, \%consumption_reaction_rates);
    return (\%production_reaction_rates, \%consumption_reaction_rates);
}

sub get_net_rate {
    my ($production, $consumption) = @_;

    my %common_processes;
    $common_processes{$_} += 1 for (grep {defined $consumption->{$_} } keys %$production);
    foreach my $item (keys %common_processes) {
        my $net_effect = $production->{$item} + $consumption->{$item};
        if ($net_effect->sum > 0) {#overall production
            $production->{$item} += $consumption->{$item};
            delete $consumption->{$item};
        } else { #overall consumption
            $production->{$item} += $production->{$item};
            delete $production->{$item};
        } 
    }
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

sub sort_data_for_plot { #create hash with production of the reactions
    my ($production_rates, $consumption_rates) = @_;
    my %production_rates = %$production_rates;
    my %consumption_rates = %$consumption_rates;
    my (@production_others, @consumption_others, @sorted_plot_data); 
    my $prod_others_max = 4e7;
    my $cons_others_max = -$prod_others_max;

    foreach my $item (keys %consumption_rates) {#sort consumption
        if ($consumption_rates{$item}->sum > $cons_others_max) { #get consumption others
            push @consumption_others, $consumption_rates{$item};
            my $cons_other_rates = cat(@consumption_others);
            $cons_other_rates = $cons_other_rates->xchg(0,1)->sumover;
            $consumption_rates{'Consumption Others'} = $cons_other_rates;
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
            push @production_others, $production_rates{$item};
            my $cons_other_rates = cat(@production_others);
            $cons_other_rates = $cons_other_rates->xchg(0,1)->sumover;
            $production_rates{'Production Others'} = $cons_other_rates;
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

    $R->set('y.max', $args->{y_max});
    $R->set('y.min', $args->{y_min});
    $R->set('y.breaks', $args->{breaks});
    $R->set('time', [@{$args->{times}}]);
    $R->run(q` my.names = c( "CH4" = "Methane", "C2H6" = "Ethane", "C3H8" = "Propane", "NC4H10" = "Butane", "IC4H10" = "2-Methylpropane", "NC5H12" = "Pentane", "IC5H12" = "2-Methylbutane", "NC6H14" = "Hexane", "NC7H16" = "Heptane", "NC8H18" = "Octane", "C2H4" = "Ethene", "C3H6" = "Propene", "BUT1ENE" = "Butene", "MEPROPENE" = "2-Methylpropene", "C5H8" = "Isoprene", "BENZENE" = "Benzene", "TOLUENE" = "Toluene", "MXYL" = "m-Xylene", "OXYL" = "o-Xylene", "PXYL" = "p-Xylene", "EBENZ" = "Ethylbenzene", "BIGALK" = "Pentane", "BIGENE" = "Butene", "ISOP" = "Isoprene", "ISO" = "Isoprene", "ETH" = "Ethane", "HC3" = "Propane", "HC5" = "Pentane", "HC8" = "Octane", "OL2" = "Ethene", "OLT" = "Propene", "OLI" = "MEPROPENE", "TOL" = "Toluene", "XYL" = "m-Xylene", "ETE" = "Ethene", "BEN" = "Benzene", "XYM" = "m-Xylene", "XYO" = "o-Xylene", "XYP" = "p-Xylene") `,
            q` my.colours = c( "Production Others" = "#696537", "NC7H16" = "#f9c600", "EBENZ" = "#76afca", "BENZENE" = "#dc3522", "BEN" = "#dc3522", "OXYL" = "#8c6238", "XYO" = "#8c6238", "PXYL" = "#9bb08f", "XYP" = "#9bb08f", "NC6H14" = "#8b1537", "IC4H10" = "#e7e85e", "C3H6" = "#0352cb", "OLT" = "#0352cb", "C2H6" = "#86b650", "ETH" = "#86b650", "MXYL" = "#6c254f", "XYL" = "#6c254f", "XYM" = "#6c254f", "C5H8" = "#ee6738","ISOP" = "#ee6738",  "ISO" = "#ee6738", "C2H4" = "#e7e85e", "OL2" = "#e7e85e", "ETE" = "#e7e85e", "NC5H12" = "#8ed6d5", "HC5" = "#8ed6d5", "BIGALK" = "#8ed6d5", "C3H8" = "#f3aa7f", "OLT" = "#f3aa7f", "TOLUENE" = "#c65d6c", "TOL" = "#c65d6c", "NC4H10" = "#888a87", "IC5H12" = "#0e5c28", "CH4" = "#b569b3", "H2O2 + hv = OH + OH" = "#2c9def", "HONO + hv = NO + OH" = "#2c9daf", "HONO + hv = OH + NO" = "#2c9daf", "O1D = OH + OH" = "#ae4903", "NO2 + OH = HNO3" = "#0c3f78", "OH + NO2 = HNO3" = "#0c3f78", "HO2 + HO2 = H2O2" = "#f8c56c", "NO + OH = HONO" = "#0e5c28", "OH + NO = HONO" = "#0e5c28", "HO2 + OH = UNITY" = "#f9c500", "OH + HO2 = UNITY" = "#0c3f78", "HNO3 + OH = NO3" = "#a67c52", "Consumption Others" = "#6c254f", "BUT1ENE" = "#86b650", "BIGENE" = "#86b650", "MEPROPENE" = "#ef6638", "OLI" = "#ef6638", "NC8H18" = "#be2448", "HC8" = "#be2448", "C3H8" = "#dc3522", "HC3" = "#dc3522" ) `,
            q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
    );

    #general plot R function
    $R->run(q` mech.plot = function(data, legend.title, plot.title, legend) { sums = colSums(data[,-1]); 
                                                                                rep.no = nrow(data); 
                                                                                regex = as.numeric(gsub("^[A-Z]+ + [A-Z]+", "", sums, perl = TRUE) ); 
                                                                                type = c(); 
                                                                                for ( i in 1:length(sums) ) { type = c(type, rep(regex[i], rep.no)) }; 
                                                                                plot.data = melt(data = data, id = names(data)[1], measured = names(data)[-1] ); 
                                                                                colnames(plot.data) = c("time", "reaction", "rate"); 
                                                                                plot.data$subset.condition = type; 
                                                                                reaction.levels = (levels(factor(plot.data$reaction))); 
                                                                                plot.data$reaction = ordered(plot.data$reaction, levels = reaction.levels); 
                                                                                plot.data = ddply( plot.data, .(reaction)); 
                                                                                cols = colorRampPalette(my.colours)(nlevels(plot.data$reaction)); plot = ggplot(data = plot.data, aes(x = time, y = rate, fill = reaction)); plot = plot + geom_area(data = subset(plot.data, subset.condition > 0), position = "stack", alpha = 1); plot = plot + geom_area(data = subset(plot.data, subset.condition > 0), position = "stack", colour = "black", show_guide = FALSE); plot = plot + geom_area(data = subset(plot.data, subset.condition < 0), position = "stack", alpha = 1); plot = plot + geom_area(data = subset(plot.data, subset.condition < 0), position = "stack", colour = "black", show_guide = FALSE); plot = plot + guides(fill = guide_legend(title = legend.title)); plot = plot + ggtitle(plot.title); plot = plot + scale_x_continuous(limits=c(0, 7), breaks=seq(0, 7, 1)); plot = plot + scale_y_continuous(limits=c(y.min, y.max), breaks=seq(y.min, y.max, y.breaks), label = scientific_10) ; plot = plot + scale_fill_manual( name = "reaction", limits = legend, values = my.colours, labels = my.names); plot = plot + theme_bw() ; plot = plot + theme(legend.key.size = unit(6, "cm")) ; plot = plot +  theme(axis.text.x = element_text(size = 170)) ; plot = plot + theme(axis.text.y = element_text(size = 170)) ; plot = plot + theme(axis.title.x = element_blank()) ; plot = plot + theme(legend.text = element_text(size = 120), legend.title = element_blank()) ; plot = plot + theme(legend.key = element_blank()) ; plot = plot + theme(axis.title.y = element_blank()) ; plot = plot + theme(plot.title = element_text(size = 230, face = "bold", vjust = 1)) ; return(plot) } `);
 
    #MCM v3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.plot.title', $args->{mcm3_1_title});
    $R->set('mcm3.1.legend', $args->{mcm3_1_legend});
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate`); 
        }
    } 
    $R->run(q` mcm3.1.plot = mech.plot(mcm3.1.data, "MCM v3.1", mcm3.1.plot.title, mcm3.1.legend) `);

     #MCM v3.2
     $R->run(q` mcm3.2.data = data.frame(time)`);
     $R->set('mcm3.2.plot.title', $args->{mcm3_2_title});
     $R->set('mcm3.2.legend', $args->{mcm3_2_legend});
     foreach my $ref (@{$args->{mcm3_2_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` mcm3.2.data[name] = rate`); 
         }
     } 
     $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, "MCM v3.2", mcm3.2.plot.title, mcm3.2.legend) `);
 
     #CRI v2
     $R->run(q` cri.data = data.frame(time)`);
     $R->set('cri.plot.title', $args->{cri_title});
     $R->set('cri.legend', $args->{cri_legend});
     foreach my $ref (@{$args->{cri_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` cri.data[name] = rate`); 
         }
     } 
     $R->run(q` cri.plot = mech.plot(cri.data, "CRI v2", cri.plot.title, cri.legend) `);
 
     #MOZART-4
     $R->run(q` mozart.data = data.frame(time)`);
     $R->set('mozart.plot.title', $args->{mozart_title});
     $R->set('mozart.legend', $args->{mozart_legend});
     foreach my $ref (@{$args->{mozart_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` mozart.data[name] = rate`); 
         }
     } 
     $R->run(q` mozart.plot = mech.plot(mozart.data, "MOZART-4", mozart.plot.title, mozart.legend) `);
 
     #RADM2
     $R->run(q` radm2.data = data.frame(time)`);
     $R->set('radm2.plot.title', $args->{radm2_title});
     $R->set('radm2.legend', $args->{radm2_legend});
     foreach my $ref (@{$args->{radm2_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` radm2.data[name] = rate`); 
         }
     } 
     $R->run(q` radm2.plot = mech.plot(radm2.data, "RADM2", radm2.plot.title, radm2.legend) `);
 
     #RACM
     $R->run(q` racm.data = data.frame(time)`);
     $R->set('racm.plot.title', $args->{racm_title});
     $R->set('racm.legend', $args->{racm_legend});
     foreach my $ref (@{$args->{racm_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` racm.data[name] = rate`); 
         }
     } 
     $R->run(q` racm.plot = mech.plot(racm.data, "RACM", racm.plot.title, racm.legend) `);
 
     #RACM2
     $R->run(q` racm2.data = data.frame(time)`);
     $R->set('racm2.plot.title', $args->{racm2_title});
     $R->set('racm2.legend', $args->{racm2_legend});
     foreach my $ref (@{$args->{racm2_data}}) {
         for my $key (keys %$ref) {
             my @value = @{ $ref->{$key} };
             my $R_name = $R->set('name', $key);
             my $R_data = $R->set('rate', [@value]);
             $R->run(q` racm2.data[name] = rate`); 
         }
     } 
     $R->run(q` racm2.plot = mech.plot(racm2.data, "RACM2", racm2.plot.title, racm2.legend) `);
 
    #CBM-IV
    $R->set('cbm4.plot.title', $args->{cbm4_title});
    $R->set('cbm4.legend', $args->{cbm4_legend});
    $R->run(q` cbm4.data = data.frame(time) `);
    foreach my $ref (@{$args->{cbm4_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cbm4.data[name] = rate`); 
        }
    } 
    $R->run(q` cbm4.plot = mech.plot(cbm4.data, "CBM-IV", cbm4.plot.title, cbm4.legend) `);

    #CB05
    $R->set('cb05.plot.title', $args->{cb05_title});
    $R->set('cb05.legend', $args->{cb05_legend});
    $R->run(q` cb05.data = data.frame(time) `);
    foreach my $ref (@{$args->{cb05_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cb05.data[name] = rate`); 
        }
    } 
    $R->run(q` cb05.plot = mech.plot(cb05.data, "CB05", cb05.plot.title, cb05.legend) `);

    $R->run(q` CairoPDF(file = "radical_overall_budgets.pdf", width = 190, height = 200) `, 
            q` y.label = textGrob(expression(bold(paste("\nRate (molecules ", cm^-3, s^-1, ")"))), rot = 90, gp = gpar(fontsize = 210), vjust = 0.6)`,
            q` main.plot = grid.arrange(y.label, 
                                        arrangeGrob(mcm3.2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    mcm3.1.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    cri.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    radm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    racm.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    racm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    mozart.plot, 
                                                    cbm4.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    cb05.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 3), 
                                       nrow = 1, ncol = 2,
                                       sub = textGrob("\nTime (days)\n", gp = gpar(fontsize = 210, fontface = "bold"), vjust = 0.2), 
                                       widths=unit.c(unit(17, "lines"), unit(1, "npc") - unit(17, "lines"))) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );

    $R->stop(); 
}
