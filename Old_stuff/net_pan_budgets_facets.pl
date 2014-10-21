#!/usr/bin/perl
# Overall PAN budgets as stacked bar plots
# Version 0 : Jane Coates 3/07/2014 

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
my $mcm_3_1_panfile = "/work/users/jco/Mechanisms/MCM/MCM-3.1/All_pans.txt";
my @mcm_3_1_pan_family = get_pans($mcm_3_1_panfile); 
$families{'PAN_mcm_3_1'} = [ @mcm_3_1_pan_family ]; 
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points 
($production_rates{'PAN_mcm_3_1'}, $consumption_rates{'PAN_mcm_3_1'}) = get_rates('PAN_mcm_3_1', $mcm_3_1_kpp, $mcm_3_1_mecca);
my ($mcm_3_1_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_mcm_3_1'}, $consumption_rates{'PAN_mcm_3_1'});
my $mcm_3_1_plot_title = "(b) MCM v3.1";

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_panfile = "/work/users/jco/Mechanisms/MCM/MCM-3.2/All_pans.txt";
my @mcm_3_2_pan_family = get_pans($mcm_3_2_panfile); 
$families{'PAN_mcm_3_2'} = [ @mcm_3_2_pan_family ]; 
($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'}) = get_rates('PAN_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run); 
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile); 
my $cri_panfile = "/work/users/jco/Mechanisms/CRI/CRI_v2_full/All_pans.txt";
my @cri_pan_family = get_pans($cri_panfile); 
$families{'PAN_cri'} = [ @cri_pan_family ]; 
($production_rates{'PAN_cri'}, $consumption_rates{'PAN_cri'}) = get_rates('PAN_cri', $cri_kpp, $cri_mecca); 
my ($cri_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_cri'}, $consumption_rates{'PAN_cri'});
my $cri_plot_title = "(c) CRI v2";

#mozart data
my $mozart_run = "/work/users/jco/MECCA/MOZART_tagging/boxmodel";
my $mozart_mecca = MECCA->new($mozart_run); 
my $mozart_eqnfile = "/work/users/jco/MECCA/MOZART_tagging/gas.eqn";
my $mozart_kpp = KPP->new($mozart_eqnfile); 
my $mozart_panfile = "/work/users/jco/Mechanisms/MOZART/MOZART/All_pans.txt";
my @mozart_pan_family = get_pans($mozart_panfile); 
$families{'PAN_mozart'} = [ @mozart_pan_family ]; 
($production_rates{'PAN_mozart'}, $consumption_rates{'PAN_mozart'}) = get_rates('PAN_mozart', $mozart_kpp, $mozart_mecca);
my ($mozart_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_mozart'}, $consumption_rates{'PAN_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_mecca = MECCA->new($radm2_run); 
my $radm2_eqnfile = "/work/users/jco/MECCA/RADM2_tagged/gas.eqn";
my $radm2_kpp = KPP->new($radm2_eqnfile); 
my $radm2_panfile = "/work/users/jco/Mechanisms/RADM2/All_pans.txt";
my @radm2_pan_family = get_pans($radm2_panfile); 
$families{'PAN_radm2'} = [ @radm2_pan_family ]; 
($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'}) = get_rates('PAN_radm2', $radm2_kpp, $radm2_mecca);
my ($radm2_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA->new($racm_run); 
my $racm_eqnfile = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqnfile); 
my $racm_panfile = "/work/users/jco/Mechanisms/RACM/All_pans.txt";
my @racm_pan_family = get_pans($racm_panfile); 
$families{'PAN_racm'} = [ @racm_pan_family ]; 
($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'}) = get_rates('PAN_racm', $racm_kpp, $racm_mecca);
my ($racm_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'});
my $racm_plot_title = "(e) RACM";

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA->new($racm2_run); 
my $racm2_eqnfile = "/work/users/jco/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqnfile); 
my $racm2_panfile = "/work/users/jco/Mechanisms/RACM2/All_pans.txt";
my @racm2_pan_family = get_pans($racm2_panfile); 
$families{'PAN_racm2'} = [ @racm2_pan_family ]; 
($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'}) = get_rates('PAN_racm2', $racm2_kpp, $racm2_mecca);
my ($racm2_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'});
my $racm2_plot_title = "(f) RACM2";

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA->new($cbm4_run);
my $cbm4_eqnfile = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqnfile); 
my $cbm4_panfile = "/work/users/jco/Mechanisms/CBM-IV/All_pans.txt";
my @cbm4_pan_family = get_pans($cbm4_panfile); 
$families{'PAN_cbm4'} = [ @cbm4_pan_family ]; 
($production_rates{'PAN_cbm4'}, $consumption_rates{'PAN_cbm4'}) = get_rates('PAN_cbm4', $cbm4_kpp, $cbm4_mecca);
my ($cbm4_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_cbm4'}, $consumption_rates{'PAN_cbm4'});
my $cbm4_plot_title = "(h) CBM-IV";

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA->new($cb05_run);
my $cb05_eqnfile = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqnfile); 
my $cb05_panfile = "/work/users/jco/Mechanisms/CB05/All_pans.txt";
my @cb05_pan_family = get_pans($cb05_panfile); 
$families{'PAN_cb05'} = [ @cb05_pan_family ]; 
($production_rates{'PAN_cb05'}, $consumption_rates{'PAN_cb05'}) = get_rates('PAN_cb05', $cb05_kpp, $cb05_mecca);
my ($cb05_sorted_plot_data) = sort_data_for_plot($production_rates{'PAN_cb05'}, $consumption_rates{'PAN_cb05'});
my $cb05_plot_title = "(i) CB05";

#Create x-axis for plot in hours
my $times = $mcm_3_1_mecca->time;
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
        mcm3_1_data     => $mcm_3_1_sorted_plot_data,
        mcm3_1_title    => $mcm_3_1_plot_title,
        mcm3_2_data     => $mcm_3_2_sorted_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title, 
        cri_data        => $cri_sorted_plot_data,
        cri_title       => $cri_plot_title,
        mozart_data     => $mozart_sorted_plot_data,
        mozart_title    => $mozart_plot_title,
        radm2_data      => $radm2_sorted_plot_data,
        radm2_title     => $radm2_plot_title,
        racm_data       => $racm_sorted_plot_data,
        racm_title      => $racm_plot_title,
        racm2_data      => $racm2_sorted_plot_data,
        racm2_title     => $racm2_plot_title,
        cbm4_data       => $cbm4_sorted_plot_data,
        cbm4_title      => $cbm4_plot_title,
        cb05_data       => $cb05_sorted_plot_data,
        cb05_title      => $cb05_plot_title,
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
            if ($parent =~ /CH4|C2H6|C3H8|NC4H10|IC4H10|NC5H12|IC5H12|NC6H14|NC7H16|NC8H18|BIGALK|ETH|HC3|HC5|HC8/) {
                $string = 'Alkanes Production'; 
            } elsif ($parent =~ /C2H4|C3H6|BUT1ENE|MEPROPENE|OL2|OLT|OLI|C5H8|ISO|ETE|BIGENE|ISOP/) {
                $string = 'Alkenes Production';
            } elsif ($parent =~ /BENZENE|TOLUENE|MXYL|PXYL|OXYL|EBENZ|TOL|XYL|BEN|XYM|XYO|XYP/) {
                $string = 'Aromatics Production';
            }
            $production_reaction_rates{$string} += $rate(1:$ntime-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        my $string;
        if (defined $parent) {
            if ($parent =~ /CH4|C2H6|C3H8|NC4H10|IC4H10|NC5H12|IC5H12|NC6H14|NC7H16|NC8H18|BIGALK|ETH|HC3|HC5|HC8/) {
                $string = 'Alkanes Consumption'; 
            } elsif ($parent =~ /C2H4|C3H6|BUT1ENE|MEPROPENE|OL2|OLT|OLI|C5H8|ISO|ETE|BIGENE|ISOP/) {
                $string = 'Alkenes Consumption';
            } elsif ($parent =~ /BENZENE|TOLUENE|MXYL|PXYL|OXYL|EBENZ|TOL|XYL|BEN|XYM|XYO|XYP/) {
                $string = 'Aromatics Consumption';
            }
            $consumption_reaction_rates{$string} += $rate(1:$ntime-2);
        }
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
    my (@sorted_plot_data, @plot_data); 

    push @sorted_plot_data, { $_ => $consumption_rates->{$_} } foreach (sort keys %$consumption_rates);
    push @sorted_plot_data, { $_ => $production_rates->{$_} } foreach (sort keys %$production_rates);

    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            my @rate_array = map { $_ } $ref->{$item}->dog; 
            push @plot_data, { $item => \@rate_array };
        }
    } 

    return (\@plot_data);
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
    $R->run(q` my.colours = c( "#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20" ) `,
            q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
    );

    #general plot R function
    $R->run(q` arrange.data = function(data, flag) {   data = ddply(data, .(time), colwise(sum)) ;
                                                       data$"Alkanes" = data$"Alkanes Production" + data$"Alkanes Consumption" ;
                                                       data$"Alkenes" = data$"Alkenes Production" + data$"Alkenes Consumption" ;
                                                       data$"Aromatics" = data$"Aromatics Production" + data$"Aromatics Consumption" ;
                                                       subset = data[1:7, c("time", "Alkanes", "Alkenes", "Aromatics")];
                                                       plot.data = melt(data = subset, id = names(subset)[1], measured = names(subset)[-1] ); 
                                                       colnames(plot.data) = c("time", "group", "rate"); 
                                                       plot.data$"Mechanism" = rep(flag, 21) ;
                                                       return(plot.data) } `);

    #MCM v3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.flag', "MCMv3.1");
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate`); 
        }
    } 
    $R->run(q` mcm3.1.plot.data = arrange.data(mcm3.1.data, mcm3.1.flag) `); 

     #MCM v3.2
     $R->run(q` mcm3.2.data = data.frame(time)`);
     $R->set('mcm3.2.flag', "MCMv3.2");
     foreach my $ref (@{$args->{mcm3_2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` mcm3.2.data[name] = rate`); 
         }
     } 
     $R->run(q` mcm3.2.plot.data = arrange.data(mcm3.2.data, mcm3.2.flag) `);

     #CRI v2
     $R->run(q` cri.data = data.frame(time)`);
     $R->set('cri.flag', "CRIv2");
     foreach my $ref (@{$args->{cri_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` cri.data[name] = rate`); 
         }
     } 
     $R->run(q` cri.plot.data = arrange.data(cri.data, cri.flag) `);

     #MOZART-4
     $R->run(q` mozart.data = data.frame(time)`);
     $R->set('mozart.flag', "MOZART-4");
     foreach my $ref (@{$args->{mozart_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` mozart.data[name] = rate`); 
         }
     } 
     $R->run(q` mozart.plot.data = arrange.data(mozart.data, mozart.flag) `);
 
     #RADM2
     $R->run(q` radm2.data = data.frame(time)`);
     $R->set('radm2.flag', "RADM2");
     foreach my $ref (@{$args->{radm2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` radm2.data[name] = rate`); 
         }
     } 
     $R->run(q` radm2.plot.data = arrange.data(radm2.data, radm2.flag) `);
 
     #RACM
     $R->run(q` racm.data = data.frame(time)`);
     $R->set('racm.flag', "RACM");
     foreach my $ref (@{$args->{racm_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` racm.data[name] = rate`); 
         }
     } 
     $R->run(q` racm.plot.data = arrange.data(racm.data, racm.flag) `);
 
     #RACM2
     $R->run(q` racm2.data = data.frame(time)`);
     $R->set('racm2.flag', "RACM2");
    foreach my $ref (@{$args->{racm2_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` racm2.data[name] = rate`); 
         }
     } 
     $R->run(q` racm2.plot.data = arrange.data(racm2.data, racm2.flag) `);
 
     #CBM-IV
    $R->run(q` cbm4.data = data.frame(time) `);
     $R->set('cbm4.flag', "CBM-IV");
    foreach my $ref (@{$args->{cbm4_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cbm4.data[name] = rate`); 
        }
    } 
     $R->run(q` cbm4.plot.data = arrange.data(cbm4.data, cbm4.flag) `);

    #CB05
    $R->run(q` cb05.data = data.frame(time) `);
     $R->set('cb05.flag', "CB05");
    foreach my $ref (@{$args->{cb05_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cb05.data[name] = rate`); 
        }
    } 
     $R->run(q` cb05.plot.data = arrange.data(cb05.data, cb05.flag) `);

     $R->run(q` plot.data = rbind(mcm3.1.plot.data, mcm3.2.plot.data, cri.plot.data, mozart.plot.data, radm2.plot.data, racm.plot.data, racm2.plot.data, cbm4.plot.data, cb05.plot.data) `);
     $R->run(q` facet_labels = list( 'Alkanes' = '(a) Alkanes', 'Alkenes' = '(b) Alkenes', 'Aromatics' = '(c) Aromatics') `);
    $R->run(q` facet_labeller = function(variable,value) {return(facet_labels[value])} `); 

    $R->run(q` plot.lines = function () { list(ylab(expression(bold(paste("\nNet Production/Loss Rate (molecules ", cm^-3, s^-1, ")")))), 
                                               geom_line(size = 3),
                                               geom_point(size = 6),
                                               facet_grid( ~ group, labeller = facet_labeller),
                                               theme_bw(), 
                                               scale_y_continuous(limits = c(-2e7, 4e7), breaks = seq(-2e7, 4e7, 1e7)),
                                               theme(axis.title.x = element_blank()), 
                                               theme(axis.title.y = element_text(size = 30, face = "bold")), 
                                               theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), legend.title = element_text(size = 30, face = "bold"), legend.key.size = unit(2, "cm"), legend.text = element_text(size = 25), legend.key = element_blank(), plot.title = element_text(size = 40, face = "bold"), strip.text.x = element_text(size = 30, face = "bold"), strip.background = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()), 
                                               scale_colour_manual(values = my.colours) ) } `); 

    $R->run(q` plot = ggplot(data = plot.data, aes(x = time, y = rate, colour = Mechanism, group = Mechanism)) `,
            q` plot = plot + plot.lines() `,
            q` CairoPDF(file = "Net_PAN_budgets.pdf", width = 30, height = 10) `,
            q` print(plot) `,
            q` dev.off() `,
    );
 

    $R->stop(); 
}
