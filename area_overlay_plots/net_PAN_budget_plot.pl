#!/usr/bin/perl
# plot net PAN budgets for each mechanism on one plot
# Version 0: Jane Coates 19/02/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA_TIM;
use KPP;
use Statistics::R;

my (%net_rates, %families, %weights);

#MCMv3.1 data
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA_TIM->new($mcm_3_1_run); 
my $mcm_3_1_eqn = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.eqn";
my $mcm_3_1_kpp = KPP->new($mcm_3_1_eqn);
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points
$families{'PAN_mcm3.1'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_NC4H10 PAN_IC4H10 PAN_NC5H12 PAN_IC5H12 PAN_NC6H14 PAN_NC7H16 PAN_NC8H18 PAN_C2H4 PAN_C3H6 PAN_BUT1ENE PAN_MEPROPENE PAN_C5H8 PAN_BENZENE PAN_TOLUENE PAN_MXYL PAN_OXYL PAN_PXYL PAN_EBENZ CH3CO3_C2H6 CH3CO3_C3H8 CH3CO3_NC4H10 CH3CO3_IC4H10 CH3CO3_NC5H12 CH3CO3_IC5H12 CH3CO3_NC6H14 CH3CO3_NC7H16 CH3CO3_NC8H18 CH3CO3_C2H4 CH3CO3_C3H6 CH3CO3_BUT1ENE CH3CO3_MEPROPENE CH3CO3_C5H8 CH3CO3_BENZENE CH3CO3_TOLUENE CH3CO3_MXYL CH3CO3_OXYL CH3CO3_PXYL CH3CO3_EBENZ ) ],
$net_rates{'PAN_mcm3.1'} = get_net_rates( 'PAN_mcm3.1', $mcm_3_1_kpp, $mcm_3_1_mecca );

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA_TIM->new($mcm_3_2_run); 
my $mcm_3_2_eqn = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqn);
$families{'PAN_mcm3.2'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_NC4H10 PAN_IC4H10 PAN_NC5H12 PAN_IC5H12 PAN_NC6H14 PAN_NC7H16 PAN_NC8H18 PAN_C2H4 PAN_C3H6 PAN_BUT1ENE PAN_MEPROPENE PAN_C5H8 PAN_BENZENE PAN_TOLUENE PAN_MXYL PAN_OXYL PAN_PXYL PAN_EBENZ CH3CO3_C2H6 CH3CO3_C3H8 CH3CO3_NC4H10 CH3CO3_IC4H10 CH3CO3_NC5H12 CH3CO3_IC5H12 CH3CO3_NC6H14 CH3CO3_NC7H16 CH3CO3_NC8H18 CH3CO3_C2H4 CH3CO3_C3H6 CH3CO3_BUT1ENE CH3CO3_MEPROPENE CH3CO3_C5H8 CH3CO3_BENZENE CH3CO3_TOLUENE CH3CO3_MXYL CH3CO3_OXYL CH3CO3_PXYL CH3CO3_EBENZ ) ],
$net_rates{'PAN_mcm3.2'} = get_net_rates( 'PAN_mcm3.2', $mcm_3_2_kpp, $mcm_3_2_mecca );

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA_TIM->new($cri_run); 
my $cri_eqn = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqn);
$families{'PAN_cri'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_NC4H10 PAN_IC4H10 PAN_NC5H12 PAN_IC5H12 PAN_NC6H14 PAN_NC7H16 PAN_NC8H18 PAN_C2H4 PAN_C3H6 PAN_BUT1ENE PAN_MEPROPENE PAN_C5H8 PAN_BENZENE PAN_TOLUENE PAN_MXYL PAN_OXYL PAN_PXYL PAN_EBENZ CH3CO3_C2H6 CH3CO3_C3H8 CH3CO3_NC4H10 CH3CO3_IC4H10 CH3CO3_NC5H12 CH3CO3_IC5H12 CH3CO3_NC6H14 CH3CO3_NC7H16 CH3CO3_NC8H18 CH3CO3_C2H4 CH3CO3_C3H6 CH3CO3_BUT1ENE CH3CO3_MEPROPENE CH3CO3_C5H8 CH3CO3_BENZENE CH3CO3_TOLUENE CH3CO3_MXYL CH3CO3_OXYL CH3CO3_PXYL CH3CO3_EBENZ ) ],
$net_rates{'PAN_cri'} = get_net_rates( 'PAN_cri', $cri_kpp, $cri_mecca );

#mozart data
my $mozart_run = "/work/users/jco/MECCA/MOZART_tagging/boxmodel";
my $mozart_mecca = MECCA_TIM->new($mozart_run); 
my $mozart_eqn = "/work/users/jco/MECCA/MOZART_tagging/gas.eqn";
my $mozart_kpp = KPP->new($mozart_eqn);
$families{'PAN_mozart'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_BIGALK PAN_C2H4 PAN_C3H6 PAN_BIGENE PAN_ISOP PAN_TOLUENE CH3CO3_C2H6 CH3CO3_C3H8 CH3CO3_BIGALK CH3CO3_C2H4 CH3CO3_C3H6 CH3CO3_BIGENE CH3CO3_ISOP CH3CO3_TOLUENE ) ],
$net_rates{'PAN_mozart'} = get_net_rates( 'PAN_mozart', $mozart_kpp, $mozart_mecca );

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_mecca = MECCA_TIM->new($radm2_run); 
my $radm2_eqn = "/work/users/jco/MECCA/RADM2_tagged/gas.eqn";
my $radm2_kpp = KPP->new($radm2_eqn);
$families{'PAN_radm2'} = [ qw( PAN_ETH PAN_HC3 PAN_HC5 PAN_HC8 PAN_OL2 PAN_OLT PAN_OLI PAN_ISO PAN_TOL PAN_XYL ACO3_ETH ACO3_HC3 ACO3_HC5 ACO3_HC8 ACO3_OL2 ACO3_OLT ACO3_OLI ACO3_ISO ACO3_TOL ACO3_XYL ) ],
$net_rates{'PAN_radm2'} = get_net_rates( 'PAN_radm2', $radm2_kpp, $radm2_mecca );

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA_TIM->new($racm_run); 
my $racm_eqn = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqn);
$families{'PAN_racm'} = [ qw( PAN_ETH PAN_HC3 PAN_HC5 PAN_HC8 PAN_ETE PAN_OLT PAN_OLI PAN_ISO PAN_TOL PAN_XYL ACO3_ETH ACO3_HC3 ACO3_HC5 ACO3_HC8 ACO3_ETE ACO3_OLT ACO3_OLI ACO3_ISO ACO3_TOL ACO3_XYL ) ],
$net_rates{'PAN_racm'} = get_net_rates( 'PAN_racm', $racm_kpp, $racm_mecca );

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA_TIM->new($racm2_run); 
my $racm2_eqn = "/work/users/jco/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqn);
$families{'PAN_racm2'} = [ qw( PAN_ETH PAN_HC3 PAN_HC5 PAN_HC8 PAN_ETE PAN_OLT PAN_OLI PAN_ISO PAN_BEN PAN_TOL PAN_XYO PAN_XYM PAN_XYP ACO3_ETH ACO3_HC3 ACO3_HC5 ACO3_HC8 ACO3_ETE ACO3_OLT ACO3_OLI ACO3_ISO ACO3_TOL ACO3_BEN ACO3_XYM ACO3_XYO ACO3_XYP ) ],
$net_rates{'PAN_racm2'} = get_net_rates( 'PAN_racm2', $racm2_kpp, $racm2_mecca );

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA_TIM->new($cbm4_run);
my $cbm4_eqn = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqn);
$families{'PAN_cbm4'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_NC4H10 PAN_IC4H10 PAN_NC5H12 PAN_IC5H12 PAN_NC6H14 PAN_NC7H16 PAN_NC8H18 PAN_C2H4 PAN_C3H6 PAN_BUT1ENE PAN_MEPROPENE PAN_C5H8 PAN_BENZENE PAN_TOLUENE PAN_MXYL PAN_OXYL PAN_PXYL PAN_EBENZ C2O3_C2H6 C2O3_C3H8 C2O3_NC4H10 C2O3_IC4H10 C2O3_NC5H12 C2O3_IC5H12 C2O3_NC6H14 C2O3_NC7H16 C2O3_NC8H18 C2O3_C2H4 C2O3_C3H6 C2O3_BUT1ENE C2O3_MEPROPENE C2O3_C5H8 C2O3_BENZENE C2O3_TOLUENE C2O3_MXYL C2O3_OXYL C2O3_PXYL C2O3_EBENZ ) ],
$net_rates{'PAN_cbm4'} = get_net_rates( 'PAN_cbm4', $cbm4_kpp, $cbm4_mecca );

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA_TIM->new($cb05_run);
my $cb05_eqn = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqn);
$families{'PAN_cb05'} = [ qw( PAN_C2H6 PAN_C3H8 PAN_NC4H10 PAN_IC4H10 PAN_NC5H12 PAN_IC5H12 PAN_NC6H14 PAN_NC7H16 PAN_NC8H18 PAN_C2H4 PAN_C3H6 PAN_BUT1ENE PAN_MEPROPENE PAN_C5H8 PAN_BENZENE PAN_TOLUENE PAN_MXYL PAN_OXYL PAN_PXYL PAN_EBENZ C2O3_C2H6 C2O3_C3H8 C2O3_NC4H10 C2O3_IC4H10 C2O3_NC5H12 C2O3_IC5H12 C2O3_NC6H14 C2O3_NC7H16 C2O3_NC8H18 C2O3_C2H4 C2O3_C3H6 C2O3_BUT1ENE C2O3_MEPROPENE C2O3_C5H8 C2O3_BENZENE C2O3_TOLUENE C2O3_MXYL C2O3_OXYL C2O3_PXYL C2O3_EBENZ ) ],
$net_rates{'PAN_cb05'} = get_net_rates( 'PAN_cb05', $cb05_kpp, $cb05_mecca );

#Create x-axis for plot in hours
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list;

net_budget_plot(\@time_axis, \%net_rates);

sub get_net_rates {
    my ($species, $kpp, $mecca) = @_;

    my (@rate_array_producers, @rate_array_consumers);
    my ($consumers, $producers, $consumer_yields, $producer_yields);
    if (exists $families{$species}) { #get family reaction numbers and yields
        $kpp->family({ 
                name    => $species,
                members => $families{$species}, 
                weights => $weights{$species},
        });
        $consumers = $kpp->consuming($species);
        $producers = $kpp->producing($species);
        $consumer_yields = $kpp->effect_on($species, $consumers); 
        $producer_yields = $kpp->effect_on($species, $producers);  
    } else { #get reaction numbers and yields
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
        push @rate_array_producers, $rate(1:$ntime-2);
    }
    
    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        push @rate_array_consumers, $rate(1:$ntime-2);
    }

    #Create piddle of all rates from piddle list
    my $rates_producer = cat(@rate_array_producers);
    my $rates_consumer = cat(@rate_array_consumers);

    #Calculate total production/consumption
    my $consumer_rate_sum = $rates_consumer->xchg(0,1)->sumover; 
    my $producer_rate_sum = $rates_producer->xchg(0,1)->sumover;

    #Calculate net production and output in array for plotting with R 
    my $net_rate_pdl = $producer_rate_sum + $consumer_rate_sum;
    my $net_rate_list = join ":", $net_rate_pdl->dog;
    my @net_rates = split /:/, $net_rate_list;

    return \@net_rates;
}

sub net_budget_plot {
    my ($time, $data) = @_;
    my %data = %$data;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(gtable) `);
    $R->run(q` library(Cairo) `);
    $R->set('time', [@$time]);   
    $R->set('rep.number', 9); #number of mechanisms

    $R->run(q` times = rep(time, rep.number) `,
            q` Species = {} `,
            q` Mechanism = {} `,
            q` Net.rates = {} `,
    );

    foreach my $item (sort keys %data) {
        my ($name, $mechanism) = split /_/, $item;
        my $R_name = $R->set('name', $name);
        if ($mechanism =~ /mcm3\.1/) {
            my $R_mech = $R->set('mechanism', 'MCM v3.1');
        } elsif ($mechanism =~ /mcm3\.2/) {
            my $R_mech = $R->set('mechanism', 'MCM v3.2');
        } elsif ($mechanism =~ /cri/) {
            my $R_mech = $R->set('mechanism', 'CRI v2');
        } elsif ($mechanism =~ /mozart/) {
            my $R_mech = $R->set('mechanism', 'MOZART-4');
        } elsif ($mechanism =~ /radm2/) {
            my $R_mech = $R->set('mechanism', 'RADM2');
        } elsif ($mechanism =~ /racm2/) {
            my $R_mech = $R->set('mechanism', 'RACM2');
        } elsif ($mechanism =~ /racm/) {
            my $R_mech = $R->set('mechanism', 'RACM');
        } elsif ($mechanism =~ /cbm4/) {
            my $R_mech = $R->set('mechanism', 'CBM-IV');
        } elsif ($mechanism =~ /cb05/) {
            my $R_mech = $R->set('mechanism', 'CB05');
        }
        my $R_data = $R->set('net.rates', [@{$data{$item}}]);
        $R->run(q` Species = cbind(Species, rep(name, length(time))) `,
                q` Mechanism = cbind(Mechanism, rep(mechanism, length(time))) `,
                q` Net.rates = cbind(Net.rates, net.rates) `,
        );
    }

    #create dataframe after converting the matrices above to vectors
    $R->run(q` Species = c(Species) `,
            q` Mechanism = c(Mechanism) `,
            q` Net.rates = c(Net.rates) `,
            q` data = data.frame(times, Species, Mechanism, Net.rates) `,
    );
    
    $R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `), #scientific label format for y-axis
    $R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);
    
    #plot
    $R->run(q` plot.lines = function () { list( geom_line(size = 3), scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), scale_y_continuous(limits = c(-2e6, 6e6), breaks = seq(-2e6, 6e6, 1e6), labels = scientific_10), theme_bw(), theme(axis.text.x = element_text(size = 50)), theme(axis.text.y = element_text(size = 50)), theme(legend.title = element_text(size = 80, face = "bold")), theme(legend.text = theme_text(size = 70)), theme(legend.key.size = unit(5, "cm")), theme(axis.title.x = element_text(size = 80, face = "bold", vjust = 0.2)), theme(axis.title.y = element_text(size = 80, face = "bold", vjust = 0.5)), scale_colour_manual(values = my.colours), theme(legend.key = element_blank()), theme(legend.position = c(0.8, 0.7)) ) } `);

    $R->run(q` plot1 = ggplot(data, aes(x = times, y = Net.rates, colour = Mechanism)) `,
            q` plot1 = plot1 + plot.lines() `, 
            q` plot1 = plot1 + ylab(expression(bold(paste("\nRate (molecules ", cm^-3, s^-1, ")")))) `,
            q` plot1 = plot1 + xlab("\nTime (days)\n") `,
            q` CairoPDF(file = "net_PAN_budget.pdf", width = 50, height = 40) `,
            q` print(plot1) `,
            q` dev.off() `,
    );

    $R->stop();
} 
