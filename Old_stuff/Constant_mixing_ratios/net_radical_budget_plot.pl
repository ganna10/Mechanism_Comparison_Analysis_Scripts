#!/usr/bin/perl
# plot net radical budgets for each mechanism on one plot
# Version 0: Jane Coates 18/02/2014

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
my $mcm_3_1_ro2file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/RO2_species.txt";
my @mcm_3_1_no2_reservoirs = get_no2_reservoirs($mcm_3_1_kpp, $mcm_3_1_ro2file);
my $mcm_3_1_radical_file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/radicals.txt";
my @mcm_3_1_radicals = get_radicals($mcm_3_1_radical_file);
$families{'radicals_mcm_3.1'} = [ @mcm_3_1_radicals, @mcm_3_1_no2_reservoirs, 'HO2NO2' ];
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points
$net_rates{'radicals_mcm3.1'} = get_net_rates( 'radicals_mcm_3.1', $mcm_3_1_kpp, $mcm_3_1_mecca );

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA_TIM->new($mcm_3_2_run); 
my $mcm_3_2_eqn = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqn);
my $mcm_3_2_ro2file = "/work/users/jco/MECCA/MCM_3.2_tagged/RO2_species.txt";
my @mcm_3_2_no2_reservoirs = get_no2_reservoirs($mcm_3_2_kpp, $mcm_3_2_ro2file);
my $mcm_3_2_radical_file = "/work/users/jco/MECCA/MCM_3.2_tagged/radicals.txt";
my @mcm_3_2_radicals = get_radicals($mcm_3_2_radical_file);
$families{'radicals_mcm_3.2'} = [ @mcm_3_2_radicals, @mcm_3_2_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_mcm3.2'} = get_net_rates( 'radicals_mcm_3.2', $mcm_3_2_kpp, $mcm_3_2_mecca );

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA_TIM->new($cri_run); 
my $cri_eqn = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqn);
my $cri_ro2file = "/work/users/jco/MECCA/CRI_tagging/RO2_species.txt";
my @cri_no2_reservoirs = get_no2_reservoirs($cri_kpp, $cri_ro2file);
my $cri_radical_file = "/work/users/jco/MECCA/CRI_tagging/radicals.txt";
my @cri_radicals = get_radicals($cri_radical_file);
$families{'radicals_cri'} = [ @cri_radicals, @cri_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_cri'} = get_net_rates( 'radicals_cri', $cri_kpp, $cri_mecca );

#mozart data
my $mozart_run = "/work/users/jco/MECCA/MOZART_tagging/boxmodel";
my $mozart_mecca = MECCA_TIM->new($mozart_run); 
my $mozart_eqn = "/work/users/jco/MECCA/MOZART_tagging/gas.eqn";
my $mozart_kpp = KPP->new($mozart_eqn);
my $mozart_ro2file = "/work/users/jco/MECCA/MOZART_tagging/RO2_species.txt";
my @mozart_no2_reservoirs = get_no2_reservoirs($mozart_kpp, $mozart_ro2file);
my $mozart_radical_file = "/work/users/jco/MECCA/MOZART_tagging/radicals.txt";
my @mozart_radicals = get_radicals($mozart_radical_file);
$families{'radicals_mozart'} = [ @mozart_radicals, @mozart_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_mozart'} = get_net_rates( 'radicals_mozart', $mozart_kpp, $mozart_mecca );

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_mecca = MECCA_TIM->new($radm2_run); 
my $radm2_eqn = "/work/users/jco/MECCA/RADM2_tagged/gas.eqn";
my $radm2_kpp = KPP->new($radm2_eqn);
my $radm2_ro2file = "/work/users/jco/MECCA/RADM2_tagged/RO2_species.txt";
my @radm2_no2_reservoirs = get_no2_reservoirs($radm2_kpp, $radm2_ro2file);
my $radm2_radical_file = "/work/users/jco/MECCA/RADM2_tagged/radicals.txt";
my @radm2_radicals = get_radicals($radm2_radical_file);
$families{'radicals_radm2'} = [ @radm2_radicals, @radm2_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_radm2'} = get_net_rates( 'radicals_radm2', $radm2_kpp, $radm2_mecca );

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_mecca = MECCA_TIM->new($racm_run); 
my $racm_eqn = "/work/users/jco/MECCA/RACM_tagging/gas.eqn";
my $racm_kpp = KPP->new($racm_eqn);
my $racm_ro2file = "/work/users/jco/MECCA/RACM_tagging/RO2_species.txt";
my @racm_no2_reservoirs = get_no2_reservoirs($racm_kpp, $racm_ro2file);
my $racm_radical_file = "/work/users/jco/MECCA/RACM_tagging/radicals.txt";
my @racm_radicals = get_radicals($racm_radical_file);
$families{'radicals_racm'} = [ @racm_radicals, @racm_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_racm'} = get_net_rates( 'radicals_racm', $racm_kpp, $racm_mecca );

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA_TIM->new($racm2_run); 
my $racm2_eqn = "/work/users/jco/MECCA/RACM2_tagged/gas.eqn";
my $racm2_kpp = KPP->new($racm2_eqn);
my $racm2_ro2file = "/work/users/jco/MECCA/RACM2_tagged/RO2_species.txt";
my @racm2_no2_reservoirs = get_no2_reservoirs($racm2_kpp, $racm2_ro2file);
my $racm2_radical_file = "/work/users/jco/MECCA/RACM2_tagged/radicals.txt";
my @racm2_radicals = get_radicals($racm2_radical_file);
$families{'radicals_racm2'} = [ @racm2_radicals, @racm2_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_racm2'} = get_net_rates( 'radicals_racm2', $racm2_kpp, $racm2_mecca );

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA_TIM->new($cbm4_run);
my $cbm4_eqn = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqn);
my $cbm4_ro2file = "/work/users/jco/MECCA/CBM4_tagging/RO2_species.txt";
my @cbm4_no2_reservoirs = get_no2_reservoirs($cbm4_kpp, $cbm4_ro2file);
my $cbm4_radical_file = "/work/users/jco/MECCA/CBM4_tagging/radicals.txt";
my @cbm4_radicals = get_radicals($cbm4_radical_file);
$families{'radicals_cbm4'} = [ @cbm4_radicals, @cbm4_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_cbm4'} = get_net_rates( 'radicals_cbm4', $cbm4_kpp, $cbm4_mecca );

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA_TIM->new($cb05_run);
my $cb05_eqn = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqn);
my $cb05_ro2file = "/work/users/jco/MECCA/CB05_tagging/RO2_species.txt";
my @cb05_no2_reservoirs = get_no2_reservoirs($cb05_kpp, $cb05_ro2file);
my $cb05_radical_file = "/work/users/jco/MECCA/CB05_tagging/radicals.txt";
my @cb05_radicals = get_radicals($cb05_radical_file);
$families{'radicals_cb05'} = [ @cb05_radicals, @cb05_no2_reservoirs, 'HO2NO2' ];
$net_rates{'radicals_cb05'} = get_net_rates( 'radicals_cb05', $cb05_kpp, $cb05_mecca );

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

sub get_radicals { #get radicals from file
    my ($file) = @_; 

    open FILE, $file or die $!; 
    my @radicals;
    for (<FILE>) {
        push @radicals, split /\s+/, $_; 
    }
    close FILE;
    return @radicals;
} 

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
    $R->run(q` plot.lines = function () { list( geom_line(size = 3), scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), scale_y_continuous(limits = c(-2e6, 9e6), breaks = seq(-2e6, 9e6, 1e6), labels = scientific_10), theme_bw(), theme(axis.text.x = element_text(size = 50)), theme(axis.text.y = element_text(size = 50)), theme(legend.title = element_text(size = 80, face = "bold")), theme(legend.text = theme_text(size = 70)), theme(legend.key.size = unit(5, "cm")), theme(axis.title.x = element_text(size = 80, face = "bold", vjust = 0.2)), theme(axis.title.y = element_text(size = 80, face = "bold", vjust = 0.5)), scale_colour_manual(values = my.colours), theme(legend.key = element_blank()), theme(legend.position = c(0.8, 0.7)) ) } `);

    $R->run(q` plot1 = ggplot(data, aes(x = times, y = Net.rates, colour = Mechanism)) `,
            q` plot1 = plot1 + plot.lines() `, 
            q` plot1 = plot1 + ylab(expression(bold(paste("\nRate (molecules ", cm^-3, s^-1, ")")))) `,
            q` plot1 = plot1 + xlab("\nTime (days)\n") `,
            q` CairoPDF(file = "net_radical_budget.pdf", width = 50, height = 40) `,
            q` print(plot1) `,
            q` dev.off() `,
    );

    $R->stop();
} 
