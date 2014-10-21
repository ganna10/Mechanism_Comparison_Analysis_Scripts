#! /usr/bin/perl
# Analyse radical budget during Pentane breakdown, compare for all mechanisms
# Version 0: Jane Coates 17/02/2014

use strict;
use diagnostics;
use KPP;
use MECCA_TIM;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);
my @species = qw( NC5H12 BIGALK HC5 );

#mcm 3.1
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA_TIM->new($mcm_3_1_run); 
my $mcm_3_1_eqnfile = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.eqn";
my $mcm_3_1_kpp = KPP->new($mcm_3_1_eqnfile); 
my $mcm_3_1_ro2file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/RO2_species.txt";
my @mcm_3_1_no2_reservoirs = get_no2_reservoirs($mcm_3_1_kpp, $mcm_3_1_ro2file);
my $mcm_3_1_radical_file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/radicals.txt";
my @mcm_3_1_radicals = get_radicals($mcm_3_1_radical_file);
$families{'radicals_mcm_3_1'} = [ @mcm_3_1_radicals, @mcm_3_1_no2_reservoirs, 'HO2NO2' ];
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points 

($production_rates{'radicals_mcm_3_1'}, $consumption_rates{'radicals_mcm_3_1'}) = get_rates($species[0], 'radicals_mcm_3_1', $mcm_3_1_kpp, $mcm_3_1_mecca);
my ($mcm_3_1_sorted_plot_data, $mcm_3_1_legend) = sort_data_for_plot($production_rates{'radicals_mcm_3_1'}, $consumption_rates{'radicals_mcm_3_1'});
my $mcm_3_1_plot_title = "(b) MCM v3.1";

#mcm 3.2
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA_TIM->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_ro2file = "/work/users/jco/MECCA/MCM_3.2_tagged/RO2_species.txt";
my @mcm_3_2_no2_reservoirs = get_no2_reservoirs($mcm_3_2_kpp, $mcm_3_2_ro2file);
my $mcm_3_2_radical_file = "/work/users/jco/MECCA/MCM_3.2_tagged/radicals.txt";
my @mcm_3_2_radicals = get_radicals($mcm_3_2_radical_file);
$families{'radicals_mcm_3_2'} = [ @mcm_3_2_radicals, @mcm_3_2_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_mcm_3_2'}, $consumption_rates{'radicals_mcm_3_2'}) = get_rates($species[0], 'radicals_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($production_rates{'radicals_mcm_3_2'}, $consumption_rates{'radicals_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#cri
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA_TIM->new($cri_run);
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile);
my $cri_ro2file = "/work/users/jco/MECCA/CRI_tagging/RO2_species.txt";
my @cri_no2_reservoirs = get_no2_reservoirs($cri_kpp, $cri_ro2file);
my $cri_radical_file = "/work/users/jco/MECCA/CRI_tagging/radicals.txt";
my @cri_radicals = get_radicals($cri_radical_file);
$families{'radicals_cri'} = [ @cri_radicals, @cri_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_cri'}, $consumption_rates{'radicals_cri'}) = get_rates($species[0], 'radicals_cri', $cri_kpp, $cri_mecca); 
my ($cri_sorted_plot_data, $cri_legend) = sort_data_for_plot($production_rates{'radicals_cri'}, $consumption_rates{'radicals_cri'});
my $cri_plot_title = "(c) CRI v2";

#mozart
my $mozart_run = '/work/users/jco/MECCA/MOZART_tagging/boxmodel';
my $mozart_mecca = MECCA_TIM->new($mozart_run);
my $mozart_eqnfile = '/work/users/jco/MECCA/MOZART_tagging/gas.eqn';
my $mozart_kpp = KPP->new($mozart_eqnfile);
my $mozart_ro2file = "/work/users/jco/MECCA/MOZART_tagging/RO2_species.txt";
my @mozart_no2_reservoirs = get_no2_reservoirs($mozart_kpp, $mozart_ro2file);
my $mozart_radical_file = "/work/users/jco/MECCA/MOZART_tagging/radicals.txt";
my @mozart_radicals = get_radicals($mozart_radical_file);
$families{'radicals_mozart'} = [ @mozart_radicals, @mozart_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_mozart'}, $consumption_rates{'radicals_mozart'}) = get_rates($species[1], 'radicals_mozart', $mozart_kpp, $mozart_mecca);
my ($mozart_sorted_plot_data, $mozart_legend) = sort_data_for_plot($production_rates{'radicals_mozart'}, $consumption_rates{'radicals_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#radm2
my $radm2_run = '/work/users/jco/MECCA/RADM2_tagged/boxmodel';
my $radm2_mecca = MECCA_TIM->new($radm2_run);
my $radm2_eqnfile = '/work/users/jco/MECCA/RADM2_tagged/gas.eqn';
my $radm2_kpp = KPP->new($radm2_eqnfile);
my $radm2_ro2file = "/work/users/jco/MECCA/RADM2_tagged/RO2_species.txt";
my @radm2_no2_reservoirs = get_no2_reservoirs($radm2_kpp, $radm2_ro2file);
my $radm2_radical_file = "/work/users/jco/MECCA/RADM2_tagged/radicals.txt";
my @radm2_radicals = get_radicals($radm2_radical_file);
$families{'radicals_radm2'} = [ @radm2_radicals, @radm2_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_radm2'}, $consumption_rates{'radicals_radm2'}) = get_rates($species[2], 'radicals_radm2', $radm2_kpp, $radm2_mecca);
my ($radm2_sorted_plot_data, $radm2_legend) = sort_data_for_plot($production_rates{'radicals_radm2'}, $consumption_rates{'radicals_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm
my $racm_run = '/work/users/jco/MECCA/RACM_tagging/boxmodel';
my $racm_mecca = MECCA_TIM->new($racm_run);
my $racm_eqnfile = '/work/users/jco/MECCA/RACM_tagging/gas.eqn';
my $racm_kpp = KPP->new($racm_eqnfile);
my $racm_ro2file = "/work/users/jco/MECCA/RACM_tagging/RO2_species.txt";
my @racm_no2_reservoirs = get_no2_reservoirs($racm_kpp, $racm_ro2file);
my $racm_radical_file = "/work/users/jco/MECCA/RACM_tagging/radicals.txt";
my @racm_radicals = get_radicals($racm_radical_file);
$families{'radicals_racm'} = [ @racm_radicals, @racm_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_racm'}, $consumption_rates{'radicals_racm'}) = get_rates($species[2], 'radicals_racm', $racm_kpp, $racm_mecca);
my ($racm_sorted_plot_data, $racm_legend) = sort_data_for_plot($production_rates{'radicals_racm'}, $consumption_rates{'radicals_racm'});
my $racm_plot_title = "(e) RACM";

#racm2
my $racm2_run = '/work/users/jco/MECCA/RACM2_tagged/boxmodel';
my $racm2_mecca = MECCA_TIM->new($racm2_run);
my $racm2_eqnfile = '/work/users/jco/MECCA/RACM2_tagged/gas.eqn';
my $racm2_kpp = KPP->new($racm2_eqnfile);
my $racm2_ro2file = "/work/users/jco/MECCA/RACM2_tagged/RO2_species.txt";
my @racm2_no2_reservoirs = get_no2_reservoirs($racm2_kpp, $racm2_ro2file);
my $racm2_radical_file = "/work/users/jco/MECCA/RACM2_tagged/radicals.txt";
my @racm2_radicals = get_radicals($racm2_radical_file);
$families{'radicals_racm2'} = [ @racm2_radicals, @racm2_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_racm2'}, $consumption_rates{'radicals_racm2'}) = get_rates($species[2], 'radicals_racm2', $racm2_kpp, $racm2_mecca);
my ($racm2_sorted_plot_data, $racm2_legend) = sort_data_for_plot($production_rates{'radicals_racm2'}, $consumption_rates{'radicals_racm2'});
my $racm2_plot_title = "(f) RACM2";

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA_TIM->new($cbm4_run);
my $cbm4_eqnfile = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqnfile);
my $cbm4_ro2file = "/work/users/jco/MECCA/CBM4_tagging/RO2_species.txt";
my @cbm4_no2_reservoirs = get_no2_reservoirs($cbm4_kpp, $cbm4_ro2file);
my $cbm4_radical_file = "/work/users/jco/MECCA/CBM4_tagging/radicals.txt";
my @cbm4_radicals = get_radicals($cbm4_radical_file);
$families{'radicals_cbm4'} = [ @cbm4_radicals, @cbm4_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_cbm4'}, $consumption_rates{'radicals_cbm4'}) = get_rates($species[0], 'radicals_cbm4', $cbm4_kpp, $cbm4_mecca);
my ($cbm4_sorted_plot_data, $cbm4_legend) = sort_data_for_plot($production_rates{'radicals_cbm4'}, $consumption_rates{'radicals_cbm4'});
my $cbm4_plot_title = "(h) CBM-IV";

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA_TIM->new($cb05_run);
my $cb05_eqnfile = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqnfile);
my $cb05_ro2file = "/work/users/jco/MECCA/CB05_tagging/RO2_species.txt";
my @cb05_no2_reservoirs = get_no2_reservoirs($cb05_kpp, $cb05_ro2file);
my $cb05_radical_file = "/work/users/jco/MECCA/CB05_tagging/radicals.txt";
my @cb05_radicals = get_radicals($cb05_radical_file);
$families{'radicals_cb05'} = [ @cb05_radicals, @cb05_no2_reservoirs, 'HO2NO2' ];

($production_rates{'radicals_cb05'}, $consumption_rates{'radicals_cb05'}) = get_rates($species[0], 'radicals_cb05', $cb05_kpp, $cb05_mecca);
my ($cb05_sorted_plot_data, $cb05_legend) = sort_data_for_plot($production_rates{'radicals_cb05'}, $consumption_rates{'radicals_cb05'});
my $cb05_plot_title = "(i) CB05";

#Create x-axis for plot in hours -> time axis is the same for all model runs
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list; 

my ($plot) = plot({#create plot
        y_max           => 14,
        y_min           => -16,
        breaks          => 2,
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

sub get_cbm4_allocations { #each parent VOC expressed as CBM-IV allocation
    my ($parent) = @_;
    my @allocations = ();
    
    if ($parent eq "C2H6" or $parent eq "C3H8" or $parent eq "NC4H10" or $parent eq "IC4H10" or $parent eq "NC5H12" or $parent eq "IC5H12" or $parent eq "NC6H14" or $parent eq "NC7H16" or $parent eq "NC8H18" or $parent eq "BENZENE"){
        @allocations = qw( PAR );
    } elsif ($parent eq "CH4" ){
        @allocations = qw( CH4 );
    } elsif ($parent eq "C3H6" or $parent eq "BUT1ENE") {
        @allocations = qw( OLE PAR );
    } elsif ($parent eq "MEPROPENE") {
        @allocations = qw( ALD2 );
    } elsif ($parent eq "C2H4") {
        @allocations = qw ( ETH );
    } elsif ($parent eq "C5H8") {
        @allocations = qw( ISOP );
    } elsif ($parent eq "TOLUENE" ) {
        @allocations = qw( TOL ) ;
    } elsif ($parent eq "MXYL" or $parent eq "OXYL" or $parent eq "PXYL") {
        @allocations = qw( XYL );
    } elsif ($parent eq "EBENZ") {
        @allocations = qw( TOL PAR );
    }
    return \@allocations;
}

sub get_cb05_allocations { #each parent VOC expressed as CB05 allocation
    my ($parent) = @_;
    my @allocations = ();
    
    if ($parent eq "C3H8" or $parent eq "NC4H10" or $parent eq "IC4H10" or $parent eq "NC5H12" or $parent eq "IC5H12" or $parent eq "NC6H14" or $parent eq "NC7H16" or $parent eq "NC8H18" or $parent eq "BENZENE"){
        @allocations = qw( PAR );
    } elsif ($parent eq "CH4" ){
        @allocations = qw( CH4 );
    } elsif ($parent eq "C3H6" or $parent eq "BUT1ENE") {
        @allocations = qw( OLE PAR );
    } elsif ($parent eq "MEPROPENE") {
        @allocations = qw( FORM PAR );
    } elsif ($parent eq "C2H4") {
        @allocations = qw ( ETH );
    } elsif ($parent eq "C5H8") {
        @allocations = qw( ISOP );
    } elsif ($parent eq "TOLUENE" ) {
        @allocations = qw( TOL ) ;
    } elsif ($parent eq "MXYL" or $parent eq "OXYL" or $parent eq "PXYL") {
        @allocations = qw( XYL );
    } elsif ($parent eq "EBENZ") {
        @allocations = qw( TOL PAR );
    } elsif ($parent eq "C2H6") {
        @allocations = qw( ETHA );
    }
    return \@allocations;
}

sub get_rates {
    my ($parent_species, $species, $kpp, $mecca) = @_;
    
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
    }
   
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq $parent_species);
		my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $string = $kpp->reaction_string($reaction);
        $string =~ s/_$parent//g; #remove tag from string
        my ($reactants, $products) = split / = /, $string;
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
        $string =~ s/_$parent//g; #remove tag from string
        my ($reactants, $products) = split / = /, $string;
        $species_consumption_rates{$reactants} += $rate(1:$ntime-2); 
    } 

    #get parent species emissions for each mechanism
    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($species =~ /cbm4/) {
        my $primary = get_cbm4_allocations($parent_species);
        foreach my $VOC (@$primary) {
            my $name = "${VOC}_${parent_species}";
            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
            $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
        }
    } elsif ($species =~ /cb05/){
        my $primary = get_cb05_allocations($parent_species);
        foreach my $VOC (@$primary) {
            my $name = "${VOC}_${parent_species}";
            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
            $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
        }
    } else {
        my $parent_source = $mecca->balance($parent_species); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
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
    my $prod_others_max = 6e-5;
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

sub plot { #create dataframe and then create plot
    my ($args) = @_;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(Cairo) `);
    $R->run(q` library(scales) `);

    $R->set('y.max', $args->{y_max});
    $R->set('y.min', $args->{y_min});
    $R->set('y.breaks', $args->{breaks});
    $R->set('time', [@{$args->{times}}]);
    $R->run(q` my.colours = c( "Production Others" = "#696537", "HCHO + hv" = "#f9c600", "CH2O + hv" = "#f9c600", "OH + PPN" = "#76afca", "C3PAN1 + OH" = "#dc3522", "PAN" = "#8c6238", "CARB13 + hv" = "#8b1537", "CARB11A + hv" = "#9bb08f", "NO + RN16O2" = "#e7e85e", "ALKO2 + NO" = "#e7e85e", "HC5P + NO" = "#e7e85e", "HC5P + HO2" = "#0352cb", "HO2 + RN16O2" = "#0352cb", "ALKO2 + HO2" = "#0352cb", "CH3O2 + HO2" = "#86b650", "OH + RN16NO3" = "#6c254f", "HO2 + RN15AO2" = "#ee6738", "CH3CO3 + HO2" = "#58691b", "HC5 + OH" = "#8ed6d5", "KET + hv" = "#86b650", "MGLY + hv" = "#c65d6c", "NO + XO2" = "#888a87", "PPN" = "#0e5c28", "MEK + hv" = "#b569b3", "ROR" = "#0e5c28", "C2O3 + NO" = "#f8c56c", "CXO3 + NO" = "#2c9def", "OH + PAR" = "#ae4903", "Consumption Others" = "#6c254f" ) `,
            q` my.names = c("PAN" = "PAN Deposition", "PPN" = "PPN Deposition") `,
    );

    #general plot R function
    $R->run(q` mech.plot = function(data, legend.title, plot.title, legend) { sums = colSums(data[,-1]); rep.no = nrow(data); regex = as.numeric(gsub("^[A-Z]+ + [A-Z]+", "", sums, perl = TRUE) ); type = c(); for ( i in 1:length(sums) ) { type = c(type, rep(regex[i], rep.no)) }; plot.data = melt(data = data, id = names(data)[1], measured = names(data)[-1] ); colnames(plot.data) = c("time", "reaction", "rate"); plot.data$subset.condition = type; reaction.levels = (levels(factor(plot.data$reaction))); plot.data$reaction = ordered(plot.data$reaction, levels = reaction.levels); plot.data = ddply( plot.data, .(reaction)); plot = ggplot(data = plot.data, aes(x = time, y = rate, fill = reaction)); plot = plot + geom_area(data = subset(plot.data, subset.condition > 0), position = "stack", alpha = 1); plot = plot + geom_area(data = subset(plot.data, subset.condition > 0), position = "stack", colour = "black", show_guide = FALSE); plot = plot + geom_area(data = subset(plot.data, subset.condition < 0), position = "stack", alpha = 1); plot = plot + geom_area(data = subset(plot.data, subset.condition < 0), position = "stack", colour = "black", show_guide = FALSE); plot = plot + guides(fill = guide_legend(title = legend.title)); plot = plot + ggtitle(plot.title); plot = plot + scale_x_continuous(limits=c(0, 7), breaks=seq(0, 7, 1)); plot = plot + scale_y_continuous(limits=c(y.min, y.max), breaks=seq(y.min, y.max, y.breaks)) ; plot = plot + scale_fill_manual( name = "reaction", limits = legend, values = my.colours, labels = my.names); plot = plot + theme_bw() ; plot = plot + theme(legend.key.size = unit(5.5, "cm")) ; plot = plot +  theme(axis.text.x = element_text(size = 170)) ; plot = plot + theme(axis.text.y = element_text(size = 170)) ; plot = plot + theme(axis.title.x = element_blank()) ; plot = plot + theme(legend.text = element_text(size = 110), legend.title = element_blank()) ; plot = plot + theme(legend.key = element_blank()) ; plot = plot + theme(axis.title.y = element_blank()) ; plot = plot + theme(plot.title = element_text(size = 230, face = "bold", vjust = 1)) ; plot = plot + theme(legend.justification = c(0.95, 0.05), legend.position = c(0.95, 0.05)) ; return(plot) } `);
 
    #MCM v3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.plot.title', $args->{mcm3_1_title});
    $R->set('mcm3.1.legend', $args->{mcm3_1_legend});
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate*1000000`); 
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
             $R->run(q` mcm3.2.data[name] = rate*1000000`); 
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
             $R->run(q` cri.data[name] = rate*1000000`); 
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
             $R->run(q` mozart.data[name] = rate*1000000`); 
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
             $R->run(q` radm2.data[name] = rate*1000000`); 
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
             $R->run(q` racm.data[name] = rate*1000000`); 
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
             $R->run(q` racm2.data[name] = rate*1000000`); 
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
            $R->run(q` cbm4.data[name] = rate*1000000`); 
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
            $R->run(q` cb05.data[name] = rate*1000000`); 
        }
    } 
    $R->run(q` cb05.plot = mech.plot(cb05.data, "CB05", cb05.plot.title, cb05.legend) `);

    $R->run(q` CairoPDF(file = "pentane_radical_budgets.pdf", width = 190, height = 200) `, 
            q` y.label = textGrob(expression(bold(paste("\nMolecules (Product) ", s^-1, "/Molecules (NMVOC) x ", 10^6))), rot = 90, gp = gpar(fontsize = 210), vjust = 0.6)`,
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
                                       nrow = 1, 
                                       sub = textGrob("\nTime (days)\n", gp = gpar(fontsize = 210, fontface = "bold"), vjust = 0.2), 
                                       widths=unit.c(unit(20, "lines"), unit(1, "npc") - unit(20, "lines"))) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );
    $R->stop(); 
}
