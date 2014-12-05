#! /usr/bin/perl
# Analyse day-time pan budget during toluene breakdown, compare for all mechanisms - stacked bar plot
# Version 0: Jane Coates 28/07/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my (%families, %weights, %production_rates, %consumption_rates);
my @species = qw( TOLUENE TOL );

#mcm 3.2
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_radical_file = "/work/users/jco/MECCA/MCM_3.2_tagged/radicals.txt";
my @mcm_3_2_radicals = get_pans($mcm_3_2_radical_file);
my $mcm_3_2_panfile = "/work/users/jco/Mechanisms/MCM/MCM-3.2/All_pans.txt";
my @mcm_3_2_pans = get_pans($mcm_3_2_panfile); 
my @mcm_3_2_pan_producers = pan_producers($mcm_3_2_kpp, @mcm_3_2_pans);
$families{'PAN_mcm_3_2'} = [ @mcm_3_2_pans , @mcm_3_2_pan_producers ]; 
my $ntime = $mcm_3_2_mecca->time->nelem; #number of time points 
($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'}) = get_rates($species[0], 'PAN_mcm_3_2', $mcm_3_2_kpp, $mcm_3_2_mecca, \@mcm_3_2_radicals, \@mcm_3_2_pans);
my ($mcm_3_2_sorted_plot_data, $mcm_3_2_legend) = sort_data_for_plot($production_rates{'PAN_mcm_3_2'}, $consumption_rates{'PAN_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#mcm 3.1
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA->new($mcm_3_1_run); 
my $mcm_3_1_eqnfile = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.eqn";
my $mcm_3_1_kpp = KPP->new($mcm_3_1_eqnfile); 
my $mcm_3_1_radical_file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/radicals.txt";
my @mcm_3_1_radicals = get_pans($mcm_3_1_radical_file); 
my $mcm_3_1_panfile = "/work/users/jco/Mechanisms/MCM/MCM-3.1/All_pans.txt";
my @mcm_3_1_pans = get_pans($mcm_3_1_panfile); 
my @mcm_3_1_pan_producers = pan_producers($mcm_3_1_kpp, @mcm_3_1_pans);
$families{'PAN_mcm_3_1'} = [ @mcm_3_1_pans , @mcm_3_1_pan_producers ]; 
($production_rates{'PAN_mcm_3_1'}, $consumption_rates{'PAN_mcm_3_1'}) = get_rates($species[0], 'PAN_mcm_3_1', $mcm_3_1_kpp, $mcm_3_1_mecca, \@mcm_3_1_radicals, \@mcm_3_1_pans);
my ($mcm_3_1_sorted_plot_data, $mcm_3_1_legend) = sort_data_for_plot($production_rates{'PAN_mcm_3_1'}, $consumption_rates{'PAN_mcm_3_1'});
my $mcm_3_1_plot_title = "(b) MCM v3.1";

#cri
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run);
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile);
my $cri_radical_file = "/work/users/jco/MECCA/CRI_tagging/radicals.txt";
my @cri_radicals = get_pans($cri_radical_file); 
my $cri_panfile = "/work/users/jco/Mechanisms/CRI/CRI_v2_full/All_pans.txt";
my @cri_pans = get_pans($cri_panfile); 
my @cri_pan_producers = pan_producers($cri_kpp, @cri_pans);
$families{'PAN_cri'} = [ @cri_pans , @cri_pan_producers ]; 
($production_rates{'PAN_cri'}, $consumption_rates{'PAN_cri'}) = get_rates($species[0], 'PAN_cri', $cri_kpp, $cri_mecca, \@cri_radicals, \@cri_pans); 
my ($cri_sorted_plot_data, $cri_legend) = sort_data_for_plot($production_rates{'PAN_cri'}, $consumption_rates{'PAN_cri'});
my $cri_plot_title = "(c) CRI v2";

#mozart
my $mozart_run = '/work/users/jco/MECCA/MOZART_tagging/boxmodel';
my $mozart_mecca = MECCA->new($mozart_run);
my $mozart_eqnfile = '/work/users/jco/MECCA/MOZART_tagging/gas.eqn';
my $mozart_kpp = KPP->new($mozart_eqnfile);
my $mozart_radical_file = "/work/users/jco/MECCA/MOZART_tagging/radicals.txt";
my @mozart_radicals = get_pans($mozart_radical_file); 
my $mozart_panfile = "/work/users/jco/Mechanisms/MOZART/MOZART/All_pans.txt";
my @mozart_pans = get_pans($mozart_panfile); 
my @mozart_pan_producers = pan_producers($mozart_kpp, @mozart_pans);
$families{'PAN_mozart'} = [ @mozart_pans , @mozart_pan_producers ]; 
($production_rates{'PAN_mozart'}, $consumption_rates{'PAN_mozart'}) = get_rates($species[0], 'PAN_mozart', $mozart_kpp, $mozart_mecca, \@mozart_radicals, \@mozart_pans);
my ($mozart_sorted_plot_data, $mozart_legend) = sort_data_for_plot($production_rates{'PAN_mozart'}, $consumption_rates{'PAN_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#radm2
my $radm2_run = '/work/users/jco/MECCA/RADM2_tagged/boxmodel';
my $radm2_mecca = MECCA->new($radm2_run);
my $radm2_eqnfile = '/work/users/jco/MECCA/RADM2_tagged/gas.eqn';
my $radm2_kpp = KPP->new($radm2_eqnfile);
my $radm2_radical_file = "/work/users/jco/MECCA/RADM2_tagged/radicals.txt";
my @radm2_radicals = get_pans($radm2_radical_file); 
my $radm2_panfile = "/work/users/jco/Mechanisms/RADM2/All_pans.txt";
my @radm2_pans = get_pans($radm2_panfile); 
my @radm2_pan_producers = pan_producers($radm2_kpp, @radm2_pans);
$families{'PAN_radm2'} = [ @radm2_pans , @radm2_pan_producers ]; 
($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'}) = get_rates($species[1], 'PAN_radm2', $radm2_kpp, $radm2_mecca, \@radm2_radicals, \@radm2_pans);
my ($radm2_sorted_plot_data, $radm2_legend) = sort_data_for_plot($production_rates{'PAN_radm2'}, $consumption_rates{'PAN_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm
my $racm_run = '/work/users/jco/MECCA/RACM_tagging/boxmodel';
my $racm_mecca = MECCA->new($racm_run);
my $racm_eqnfile = '/work/users/jco/MECCA/RACM_tagging/gas.eqn';
my $racm_kpp = KPP->new($racm_eqnfile);
my $racm_radical_file = "/work/users/jco/MECCA/RACM_tagging/radicals.txt";
my @racm_radicals = get_pans($racm_radical_file); 
my $racm_panfile = "/work/users/jco/Mechanisms/RACM/All_pans.txt";
my @racm_pans = get_pans($racm_panfile); 
my @racm_pan_producers = pan_producers($racm_kpp, @racm_pans);
$families{'PAN_racm'} = [ @racm_pans , @racm_pan_producers ]; 
($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'}) = get_rates($species[1], 'PAN_racm', $racm_kpp, $racm_mecca, \@racm_radicals, \@racm_pans);
my ($racm_sorted_plot_data, $racm_legend) = sort_data_for_plot($production_rates{'PAN_racm'}, $consumption_rates{'PAN_racm'});
my $racm_plot_title = "(e) RACM";

#racm2
my $racm2_run = '/work/users/jco/MECCA/RACM2_tagged/boxmodel';
my $racm2_mecca = MECCA->new($racm2_run);
my $racm2_eqnfile = '/work/users/jco/MECCA/RACM2_tagged/gas.eqn';
my $racm2_kpp = KPP->new($racm2_eqnfile);
my $racm2_radical_file = "/work/users/jco/MECCA/RACM2_tagged/radicals.txt";
my @racm2_radicals = get_pans($racm2_radical_file); 
my $racm2_panfile = "/work/users/jco/Mechanisms/RACM2/All_pans.txt";
my @racm2_pans = get_pans($racm2_panfile); 
my @racm2_pan_producers = pan_producers($racm2_kpp, @racm2_pans);
$families{'PAN_racm2'} = [ @racm2_pans , @racm2_pan_producers ]; 
($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'}) = get_rates($species[1], 'PAN_racm2', $racm2_kpp, $racm2_mecca, \@racm2_radicals, \@racm2_pans);
my ($racm2_sorted_plot_data, $racm2_legend) = sort_data_for_plot($production_rates{'PAN_racm2'}, $consumption_rates{'PAN_racm2'});
my $racm2_plot_title = "(f) RACM2";

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA->new($cbm4_run);
my $cbm4_eqnfile = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqnfile);
my $cbm4_radical_file = "/work/users/jco/MECCA/CBM4_tagging/radicals.txt";
my @cbm4_radicals = get_pans($cbm4_radical_file); 
my $cbm4_panfile = "/work/users/jco/Mechanisms/CBM-IV/All_pans.txt";
my @cbm4_pans = get_pans($cbm4_panfile); 
my @cbm4_pan_producers = pan_producers($cbm4_kpp, @cbm4_pans);
$families{'PAN_cbm4'} = [ @cbm4_pans , @cbm4_pan_producers ]; 
($production_rates{'PAN_cbm4'}, $consumption_rates{'PAN_cbm4'}) = get_rates($species[0], 'PAN_cbm4', $cbm4_kpp, $cbm4_mecca, \@cbm4_radicals, \@cbm4_pans);
my ($cbm4_sorted_plot_data, $cbm4_legend) = sort_data_for_plot($production_rates{'PAN_cbm4'}, $consumption_rates{'PAN_cbm4'});
my $cbm4_plot_title = "(h) CBM-IV";

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA->new($cb05_run);
my $cb05_eqnfile = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqnfile);
my $cb05_radical_file = "/work/users/jco/MECCA/CB05_tagging/radicals.txt";
my @cb05_radicals = get_pans($cb05_radical_file); 
my $cb05_panfile = "/work/users/jco/Mechanisms/CB05/All_pans.txt";
my @cb05_pans = get_pans($cb05_panfile); 
my @cb05_pan_producers = pan_producers($cb05_kpp, @cb05_pans);
$families{'PAN_cb05'} = [ @cb05_pans , @cb05_pan_producers ]; 
($production_rates{'PAN_cb05'}, $consumption_rates{'PAN_cb05'}) = get_rates($species[0], 'PAN_cb05', $cb05_kpp, $cb05_mecca, \@cb05_radicals, \@cb05_pans);
my ($cb05_sorted_plot_data, $cb05_legend) = sort_data_for_plot($production_rates{'PAN_cb05'}, $consumption_rates{'PAN_cb05'});
my $cb05_plot_title = "(i) CB05";

#Create x-axis for plot in hours -> time axis is the same for all model runs
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 3400;
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

my ($plot) = plot({#create plot
        times           => \@time_blocks,
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

sub get_pans { #get radicals from file
    my ($file) = @_; 

    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @lines = <$in>;
    close $in;
    my %hash;
    my @separate = map { split /\s+/, $_} @lines;
    $hash{$_} += 1 foreach (@separate);
    return keys %hash;
} 

sub get_rates {
    my ($parent_species, $species, $kpp, $mecca, $radicals, $pans) = @_;
    
    print "$species\n";
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
        my ($reactants, $products) = split / = /, $string;
        $reactants =~ s/^\s+|\s+$//g;
        $reactants =~ s/_${parent}//g; 
        my ($category) = get_category($reactants, $radicals, $pans, 'production');
        $species_production_rates{$category} += $rate(1:$ntime-2);
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
        $reactants =~ s/^\s+|\s+$//g;
        $reactants =~ s/_${parent}//g; 
        my ($category) = get_category($reactants, $radicals, $pans, 'consumption');
        $species_consumption_rates{$category} += $rate(1:$ntime-2); 
    } 

    #get parent species emissions for each mechanism
    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($species =~ /cbm4/ or $species =~ /cb05/) {
        my $name = "TOL_TOLUENE";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt ;
    } else {
        my $parent_source = $mecca->balance($parent_species); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $species_production_rates{$_} /= $parent_emissions foreach (sort keys %species_production_rates);
    $species_consumption_rates{$_} /= $parent_emissions foreach (sort keys %species_consumption_rates);
    return (\%species_production_rates, \%species_consumption_rates);
}

sub get_category {
    my ($reactants, $radicals, $pans, $flag) = @_;
    my $category;
    #print "$reactants\n";
    if ($reactants =~ /hv/) {
        $category = 'Photolysis' ;
    } elsif ($reactants ~~ @$pans) {
        $category = 'Deposition';
    } elsif ($reactants ~~ @$radicals and $flag eq 'consumption') {
        $category = 'Radical Consumption';
    } elsif ($reactants ~~ @$radicals and $flag eq 'production') {
        $category = 'Radical Production';
    } elsif ($reactants =~ /(.*?)\+(.*?)/) {
        my ($one, $two) = split / \+ /, $reactants;
        if ($one ~~ @$pans and $two eq 'OH') {
            $category = 'PAN + OH Reactions';
        } elsif ($one eq 'OH' and $two ~~ @$pans) {
            $category = 'PAN + OH Reactions';
        } elsif ($one eq 'OH' or $two eq 'OH') {
            $category = 'Aldehyde + OH Reactions';
        } elsif ($one eq 'NO3' or $two eq 'NO3') {
            if ($flag eq 'production') {
                $category = 'Production: Nitrate Reaction';
            } elsif ($flag eq 'consumption') {
                $category = 'Consumption: Nitrate Reaction';
            }
        } elsif ($one eq 'NO' or $two eq 'NO') {
            if ($flag eq 'production') {
                $category = 'Production: NO Reaction';
            } elsif ($flag eq 'consumption') {
                $category = 'Consumption: NO Reaction';
            }
        } elsif ($one eq 'HO2' or $two eq 'HO2') {
            $category = 'HO2 Reaction';
        } elsif ($one eq 'O3' or $two eq 'O3') {
            $category = 'Ozone Reaction';
        } elsif ($one eq 'O' or $two eq 'O') {
            $category = 'Atomic Oxygen Reaction';
        } else {
            $category = 'other';
            print "$reactants\n";
        }
    } else {
        $category = 'other';
        print "$reactants\n";
    }
    return $category;
}

sub sort_data_for_plot { #create hash with production of the reactions
    my ($production_rates, $consumption_rates) = @_;
    my @sorted_plot_data; 

    my $sort_function = sub { $_[0]->sum };
    my @sorted_cons = reverse sort { &$sort_function($consumption_rates->{$b}) <=> &$sort_function($consumption_rates->{$a}) } keys %$consumption_rates; 
    push @sorted_plot_data, { $_ => $consumption_rates->{$_} } foreach (@sorted_cons); #sum up rates of reactions, starting with reaction with lowest sum, consumption others added separately 

    my @sorted_prod = sort { &$sort_function($production_rates->{$b}) <=> &$sort_function($production_rates->{$a}) } keys %$production_rates; 
    push @sorted_plot_data, { $_ => $production_rates->{$_} } foreach (@sorted_prod);#sum up rates of reactions, starting with reaction with lowest sum, production others added separately 

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
    $R->run(q` library(gtable) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(Cairo) `);
    $R->run(q` library(scales) `);

    $R->set('time', [@{$args->{times}}]);
    $R->run(q` my.colours = c( "Photolysis" = "#76afca", "Radical Production" = "#ef6638", "Deposition" = "#898989", "PAN + OH Reactions" = "#b569b3", "Aldehyde + OH Reactions" = "#8c6238", "Radical Consumption" = "#f9c500", "Consumption: NO Reaction" = "#8d1435" , "Production: NO Reaction" = "#dc3522", "HO2 Reaction" = "#ba8b01", "Ozone Reaction" = "#0c3f76", "Atomic Oxygen Reaction" = "#4b9483" ) `,
            q` my.legend = c( "Ozone Reaction", "Production: NO Reaction", "Photolysis", "Radical Production", "Aldehyde + OH Reactions", "Consumption: NO Reaction", "Deposition", "Radical Consumption", "PAN + OH Reactions", "HO2 Reaction" ) `,
    );

    #general plot R function
    $R->run(q` mech.plot = function(data, plot.title, plot.legend) {data = ddply(data, .(time), colwise(sum)) ;
                                                                    subset = data[1:7,] ;
                                                                    plot.data = melt( data = subset, id = names(subset)[1], measured = names(subset)[-1] ); 
                                                                    colnames(plot.data) = c("time", "process", "rate"); 
                                                                    process.levels = (levels(factor(plot.data$process))); 
                                                                    plot.data$process = ordered(plot.data$process, levels = process.levels); 
                                                                    plot.data = ddply( plot.data, .(process)); 
                                                                    plot = ggplot(data = plot.data, aes(x = time, y = rate)); 
                                                                    plot = plot + geom_bar(data = subset(plot.data, rate > 0), aes(fill = process), stat = "identity", width = 0.5) ;
                                                                    plot = plot + geom_bar(data = subset(plot.data, rate < 0), aes(fill = process), stat = "identity", width = 0.5) ;
                                                                    plot = plot + ggtitle(plot.title); 
                                                                    plot = plot + scale_fill_manual( limits = my.legend, values = my.colours, guide = guide_legend(label.theme = element_text(size = 160, angle = 0, face = "bold"), label.hjust = 0.0, label.vjust = 0.5, nrow = 2)); 
                                                                    plot = plot + theme_bw() ; 
                                                                    plot = plot + theme(axis.text.x = element_text(size = 140)) ; 
                                                                    plot = plot + theme(axis.text.y = element_text(size = 150)) ; 
                                                                    plot = plot + theme(axis.title.x = element_blank()) ; 
                                                                    plot = plot + theme(axis.title.y = element_blank()) ; 
                                                                    plot = plot + theme(plot.title = element_text(size = 200, face = "bold", vjust = 1)) ; 
                                                                    plot = plot + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) ;
                                                                    plot = plot + scale_y_continuous(limits = c(-50, 60), breaks = seq(-50, 60, 20)) ;
                                                                    plot = plot + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) ;
                                                                    plot = plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.height = unit(20, "cm"), legend.key.width = unit(20, "cm"), legend.key = element_blank()); 
                                                                    return(plot) } `);
 
    #MCM v3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.plot.title', $args->{mcm3_1_title});
    $R->set('mcm3.1.legend', $args->{mcm3_1_legend});
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate*100000`); 
        }
    } 
    $R->run(q` mcm3.1.plot = mech.plot(mcm3.1.data, mcm3.1.plot.title, mcm3.1.legend) `);

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
 
     #CRI v2
     $R->run(q` cri.data = data.frame(time)`);
     $R->set('cri.plot.title', $args->{cri_title});
     $R->set('cri.legend', $args->{cri_legend});
     foreach my $ref (@{$args->{cri_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` cri.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` cri.plot = mech.plot(cri.data, cri.plot.title, cri.legend) `);
 
     #MOZART-4
     $R->run(q` mozart.data = data.frame(time)`);
     $R->set('mozart.plot.title', $args->{mozart_title});
     $R->set('mozart.legend', $args->{mozart_legend});
     foreach my $ref (@{$args->{mozart_data}}) {
         for my $key (keys %$ref) {
             $R->set('name', $key);
             $R->set('rate', [@{$ref->{$key}}]);
             $R->run(q` mozart.data[name] = rate*100000`); 
         }
     } 
     $R->run(q` mozart.plot = mech.plot(mozart.data, mozart.plot.title, mozart.legend) `);
 
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
 
    #CBM-IV
    $R->set('cbm4.plot.title', $args->{cbm4_title});
    $R->set('cbm4.legend', $args->{cbm4_legend});
    $R->run(q` cbm4.data = data.frame(time) `);
    foreach my $ref (@{$args->{cbm4_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cbm4.data[name] = rate*100000`); 
        }
    } 
     $R->run(q` cbm4.plot = mech.plot(cbm4.data, cbm4.plot.title, cbm4.legend) `);

    #CB05
    $R->set('cb05.plot.title', $args->{cb05_title});
    $R->set('cb05.legend', $args->{cb05_legend});
    $R->run(q` cb05.data = data.frame(time) `);
    foreach my $ref (@{$args->{cb05_data}}) {
        for my $key (keys %$ref) {
            $R->set('name', $key);
            $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cb05.data[name] = rate*100000`); 
        }
    } 
    $R->run(q` cb05.plot = mech.plot(cb05.data, cb05.plot.title, cb05.legend) `);

     $R->run(q` CairoPDF(file = "toluene_pan_budgets.pdf", width = 210, height = 200) `, 
            q` y.label = textGrob(expression(bold(paste("\nMolecules (Product) ", s^-1, "/Molecules (NMVOC) x ", 10^5))), rot = 90, gp = gpar(fontsize = 210), vjust = 0.6)`,
            q` legend = gtable_filter(ggplot_gtable(ggplot_build(mcm3.2.plot)), "guide-box") `,
            q` main.plot = grid.arrange(y.label, 
                                        arrangeGrob(mcm3.2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none"), 
                                                    mcm3.1.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    cri.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    radm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none"), 
                                                    racm.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    racm2.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    mozart.plot + theme(legend.position = "none"), 
                                                    cbm4.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    cb05.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"), 
                                                    nrow = 3), 
                                       blank = grid.rect(gp = gpar(col = "white")), legend,
                                       nrow = 2, ncol = 2,
                                       sub = textGrob("\n", gp = gpar(fontsize = 210, fontface = "bold"), vjust = 0.2), 
                                       widths=unit.c(unit(25, "lines"), unit(1, "npc") - unit(25, "lines")),
                                       heights = c(180, 30)) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );
    $R->stop(); 
}
