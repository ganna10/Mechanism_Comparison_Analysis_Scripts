#! /usr/bin/env perl
# Analyse number of Carbon atoms per radical produced during Pentane breakdown, compare for all mechanisms - stacked bar plot
# Version 0: Jane Coates 13/02/2014
# Version 1: Jane Coates 08/05/2014 changed plots to stacked bar charts and corrected CBMs emissions
# Version 2: Jane Coates 03/06/2014 removing night-time production as it is insignificant compared to day-time production
# Version 3: Jane Coates 19/06/2014 XO2 included as its own category
# Version 4: Jane Coates 8/7/2014 Attributing HO2 and XO2 to carbon numbers of their reactants

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my (%n_carbon, %families, %weights, %production_rates);
$families{'HO2x'} = [ qw( HO2 HO2NO2 ) ];

my $species = "NC5H12";
my $moz_species = "BIGALK";
my $RACM_species = "HC5";

#mcm 3.1
my $mcm_3_1_species_file = "/work/users/jco/Mechanisms/MCM/MCM-3.1/species_smiles.txt";
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA->new($mcm_3_1_run); 
my $mcm_3_1_eqnfile = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.eqn";
my $mcm_3_1_kpp = KPP->new($mcm_3_1_eqnfile); 
my $mcm_3_1_ro2file = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/RO2_species.txt";
my @mcm_3_1_no2_reservoirs = get_no2_reservoirs($mcm_3_1_kpp, $mcm_3_1_ro2file);
$families{'Ox_mcm_3_1'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @mcm_3_1_no2_reservoirs ];
$weights{'Ox_mcm_3_1'} = { NO3 => 2, N2O5 => 3};

my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points

($n_carbon{'MCM 3.1'}) = mcm_n_carbon($mcm_3_1_species_file); 
($production_rates{'Ox_mcm_3_1'}) = get_production_rates($species, 'Ox_mcm_3_1', $n_carbon{'MCM 3.1'}, $mcm_3_1_kpp, $mcm_3_1_mecca);
my ($mcm_3_1_plot_data) = get_plot_hash('Ox_mcm_3_1', $production_rates{'Ox_mcm_3_1'});
my $mcm_3_1_plot_title = "(b) MCM v3.1";

#mcm 3.2
my $mcm_3_2_species_file = "/work/users/jco/Mechanisms/MCM/MCM-3.2/smiles.out";
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $mcm_3_2_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $mcm_3_2_kpp = KPP->new($mcm_3_2_eqnfile); 
my $mcm_3_2_ro2file = "/work/users/jco/MECCA/MCM_3.2_tagged/RO2_species.txt";
my @mcm_3_2_no2_reservoirs = get_no2_reservoirs($mcm_3_2_kpp, $mcm_3_2_ro2file);
$families{'Ox_mcm_3_2'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @mcm_3_2_no2_reservoirs ];
$weights{'Ox_mcm_3_2'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'MCM 3.2'}) = mcm_n_carbon($mcm_3_2_species_file); 
($production_rates{'Ox_mcm_3_2'}) = get_production_rates($species, 'Ox_mcm_3_2', $n_carbon{'MCM 3.2'}, $mcm_3_2_kpp, $mcm_3_2_mecca);
my ($mcm_3_2_plot_data) = get_plot_hash('Ox_mcm_3_2', $production_rates{'Ox_mcm_3_2'});
my $mcm_3_2_plot_title = "(a) MCM v3.2";

#mozart
my $mozart_species_file = "/work/users/jco/Mechanisms/MOZART/MOZART/chem_mech.in";
my $mozart_run = '/work/users/jco/MECCA/MOZART_tagging/boxmodel';
my $mozart_mecca = MECCA->new($mozart_run);
my $mozart_eqnfile = '/work/users/jco/MECCA/MOZART_tagging/gas.eqn';
my $mozart_kpp = KPP->new($mozart_eqnfile);
my $mozart_ro2file = "/work/users/jco/MECCA/MOZART_tagging/RO2_species.txt";
my @mozart_no2_reservoirs = get_no2_reservoirs($mozart_kpp, $mozart_ro2file);
$families{'Ox_mozart'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @mozart_no2_reservoirs ];
$weights{'Ox_mozart'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'MOZART'}) = mozart_n_carbon($mozart_species_file);
($production_rates{'Ox_mozart'}) = get_production_rates($moz_species, 'Ox_mozart', $n_carbon{'MOZART'}, $mozart_kpp, $mozart_mecca);
my ($mozart_plot_data) = get_plot_hash('Ox_mozart', $production_rates{'Ox_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#cri
my $cri_species_file = "/work/users/jco/Mechanisms/CRI/CRI_v2_full/carbons.txt";
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA->new($cri_run);
my $cri_eqnfile = "/work/users/jco/MECCA/CRI_tagging/gas.eqn";
my $cri_kpp = KPP->new($cri_eqnfile);
my $cri_ro2file = "/work/users/jco/MECCA/CRI_tagging/RO2_species.txt";
my @cri_no2_reservoirs = get_no2_reservoirs($cri_kpp, $cri_ro2file);
$families{'Ox_cri'} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5), @cri_no2_reservoirs ];
$weights{'Ox_cri'} = { NO3 => 2, N2O5 => 3 };

($n_carbon{'CRI'}) = carbons_others($cri_species_file);
($production_rates{'Ox_cri'}) = get_production_rates($species, 'Ox_cri', $n_carbon{'CRI'}, $cri_kpp, $cri_mecca);
my ($cri_plot_data) = get_plot_hash('Ox_cri', $production_rates{'Ox_cri'});
my $cri_plot_title = "(c) CRI v2";

#radm2
my $radm2_species_file = "/work/users/jco/Mechanisms/RADM2/carbon_numbers.txt";
my $radm2_run = '/work/users/jco/MECCA/RADM2_tagged/boxmodel';
my $radm2_mecca = MECCA->new($radm2_run);
my $radm2_eqnfile = '/work/users/jco/MECCA/RADM2_tagged/gas.eqn';
my $radm2_kpp = KPP->new($radm2_eqnfile);
my $radm2_ro2file = "/work/users/jco/MECCA/RADM2_tagged/RO2_species.txt";
my @radm2_no2_reservoirs = get_no2_reservoirs($radm2_kpp, $radm2_ro2file);
$families{'Ox_radm2'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @radm2_no2_reservoirs ];
$weights{'Ox_radm2'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'RADM2'}) = carbons_others($radm2_species_file);
($production_rates{'Ox_radm2'}) = get_production_rates($RACM_species, 'Ox_radm2', $n_carbon{'RADM2'}, $radm2_kpp, $radm2_mecca);
my ($radm2_plot_data) = get_plot_hash('Ox_radm2', $production_rates{'Ox_radm2'});
my $radm2_plot_title = "(d) RADM2";

#racm
my $racm_species_file = "/work/users/jco/Mechanisms/RACM/carbon_numbers.txt";
my $racm_run = '/work/users/jco/MECCA/RACM_tagging/boxmodel';
my $racm_mecca = MECCA->new($racm_run);
my $racm_eqnfile = '/work/users/jco/MECCA/RACM_tagging/gas.eqn';
my $racm_kpp = KPP->new($racm_eqnfile);
my $racm_ro2file = "/work/users/jco/MECCA/RACM_tagging/RO2_species.txt";
my @racm_no2_reservoirs = get_no2_reservoirs($racm_kpp, $racm_ro2file);
$families{'Ox_racm'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @racm_no2_reservoirs ];
$weights{'Ox_racm'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'RACM'}) = carbons_others($racm_species_file);
($production_rates{'Ox_racm'}) = get_production_rates($RACM_species, 'Ox_racm', $n_carbon{'RACM'}, $racm_kpp, $racm_mecca);
my ($racm_plot_data) = get_plot_hash('Ox_racm', $production_rates{'Ox_racm'});
my $racm_plot_title = "(e) RACM";

#racm2
my $racm2_species_file = "/work/users/jco/Mechanisms/RACM2/carbon_numbers.txt";
my $racm2_run = '/work/users/jco/MECCA/RACM2_tagged/boxmodel';
my $racm2_mecca = MECCA->new($racm2_run);
my $racm2_eqnfile = '/work/users/jco/MECCA/RACM2_tagged/gas.eqn';
my $racm2_kpp = KPP->new($racm2_eqnfile);
my $racm2_ro2file = "/work/users/jco/MECCA/RACM2_tagged/RO2_species.txt";
my @racm2_no2_reservoirs = get_no2_reservoirs($racm2_kpp, $racm2_ro2file);
$families{'Ox_racm2'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @racm2_no2_reservoirs ];
$weights{'Ox_racm2'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'RACM2'}) = carbons_others($racm2_species_file);
($production_rates{'Ox_racm2'}) = get_production_rates($RACM_species, 'Ox_racm2', $n_carbon{'RACM2'}, $racm2_kpp, $racm2_mecca);
my ($racm2_plot_data) = get_plot_hash('Ox_racm2', $production_rates{'Ox_racm2'});
my $racm2_plot_title = "(f) RACM2";

#CBM-IV
my $cbm4_species_file = "/work/users/jco/Mechanisms/CBM-IV/carbons.txt";
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA->new($cbm4_run);
my $cbm4_eqnfile = "/work/users/jco/MECCA/CBM4_tagging/gas.eqn";
my $cbm4_kpp = KPP->new($cbm4_eqnfile);
my $cbm4_ro2file = "/work/users/jco/MECCA/CBM4_tagging/RO2_species.txt";
my @cbm4_no2_reservoirs = get_no2_reservoirs($cbm4_kpp, $cbm4_ro2file);
$families{'Ox_cbm4'} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5), @cbm4_no2_reservoirs ];
$weights{'Ox_cbm4'} = { NO3 => 2, N2O5 => 3 };

($n_carbon{'CBM-IV'}) = carbons_others($cbm4_species_file);
($production_rates{'Ox_cbm4'}) = get_production_rates($species, 'Ox_cbm4', $n_carbon{'CBM-IV'}, $cbm4_kpp, $cbm4_mecca);
my ($cbm4_plot_data) = get_plot_hash('Ox_cbm4', $production_rates{'Ox_cbm4'});
my $cbm4_plot_title = "(h) CBM-IV";

#CB05
my $cb05_species_file = "/work/users/jco/Mechanisms/CB05/carbons.txt";
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA->new($cb05_run);
my $cb05_eqnfile = "/work/users/jco/MECCA/CB05_tagging/gas.eqn";
my $cb05_kpp = KPP->new($cb05_eqnfile);
my $cb05_ro2file = "/work/users/jco/MECCA/CB05_tagging/RO2_species.txt";
my @cb05_no2_reservoirs = get_no2_reservoirs($cb05_kpp, $cb05_ro2file);
$families{'Ox_cb05'} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5), @cb05_no2_reservoirs ];
$weights{'Ox_cb05'} = { NO3 => 2, N2O5 => 3 };

($n_carbon{'CB05'}) = carbons_others($cb05_species_file);
($production_rates{'Ox_cb05'}) = get_production_rates($species, 'Ox_cb05', $n_carbon{'CB05'}, $cb05_kpp, $cb05_mecca);
my ($cb05_plot_data) = get_plot_hash('Ox_cb05', $production_rates{'Ox_cb05'});
my $cb05_plot_title = "(i) CB05";

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
        mcm3_1_data     => $mcm_3_1_plot_data,
        mcm3_1_title    => $mcm_3_1_plot_title,
        mcm3_2_data     => $mcm_3_2_plot_data,
        mcm3_2_title    => $mcm_3_2_plot_title,
        cri_data        => $cri_plot_data,
        cri_title       => $cri_plot_title,
        mozart_data     => $mozart_plot_data,
        mozart_title    => $mozart_plot_title,
        radm2_data      => $radm2_plot_data,
        radm2_title     => $radm2_plot_title,
        racm_data       => $racm_plot_data,
        racm_title      => $racm_plot_title,
        racm2_data      => $racm2_plot_data,
        racm2_title     => $racm2_plot_title, 
        cbm4_data       => $cbm4_plot_data,
        cbm4_title      => $cbm4_plot_title, 
        cb05_data       => $cb05_plot_data,
        cb05_title      => $cb05_plot_title, 
});

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub mozart_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my $words = join ',', (<$in>);
    close $in;
    my ($string) = $words =~ /Solution(.*?)End\sSolution/s;
    $string =~ s/^\s+,//;
    $string =~ s/,\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/RO2(.*?)->//;
    $string =~ s/ROOH(.*?)->//;
    my @species = split ',', $string;
    my %carbons;
    foreach my $species (@species) {
        $species =~ s/^\s+|\s+$//g;
        my $C_number = 0;
        if ($species !~ /->/ and $species !~ /(C[0-9])/) {
            $C_number ++ while ($species =~ m/C/g);
            $carbons{$species} = $C_number;
        } elsif ($species !~ /->/ and $species =~ /(C[0-9])/) { 
            my ($c_nr) = $species =~ /(C[0-9]+)/s;
            $c_nr =~ s/C//; 
            $C_number = $c_nr;
            $carbons{$species} = $C_number;
        } else {
            my ($mech, $molecule) = split ' -> ', $species;
            $mech =~ s/^\s+|\s+$//g;
            if ($molecule =~ /(C[0-9]+)/) { 
                my ($c_nr) = $molecule =~ /(C[0-9]+)/s;
                $c_nr =~ s/C//; 
                $C_number = $c_nr;
                $carbons{$mech} = $C_number;
            } else {
                $C_number ++ while ($molecule =~ m/C/g);
                $carbons{$mech} = $C_number;
            }
        }
    } 
    return \%carbons;
}

sub carbons_others { #get C-number from file names that have species and C# separated by space
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Cannot open file $file: $!";
    my (@lines) = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
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

sub get_production_rates { #get production rates normalised by emissions of the reactions of the intermediate species of a VOC
    my ($parent_species, $species, $n_carbons, $kpp, $mecca) = @_;
    my %carbons = %$n_carbons;
    
    my @families = ($species, 'HO2x');
    my ($producers, $producer_yields, %carbons_production_rates);
    foreach my $item (@families) {
        if (exists $families{$item}) { #get family reaction numbers and yields
            $kpp->family({ 
                    name    => $item,
                    members => $families{$item},
                    weights => $weights{$item},
            });
            $producers = $kpp->producing($item);
            $producer_yields = $kpp->effect_on($item, $producers);  
        } else { 
            die "No family found!\n";
        }

        #check that species reactions are found
        die "No producers found for $item\n" if (@$producers == 0);

        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $parent_species);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $species =~ /radm2|racm|racm2|cbm4|cb05/);
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) { 
                        $carbons_production_rates{"C$carbons{$lookup}"} += $rate(1:$ntime-2);
                    } elsif ($lookup =~ /CO/) {
                        $carbons_production_rates{'C1'} += $rate(1:$ntime-2);
                    } else {
                        print "nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $carbons_production_rates{"C$carbons{$lookup}"} += $rate(1:$ntime-2);
                }
            } 
        } 
    }

    #operator allocation for those mechanisms that use it: RADM2, RACM, RACM2, CBM4, CB05 -> XO2
    if ($species =~ /radm2|racm|racm2|cbm4|cb05/) {
        my $operator = "XO2_" . $parent_species;
        my $op_producers = $kpp->producing($operator);
        my $op_producer_yields = $kpp->effect_on($operator, $op_producers); 
        die "No producers found for $operator\n" if (@$op_producers == 0);

        for (0..$#$op_producers) { #get rates for all producing reactions
            my $reaction = $op_producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) {
                    $carbons_production_rates{"C$carbons{$lookup}"} += $rate(1:$ntime-2);
                    } else {
                        print "$species => nothing found for $lookup\n";
                    }
                }
            } 
        } 
    }

    #get parent species emissions for each mechanism
    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($species =~ /cbm4/ or $species =~ /cb05/) {
        my $name = "PAR_NC5H12";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
    } else {
        my $parent_source = $mecca->balance($parent_species); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (Product) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $carbons_production_rates{$_} /= $parent_emissions foreach (sort keys %carbons_production_rates);
    
    return \%carbons_production_rates;
}

sub get_plot_hash { #create hash with production of the reactions
    my ($species, $production_rates) = @_;
    my (@plot_data, @sorted_plot_data);

    my @sorted_production = sort { $a cmp $b } keys %$production_rates; #sort data 
    push @sorted_plot_data, { $_ => $production_rates->{$_} } foreach (@sorted_production); #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 

    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    }
    return \@plot_data;
}

sub plot { #create dataframe and then create plot
    my ($args) = @_;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(Cairo) `);
    $R->run(q` library(RColorBrewer) `);

    $R->set('y.max', $args->{y_max});
    $R->set('y.breaks', $args->{breaks});
    $R->set('time', [@{$args->{times}}]);
    
    $R->run(q` my.colours = c("C8" = "#6d6537", "C7.75" = "#4b9483", "C7.1" = "#6c254f", "C7" = "#ef6638", "C6.6" = "#6db875", "C6" = "#76afca", "C5.6" = "#8d1435", "C5" = "#f9c500", "C4.8" = "#603912", "C4.5" = "#86b650", "C4.2" = "#f8c56c", "C4" = "#2c9daf", "C3.9" = "#58691b", "C3.6" = "#b569b3", "C3.5" = "#c9a415", "C3" = "#1c3e3d", "C2.9" = "#cc6329", "C2.4" = "#e7e85e", "C2" = "#dc3522", "C1" = "#0c3f78") `);

    #general R plot function
    $R->run(q` mech.plot = function(data, plot.title) { data = ddply(data, .(time), colwise(sum)) ;
                                                        subset = data[1:7,] ;
                                                        plot.data = melt( data = subset, id = names(subset)[1], measured = names(subset)[-1] ); 
                                                        colnames(plot.data) = c("time", "C.number", "rate"); 
                                                        C.number.levels = (levels(factor(plot.data$C.number))); 
                                                        plot.data$C.number = ordered(plot.data$C.number, levels = C.number.levels); 
                                                        plot.data = ddply( plot.data, .(C.number)); 
                                                        plot = ggplot(data = plot.data, aes(x = time, y = rate)); 
                                                        plot = plot + geom_bar(aes(fill = C.number), stat = "identity", width = 0.5) ;
                                                        plot = plot + ggtitle(plot.title); 
                                                        plot = plot + scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) ;
                                                        plot = plot + scale_fill_manual( limits = rev(C.number.levels), values = my.colours); 
                                                        plot = plot + theme_bw() ; 
                                                        plot = plot + theme(legend.key.size = unit(8.0, "cm")) ; 
                                                        plot = plot + theme(axis.text.x = element_text(size = 120)) ; 
                                                        plot = plot + theme(axis.text.y = element_text(size = 130)) ; 
                                                        plot = plot + theme(axis.title.x = element_blank()) ; 
                                                        plot = plot + theme(legend.text = element_text(size = 150)) ; 
                                                        plot = plot + theme(legend.key = element_blank()) ; 
                                                        plot = plot + theme(axis.title.y = element_blank()) ; 
                                                        plot = plot + theme(plot.title = element_text(size = 160, face = "bold", vjust = 1)) ; 
                                                        plot = plot + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) ;
                                                        plot = plot + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) ;
                                                        plot = plot + theme(legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95)) ; 
                                                        plot = plot + theme(legend.title = element_blank()); 
                                                        return(plot) } `);

    #MCM 3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.plot.title', $args->{mcm3_1_title});
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate*10000 `); 
        }
    } 
    $R->run(q` mcm3.1.plot = mech.plot(mcm3.1.data, mcm3.1.plot.title) `);

    #MCM v3.2
    $R->run(q` mcm3.2.data = data.frame(time)`);
    $R->set('mcm3.2.plot.title', $args->{mcm3_2_title});
    foreach my $ref (@{$args->{mcm3_2_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.2.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` mcm3.2.plot = mech.plot(mcm3.2.data, mcm3.2.plot.title) `);

    #CRI
    $R->set('cri.plot.title', $args->{cri_title});
    $R->run(q` cri.data = data.frame(time) `);
    foreach my $ref (@{$args->{cri_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cri.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` cri.plot = mech.plot(cri.data, cri.plot.title) `);

    #MOZART
    $R->set('mozart.plot.title', $args->{mozart_title});
    $R->run(q` mozart.data = data.frame(time) `);
    foreach my $ref (@{$args->{mozart_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mozart.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` mozart.plot = mech.plot(mozart.data, mozart.plot.title) `);

    #RADM2
    $R->set('radm2.plot.title', $args->{radm2_title});
    $R->run(q` radm2.data = data.frame(time) `);
    foreach my $ref (@{$args->{radm2_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` radm2.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` radm2.plot = mech.plot(radm2.data, radm2.plot.title) `);

    #RACM
    $R->set('racm.plot.title', $args->{racm_title});
    $R->run(q` racm.data = data.frame(time) `);
    foreach my $ref (@{$args->{racm_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` racm.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` racm.plot = mech.plot(racm.data, racm.plot.title) `);

    #RACM2
    $R->set('racm2.plot.title', $args->{racm2_title});
    $R->run(q` racm2.data = data.frame(time) `);
    foreach my $ref (@{$args->{racm2_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` racm2.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` racm2.plot = mech.plot(racm2.data, racm2.plot.title) `);

    #CBM-IV
    $R->set('cbm4.plot.title', $args->{cbm4_title});
    $R->run(q` cbm4.data = data.frame(time) `);
    foreach my $ref (@{$args->{cbm4_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cbm4.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` cbm4.plot = mech.plot(cbm4.data, cbm4.plot.title) `);

    #CB05
    $R->set('cb05.plot.title', $args->{cb05_title});
    $R->run(q` cb05.data = data.frame(time) `);
    foreach my $ref (@{$args->{cb05_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` cb05.data[name] = rate*10000`); 
        }
    } 
    $R->run(q` cb05.plot = mech.plot(cb05.data, cb05.plot.title) `);

    $R->run(q` CairoPDF(file = "pentane_carbon_breakdown.pdf", width = 130, height = 130) `,
            q` y.label = textGrob(expression(bold(paste("Molecules (Product) ", s^-1, "/Molecules (NMVOC) x ", 10^4))), rot = 90, gp = gpar(fontsize = 150), vjust = 0.5)`,
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
                                       sub = textGrob("\n", gp = gpar(fontsize = 20)), 
                                       widths=unit.c(unit(16, "lines"), unit(1, "npc") - unit(16, "lines"))) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );
    $R->stop(); 
}
