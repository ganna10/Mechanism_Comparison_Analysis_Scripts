#! /usr/bin/perl
# Analyse number of Carbon atoms per radical produced during Toluene breakdown, compare for all mechanisms
# Version 0: Jane Coates 13/02/2014

use strict;
use diagnostics;
use KPP;
use MECCA_TIM;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my (%n_carbon, %families, %weights, %production_rates);
$families{'HO2x'} = [ qw( HO2 HO2NO2 ) ];

my $species = "TOLUENE";
my $RACM_species = "TOL";

#mcm 3.1
my $mcm_3_1_species_file = "/work/users/jco/Mechanisms/MCM/MCM-3.1/species_smiles.txt";
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_mecca = MECCA_TIM->new($mcm_3_1_run); 
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
my $mcm_3_2_mecca = MECCA_TIM->new($mcm_3_2_run); 
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
my $mozart_mecca = MECCA_TIM->new($mozart_run);
my $mozart_eqnfile = '/work/users/jco/MECCA/MOZART_tagging/gas.eqn';
my $mozart_kpp = KPP->new($mozart_eqnfile);
my $mozart_ro2file = "/work/users/jco/MECCA/MOZART_tagging/RO2_species.txt";
my @mozart_no2_reservoirs = get_no2_reservoirs($mozart_kpp, $mozart_ro2file);
$families{'Ox_mozart'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @mozart_no2_reservoirs ];
$weights{'Ox_mozart'} = { NO3 => 2, N2O5 => 3};

($n_carbon{'MOZART'}) = mozart_n_carbon($mozart_species_file);
($production_rates{'Ox_mozart'}) = get_production_rates($species, 'Ox_mozart', $n_carbon{'MOZART'}, $mozart_kpp, $mozart_mecca);
my ($mozart_plot_data) = get_plot_hash('Ox_mozart', $production_rates{'Ox_mozart'});
my $mozart_plot_title = "(g) MOZART-4";

#cri
my $cri_species_file = "/work/users/jco/Mechanisms/CRI/CRI_v2_full/carbons.txt";
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_mecca = MECCA_TIM->new($cri_run);
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
my $radm2_mecca = MECCA_TIM->new($radm2_run);
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
my $racm_mecca = MECCA_TIM->new($racm_run);
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
my $racm2_mecca = MECCA_TIM->new($racm2_run);
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
my $cbm4_mecca = MECCA_TIM->new($cbm4_run);
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
my $cb05_mecca = MECCA_TIM->new($cb05_run);
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
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list; 

my ($plot) = plot({#create plot
        y_max           => 20,
        breaks          => 2,
        times           => \@time_axis,
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

    open FILE, "<$file" or die $!;
    my @lines = (<FILE>);
    close FILE;
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

    open FILE, "<$file" or die $!;
    my $words = join ',', (<FILE>);
    close FILE;
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

    open FILE, "<$file" or die "Cannot open file $file: $!";
    my (@lines) = (<FILE>);
    close FILE;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
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
                    next if ($_ =~ /XO2/);
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup} and $item =~ /Ox/) {
                        $carbons_production_rates{"C$carbons{$lookup}"} += $rate(1:$ntime-2);
                    } elsif (defined $carbons{$lookup} and $item =~ /HO2x/) {
                        $carbons_production_rates{'HO2'} += $rate(1:$ntime-2);
                    } else {
                        print "nothing found for $lookup\n";
                    }
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
    $carbons_production_rates{$_} /= $parent_emissions foreach (sort keys %carbons_production_rates);
    
    return \%carbons_production_rates;
}

sub get_plot_hash { #create hash with production of the reactions
    my ($species, $production_rates) = @_;
    my %production_rates = %$production_rates;
    my @plot_data;

    #sort data
    my $HO2_data = $production_rates{'HO2'};
    delete $production_rates{'HO2'};
    my @sorted_production = sort { $a cmp $b } keys %production_rates;
    my @sorted_plot_data;
    push @sorted_plot_data, { $_ => $production_rates{$_} } foreach (@sorted_production); #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
    unshift @sorted_plot_data, { 'HO2' => $HO2_data };

    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            my $rate_list = join ":", $ref->{$item}->dog;
            my @rate_array = split /:/, $rate_list;
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
    $R->run(q` my.colours = c("C8" = "#6d6537", "C7.75" = "#4b9483", "C7.1" = "#6c254f", "C7" = "#ef6638", "C6.6" = "#6db875", "C6" = "#76afca", "C5.6" = "#8d1435", "C5" = "#f9c500", "C4.8" = "#603912", "C4.5" = "#86b650", "C4.2" = "#f8c56c", "C4" = "#2c9daf", "C3.9" = "#58691b", "C3.6" = "#b569b3", "C3.5" = "#c9a415", "C3" = "#1c3e3d", "C2.9" = "#cc6329", "C2" = "#e7e85e", "C1" = "#dc3522", "HO2" = "#0c3f78") `);

    #general R plot function
    $R->run(q` mech.plot = function(data, plot.title) { plot.data = melt( data = data, id = names(data)[1], measured = names(data)[-1] ); colnames(plot.data) = c("time", "C.number", "rate"); C.number.levels = (levels(factor(plot.data$C.number))); plot.data$C.number = ordered(plot.data$C.number, levels = C.number.levels); plot.data = ddply( plot.data, .(C.number)); plot = ggplot(data = plot.data, aes(x = time, y = rate, fill = C.number)); plot = plot + geom_area(position = "stack", alpha = 1); plot = plot + geom_line(position = "stack"); plot = plot + ggtitle(plot.title); plot = plot + scale_x_continuous(limits=c(0, 7), breaks=seq(0, 7, 1)); plot = plot + scale_y_continuous(limits=c(0, y.max), breaks=seq(0, y.max, y.breaks)); plot = plot + scale_fill_manual( limits = rev(C.number.levels), values = my.colours); plot = plot + theme_bw() ; plot = plot + theme(legend.key.size = unit(4.5, "cm")) ; plot = plot +  theme(axis.text.x = element_text(size = 120)) ; plot = plot + theme(axis.text.y = element_text(size = 120)) ; plot = plot + theme(axis.title.x = element_blank()) ; plot = plot + theme(legend.text = element_text(size = 100)) ; plot = plot + theme(legend.key = element_blank()) ; plot = plot + theme(axis.title.y = element_blank()) ; plot = plot + theme(plot.title = element_text(size = 130, face = "bold", vjust = 1)) ; plot = plot + theme(legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95)) ; plot = plot + theme(legend.title = element_blank()); return(plot) } `);

    #MCM 3.1
    $R->run(q` mcm3.1.data = data.frame(time)`);
    $R->set('mcm3.1.plot.title', $args->{mcm3_1_title});
    foreach my $ref (@{$args->{mcm3_1_data}}) {
        for my $key (keys %$ref) {
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@{$ref->{$key}}]);
            $R->run(q` mcm3.1.data[name] = rate*100000 `); 
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
            $R->run(q` mcm3.2.data[name] = rate*100000`); 
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
            $R->run(q` cri.data[name] = rate*100000`); 
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
            $R->run(q` mozart.data[name] = rate*100000`); 
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
            $R->run(q` radm2.data[name] = rate*100000`); 
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
            $R->run(q` racm.data[name] = rate*100000`); 
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
            $R->run(q` racm2.data[name] = rate*100000`); 
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
            $R->run(q` cbm4.data[name] = rate*100000`); 
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
            $R->run(q` cb05.data[name] = rate*100000`); 
        }
    } 
    $R->run(q` cb05.plot = mech.plot(cb05.data, cb05.plot.title) `);

    $R->run(q` CairoPDF(file = "toluene_carbon_breakdown.pdf", width = 100, height = 90) `,
            q` y.label = textGrob(expression(bold(paste("Molecules (Product) ", s^-1, "/Molecules (NMVOC) x ", 10^5))), rot = 90, gp = gpar(fontsize = 120), vjust = 0.5)`,
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
                                       sub = textGrob("\nTime (days)", gp = gpar(fontsize = 120, fontface = "bold"), vjust = 0), 
                                       widths=unit.c(unit(12, "lines"), unit(1, "npc") - unit(12, "lines"))) `,
            q` print(main.plot) `,
            q` dev.off() `,
    );
    $R->stop(); 
}
