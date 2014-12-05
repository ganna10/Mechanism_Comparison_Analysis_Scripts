#!/usr/bin/perl
# plot mixing ratio of O3, NO, NO2, OH, PAN, Radicals
# Version 0: Jane Coates 11/02/2014
# Version 1: Jane Coates 05/05/2014 Changed to ppb not mol/mol

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use PDL::NetCDF;
use Statistics::R;

my %concentration;
my @species = qw( O3 );

#MCMv3.1 data
my $mcm_3_1_run = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/boxmodel";
my $mcm_3_1_species = "/work/users/jco/MECCA/MCM_3.1_tagged_3.2rates/gas.spc";
my $mcm_3_1_mecca = MECCA->new($mcm_3_1_run); 
my $ntime = $mcm_3_1_mecca->time->nelem; #number of time points

#MCMv3.2 data
my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_species = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.spc";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 

#cri data
my $cri_run = "/work/users/jco/MECCA/CRI_tagging/boxmodel";
my $cri_species = "/work/users/jco/MECCA/CRI_tagging/gas.spc";
my $cri_mecca = MECCA->new($cri_run); 

#mozart data
my $mozart_run = "/work/users/jco/MECCA/MOZART_tagging/boxmodel";
my $mozart_species = "/work/users/jco/MECCA/MOZART_tagging/gas.spc";
my $mozart_mecca = MECCA->new($mozart_run); 

#radm2 data
my $radm2_run = "/work/users/jco/MECCA/RADM2_tagged/boxmodel";
my $radm2_species = "/work/users/jco/MECCA/RADM2_tagged/gas.spc";
my $radm2_mecca = MECCA->new($radm2_run); 

#racm data
my $racm_run = "/work/users/jco/MECCA/RACM_tagging/boxmodel";
my $racm_species = "/work/users/jco/MECCA/RACM_tagging/gas.spc";
my $racm_mecca = MECCA->new($racm_run); 

#racm2 data
my $racm2_run = "/work/users/jco/MECCA/RACM2_tagged/boxmodel";
my $racm2_mecca = MECCA->new($racm2_run); 
my $racm2_species = "/work/users/jco/MECCA/RACM2_tagged/gas.spc";

#CBM-IV
my $cbm4_run = "/work/users/jco/MECCA/CBM4_tagging/boxmodel";
my $cbm4_mecca = MECCA->new($cbm4_run);
my $cbm4_species = "/work/users/jco/MECCA/CBM4_tagging/gas.spc";

#CB05
my $cb05_run = "/work/users/jco/MECCA/CB05_tagging/boxmodel";
my $cb05_mecca = MECCA->new($cb05_run);
my $cb05_species = "/work/users/jco/MECCA/CB05_tagging/gas.spc";

foreach my $species (@species) {
    my $mcm_3_1_proper_species = get_tagged_species($species, $mcm_3_1_species);
    ($concentration{"${species}_MCM v3.1"}) = get_concentration($mcm_3_1_proper_species, $mcm_3_1_mecca);
    my $mcm_3_2_proper_species = get_tagged_species($species, $mcm_3_2_species);
    ($concentration{"${species}_MCM v3.2"}) = get_concentration($mcm_3_2_proper_species, $mcm_3_2_mecca);
    my $cri_proper_species = get_tagged_species($species, $cri_species);
    ($concentration{"${species}_CRI v2"}) = get_concentration($cri_proper_species, $cri_mecca);
    my $mozart_proper_species = get_tagged_species($species, $mozart_species);
    ($concentration{"${species}_MOZART-4"}) = get_concentration($mozart_proper_species, $mozart_mecca);
    my $radm2_proper_species = get_tagged_species($species, $radm2_species);
    ($concentration{"${species}_RADM2"}) = get_concentration($radm2_proper_species, $radm2_mecca);
    my $racm_proper_species = get_tagged_species($species, $racm_species);
    ($concentration{"${species}_RACM"}) = get_concentration($racm_proper_species, $racm_mecca);
    my $racm2_proper_species = get_tagged_species($species, $racm2_species);
    ($concentration{"${species}_RACM2"}) = get_concentration($racm2_proper_species, $racm2_mecca);
    my $cbm4_proper_species = get_tagged_species($species, $cbm4_species);
    ($concentration{"${species}_CBM-IV"}) = get_concentration($cbm4_proper_species, $cbm4_mecca) ;
    my $cb05_proper_species = get_tagged_species($species, $cb05_species);
    ($concentration{"${species}_CB05"}) = get_concentration($cb05_proper_species, $cb05_mecca);
}

#Create x-axis for plot in hours
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list;

concentration_plot(\@time_axis, \%concentration);

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

sub get_concentration {
    my ($species, $mecca) = @_;
    
    my ($concs, @conc_array);
#    if (exists $families{$species}) { #radicals family
#        foreach my $item (@{$families{$species}}){ 
#            next unless (exists $tracers{$item});
#            $concs += $mecca->tracer($item);
#        }
#    } else {
        foreach my $item (@$species){ 
            #next unless (exists $tracers{$item});
            $concs += $mecca->tracer($item);
        }
#    }
    if (defined $concs) {
        $concs = $concs(1:$ntime-2) * 1e9;
        @conc_array = map { $_ } $concs->dog;
    }
    return \@conc_array;
}

sub get_tagged_species {
    my ($species, $spc_file) = @_;

    open FILE, "<$spc_file" or die $!;
    my @species = grep { $_ =~ /^[A-Z]/ } (<FILE>);
    close FILE; 
    my @all_species = map { $_ =~ /^(.*?)\s/ } @species;
    my @tagged_species = grep { $_ =~ /^${species}_/ } @all_species;
    push @tagged_species, $species;
    return \@tagged_species;
}

sub concentration_plot {
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
            q` Mixing.ratio = {} `,
    );

    foreach my $item (sort keys %data) {
        my ($name, $mechanism) = split /_/, $item;
        my $R_name = $R->set('name', $name);
        my $R_mech = $R->set('mechanism', $mechanism);
        my $R_data = $R->set('mixing.ratio', [@{$data{$item}}]);
        $R->run(q` Species = cbind(Species, rep(name, length(time))) `,
                q` Mechanism = cbind(Mechanism, rep(mechanism, length(time))) `,
                q` Mixing.ratio = cbind(Mixing.ratio, mixing.ratio) `,
        );
    }

    #create dataframe after converting the matrices above to vectors
    $R->run(q` Species = c(Species) `,
            q` Mechanism = c(Mechanism) `,
            q` Mixing.ratio = c(Mixing.ratio) `,
            q` data = data.frame(times, Species, Mechanism, Mixing.ratio) `,
    );
    
    $R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `), #scientific label format for y-axis
    $R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);
    
    #plot
    $R->run(q` plot.lines = function () { list( geom_line(size = 3), 
                                                scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), 
                                                theme_bw(), 
                                                theme(axis.text.x = element_text(size = 40)), 
                                                theme(axis.text.y = element_text(size = 40)), 
                                                theme(legend.title = element_text(size = 60, face = "bold")), 
                                                theme(legend.text = element_text(size = 50)), 
                                                theme(legend.key.size = unit(4, "cm")), 
                                                theme(axis.title.x = element_blank()), 
                                                theme(axis.title.y = element_blank()), 
                                                scale_colour_manual(values = my.colours), 
                                                theme(legend.key = element_blank()), 
                                                theme(plot.title = element_text(size = 60, face = "bold", vjust = 1)) ) } `,
    );

    $R->run(q` plot5 = ggplot(data = subset(data, Species == "O3"), aes(x = times, y = Mixing.ratio, colour = Mechanism)) `,
            q` plot5 = plot5 + plot.lines() `, 
            q` plot5 = plot5 + ggtitle("O3") `,

            q` CairoPDF(file = "mixing_ratios.pdf", width = 40, height = 20) `,
            q` print(plot5) `,
            q` dev.off() `,
    );

    $R->stop();
} 
