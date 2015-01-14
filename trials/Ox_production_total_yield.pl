#! /usr/bin/env perl
# Plot Total Ox yield for each alkane 
# Version 0: Jane Coates 14/1/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
#my @mechanisms = qw( RADM2 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    ($data{$mechanism}) = get_data($kpp, $mecca, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
);

foreach my $mechanism (sort keys %data) {
    print "$mechanism\n";
    foreach my $alk (sort keys %{$data{$mechanism}}) {
        print "\t$alk => $data{$mechanism}{$alk}\n";
    }
}

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, $producer_yields, %production_rates, $consumers, $consumer_yields, %consumption_rates);
    my @families = ("Ox_${mechanism}", "HO2x");
    my %yield;
    foreach my $family (@families) {
        if (exists $families{$family}) { 
            $kpp->family({ 
                    name    => $family,
                    members => $families{$family},
                    weights => $weights{$family},
            });
            $producers = $kpp->producing($family);
            $producer_yields = $kpp->effect_on($family, $producers);  
            $consumers = $kpp->consuming($family);
            $consumer_yields = $kpp->effect_on($family, $consumers);  
        } else {
            print "No family found for $family\n";
        } 
        die "No producers found for $family\n" if (@$producers == 0);
        die "No consumers found for $family\n" if (@$consumers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent =~ /C2H6|C3H8|NC|IC|BIGALK|HC|ETH/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $yield{$parent} += $rate(1:$NTIME-2);
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent =~ /C2H6|C3H8|NC|IC|BIGALK|HC|ETH/);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $yield{$parent} += $rate(1:$NTIME-2);
        }
    }

    foreach my $alkane (sort keys %yield) {
        my $reshape = $yield{$alkane}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $yield{$alkane} = $integrate->at(0);
    }
    return \%yield;
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
