#! /usr/bin/env perl
# Plot Ox production per yield of each size fragment of degradation products pentane 
# Version 0: Jane Coates 2/1/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2",  "CBM-IV", "CB05" );
#my @mechanisms = qw(CB05 CBM-IV);
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $ro2file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    my @parents = qw( Pentane );
    foreach my $NMVOC (@parents) {
        my $parent = get_mechanism_species($NMVOC, $mechanism);
        ($data{$mechanism}{$NMVOC}) = get_data($kpp, $mecca, $mechanism, $parent);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(scales) `,
        q` library(grid) `,
);

$R->run(q` data = data.frame(Mechanism = as.numeric(0), VOC = as.numeric(0), Yield = as.numeric(0), Day1 = as.numeric(0), Total = as.numeric(0)) `);

foreach my $run (sort keys %data) {
    foreach my $VOC (sort keys %{$data{$run}}) {
        $R->set('mechanism', $run);
        $R->set('Parent', $VOC); 
        foreach my $yield (keys %{$data{$run}{$VOC}}) {
            my $sum = $data{$run}{$VOC}{$yield}->sumover;
            $R->set('yield', $yield);
            $R->set('day.1', $data{$run}{$VOC}{$yield}->at(0));
            $R->set('cumulative', $sum->at(0));
        }
        $R->run(q` data = rbind(data, c(mechanism, Parent, yield, day.1, cumulative)) `);
    }
}
$R->run(q` data = data[-1, ] `);
$R->run(q` data = arrange(data, desc(Yield)) `);
my $p = $R->run(q` print(data) `);
print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Yield, y = Day1, colour = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + facet_wrap( ~ VOC ) `,
);

$R->run(q` CairoPDF(file = "Total_Yield_vs_Ox_day1_Production.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Yield, y = Total, colour = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + facet_wrap( ~ VOC ) `,
);

$R->run(q` CairoPDF(file = "Total_Yield_vs_Ox_Total_Production.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC) = @_;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, $producer_yields, $production_rate, $yield, $consumers, $consumer_yields, $consumption_rate);
    my @families = ("Ox_${mechanism}_$VOC", "HO2x_${mechanism}_$VOC");
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
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $yield += $producer_yields->[$_];
            $production_rate += $rate(1:$NTIME-2);
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $yield += $consumer_yields->[$_];
            $consumption_rate += $rate(1:$NTIME-2);
        }
    }

    $production_rate += $consumption_rate;

    if ($mechanism =~ /MOZ/) {
        $production_rate *= 0.146;
    } elsif ($mechanism =~ /RADM2/ or $mechanism =~ /RACM/) {
        $production_rate *= 0.264;
    } elsif ($mechanism =~ /CB/) {
        $production_rate /= 5;
    } 
    my $reshape = $production_rate->reshape($N_PER_DAY, $N_DAYS);
    my $integrate = $reshape->sumover;
    $integrate = $integrate(0:13:2);
    $production_rate = $integrate;
    my %data;
    $data{$yield} = $production_rate;
    return \%data;
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

sub get_mechanism_species {
    my ($NMVOC, $run) = @_;

    my $mechanism_species;
    if ($NMVOC eq "Pentane") {
        if ($run =~ /MCM|CRI|CB/) {
            $mechanism_species = "NC5H12";
        } elsif ($run =~ /MOZART/) {
            $mechanism_species = "BIGALK";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "HC5";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } elsif ($NMVOC eq "Toluene") {
        if ($run =~ /MCM|CRI|MOZART|CB/) {
            $mechanism_species = "TOLUENE";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "TOL";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } else {
        print "No $NMVOC data\n";
    }
    return $mechanism_species;
}
