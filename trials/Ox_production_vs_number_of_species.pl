#! /usr/bin/env perl
# Correlate alkane total Ox production to number of Ox species and yield
# Version 0: Jane Coates 26/12/2014

use strict;
use diagnostics;
use KPP;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 5 };
    $data{$mechanism} = get_data($kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` data = data.frame(Mechanism = as.numeric(0), VOC = as.numeric(0), Yield = as.numeric(0)) `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $parent (sort keys %{$data{$mechanism}}) {
        $R->set('voc', $parent);
        $R->set('yield', $data{$mechanism}{$parent});
        $R->run(q` data = rbind(data, c(mechanism, voc, yield)) `);
    }
}
$R->run(q` data = data[-1,] `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Yield)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + facet_wrap( ~ VOC )`,
        q` plot = plot + coord_flip() `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05"))) `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_total_yield.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
    my @loop = ("Ox_$mechanism", "HO2x");
    my (%production, %consumption);

    foreach my $species (@loop) {
        my ($producers, $producer_yields, $consumers, $consumer_yields);
        if (exists $families{$species}) {
            $kpp->family({
                    name    => $species,
                    members => $families{$species},
                    weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        } else {
            print "No family found for $species\n";
        }
        print "No producers found for $species\n" if (@$producers == 0);
        print "No consumers found for $species\n" if (@$consumers == 0);
        
        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            next unless ($reaction =~ /NC|IC|C2H6|C3H8|BIGALK|HC|ETH/);
            my ($number, $parent) = split /_/, $reaction;
            my $yield = $producer_yields->[$_];
            $production{$parent} += $yield;
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            next unless ($reaction =~ /NC|IC|C2H6|C3H8|BIGALK|HC|ETH/);
            my ($number, $parent) = split /_/, $reaction;
            my $yield = $consumer_yields->[$_];
            $consumption{$parent} += $yield;
        }
    }

    foreach my $parent (keys %production) {
        $production{$parent} += $consumption{$parent};
        if ($mechanism eq "MOZART-4") {
            if ($parent eq "BIGALK") {
                $production{"NC4H10"} = 0.285 * $production{$parent};
                $production{"IC4H10"} = 0.151 * $production{$parent};
                $production{"NC5H12"} = 0.146 * $production{$parent};
                $production{"IC5H12"} = 0.340 * $production{$parent};
                $production{"NC6H14"} = 0.048 * $production{$parent};
                $production{"NC7H16"} = 0.020 * $production{$parent};
                $production{"NC8H18"} = 0.010 * $production{$parent};
                delete $production{$parent};
            }
        } elsif ($mechanism eq "RADM2" or $mechanism =~ /RACM/) {
            if ($parent eq "ETH") {
                $production{"C2H6"} = $production{$parent};
                delete $production{$parent};
            } elsif ($parent eq "HC8") {
                $production{"NC8H18"} = $production{$parent};
                delete $production{$parent};
            } elsif ($parent eq "HC3") {
                $production{"C3H8"} = 0.628 * $production{$parent};
                $production{"NC4H10"} = 0.243 * $production{$parent};
                $production{"IC4H10"} = 0.129 * $production{$parent};
                delete $production{$parent};
            } elsif ($parent eq "HC5") {
                $production{"NC5H12"} = 0.264 * $production{$parent};
                $production{"IC5H12"} = 0.615 * $production{$parent};
                $production{"NC6H14"} = 0.086 * $production{$parent};
                $production{"NC7H16"} = 0.035 * $production{$parent};
                delete $production{$parent};
            } 
        }
    }

    if ($mechanism eq "CBM-IV") {
        $production{"C2H6"} /= 0.4;
    }
    if ($mechanism =~ /CB/) {
        $production{"C3H8"} /= 1.5;
        $production{"NC4H10"} /= 4;
        $production{"IC4H10"} /= 4;
        $production{"NC5H12"} /= 5;
        $production{"IC5H12"} /= 5;
        $production{"NC6H14"} /= 6;
        $production{"NC7H16"} /= 7;
        $production{"NC8H18"} /= 8;
    }
    return \%production;
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
