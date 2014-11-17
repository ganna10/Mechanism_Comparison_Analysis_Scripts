#! /usr/bin/env perl
# Calculate OxPE of VOCs in each mechanism by normalising over total Ox consumption
# Version 0: Jane Coates 14/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT ;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#my @runs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "MCM v3.2", "MCM v3.1", "CRI v2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my $index = 0;
my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $VOC (sort keys %{$plot_data{$mechanism}}) {
        my $chemical = get_species_name($VOC);
        $R->set('voc', $chemical);
        $R->set('oxpe', [map { $_ } $plot_data{$mechanism}{$VOC}->dog] );
        $R->run(q` pre[voc] = oxpe `);
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "VOC", value.name = "OxPE") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRI v2" = "#ef6638", "MCM v3.1" = "#000000", "MCM v3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap(~ VOC, scales = "free") `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "VOC_OxPEs.pdf", width = 10, height = 10) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_species_name {
    my ($VOC) = @_;
    if ($VOC =~ /C2H6|ETH/) {
        $VOC = "Ethane";
    } elsif ($VOC =~ /C3H8|HC3/) {
        $VOC = "Propane";
    } elsif ($VOC =~ /NC4H10/) {
        $VOC = "Butane";
    } elsif ($VOC =~ /IC4H10/) {
        $VOC = "2-Methylpropane";
    } elsif ($VOC =~ /NC5H12|HC5|BIGALK/) {
        $VOC = "Pentane";
    } elsif ($VOC =~ /IC5H12/) {
        $VOC = "2-Methylbutane";
    } elsif ($VOC =~ /NC6H14/) {
        $VOC = "Hexane";
    } elsif ($VOC =~ /NC7H16/) {
        $VOC = "Heptane";
    } elsif ($VOC =~ /NC8H18|HC8/) {
        $VOC = "Octane";
    } elsif ($VOC =~ /C2H4|OL2|ETE/) {
        $VOC = "Ethene";
    } elsif ($VOC =~ /C3H6|OLT/) {
        $VOC = "Propene";
    } elsif ($VOC =~ /BUT1ENE|BIGENE/) {
        $VOC = "Butene";
    } elsif ($VOC =~ /MEPROPENE|OLI/) {
        $VOC = "2-Methylpropene";
    } elsif ($VOC =~ /C5H8|ISO/) {
        $VOC = "Isoprene";
    } elsif ($VOC =~ /BENZENE|BEN/) {
        $VOC = "Benzene";
    } elsif ($VOC =~ /TOLUENE|TOL/) {
        $VOC = "Toluene";
    } elsif ($VOC =~ /XYO|OXYL/) {
        $VOC = "o-Xylene";
    } elsif ($VOC =~ /XYP|PXYL/) {
        $VOC = "p-Xylene";
    } elsif ($VOC =~ /MXYL|XYL|XYM/) {
        $VOC = "m-Xylene";
    } elsif ($VOC =~ /EBENZ/) {
        $VOC = "Ethylbenzene";
    } elsif ($VOC =~ /CH4/) {
        $VOC = "Methane";
    } elsif ($VOC =~ /Inorganic/) {
        $VOC = $VOC;
    } else {
        print "No VOC for $VOC\n";
    }
    return $VOC;
}

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
    my @loop = ("Ox_$mechanism", "HO2x");
    my (%production, $consumption);

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

        foreach (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($r_number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production{$parent} += $rate(1:$NTIME-2);
            } else {
                $production{"Inorganic"} += $rate(1:$NTIME-2);
            }
        }

        foreach (0..$#$consumers) {
            next if ($species eq "HO2x");
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            $consumption += $rate(1:$NTIME-2);
        }
    }

    $production{$_} /= -$consumption foreach (keys %production); #normalise each processes' Ox production by total consumption

    foreach my $VOC (keys %production) {
        my $reshape = $production{$VOC}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$VOC} = $integrate;
    }
    return \%production;
}

sub get_no2_reservoirs {
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
