#! /usr/bin/env perl
# Calculate daily TOPP of all VOCs in each mechanism for constant emission runs
# Version 0: Jane Coates 5/12/2014
# Version 1: Jane Coates 17/12/2014 de-lumping VOC into MCM species for RADM2, RACM, RACM2 and MOZART-4
# Version 2: Jane Coates 6/1/2015 re-doing TOPP calculation for lumped species
# Version 3: Jane Coates 13/1/2015 treating de-lumped VOC TOPP values by scaling by carbon number (similar to CB)

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

#my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my @mechanisms = qw( MCMv3.2 );
my (%families, %weights, %data); 
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    ($data{$mechanism}) = get_TOPPs($mecca, $kpp, $mechanism);
}

my $output = "_TOPP_values.txt";
foreach my $mechanism (sort keys %data) {
    my $file = $mechanism . $output;
    open my $out, '>:encoding(utf-8)', $file or die "Can't open $file : $!";
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        my @topps = $data{$mechanism}{$VOC}->dog;
        print $out "$VOC => @topps\n";
    }
    close $out;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` my.colours = c( "Ethane " = "#696537", "Propane " = "#f9c600", "Butane " = "#76afca", "2-Methylpropane " = "#dc3522", "Pentane " = "#8c6238", "2-Methylbutane " = "#9bb08f", "Hexane " = "#8b1537", "Heptane " = "#ba8b01", "Octane " = "#0352cb", "Ethene " = "#86b650", "Propene " = "#6c254f", "Butene " = "#ee6738", "2-Methylpropene " = "#58691b", "Isoprene " = "#8ed6d5", "Benzene " = "#f3aa7f", "Toluene " = "#c65d6c", "m-Xylene " = "#888a87", "o-Xylene " = "#0e5c28", "p-Xylene " = "#b569b3", "Ethylbenzene " = "#2c9def" ) `, 
);
$R->set('Time', [ (1..7) ]);
$R->run(q` plotting = function (data, title, filename) { plot = ggplot(data, aes( x = Time, y = TOPP, colour = VOC, group = VOC));
                                                         plot = plot + geom_line();
                                                         plot = plot + geom_point();
                                                         plot = plot + ggtitle(title);
                                                         plot = plot + theme_bw();
                                                         plot = plot + scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, 1));
                                                         plot = plot + xlab("Time (days)");
                                                         plot = plot + ylab("TOPP (molecules (Ox) / molecules (VOC))");
                                                         plot = plot + theme(axis.title = element_text(face = "bold"));
                                                         plot = plot + theme(plot.title = element_text(face = "bold"));
                                                         plot = plot + theme(panel.grid = element_blank());
                                                         plot = plot + theme(legend.key = element_blank());
                                                         plot = plot + theme(legend.title = element_blank());
                                                         plot = plot + theme(panel.border = element_rect(colour = "black"));
                                                         plot = plot + scale_colour_manual(values = my.colours);
                                                         CairoPDF(file = filename);
                                                         print(plot);
                                                         dev.off() } `);

foreach my $mechanism (sort keys %data) {
    $R->set('daily.filename', "${mechanism}_TOPP_daily.pdf");
    $R->set('daily.title', "$mechanism : TOPP Daily Values");
    $R->set('cumulative.filename', "${mechanism}_TOPP_cumulative.pdf");
    $R->set('cumulative.title', "$mechanism : TOPP Cumulative Values");
    $R->run(q` data = data.frame(Time) `);
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        next if ($VOC eq "CH4");
        my $name = get_chemical_name($VOC);
        $R->set('voc', $name);
        $R->set('topp', [ map { $_ } $data{$mechanism}{$VOC}->dog ]);
        $R->run(q` data[voc] = topp `);
    }
    $R->run(q` daily.data = gather(data, VOC, TOPP, -Time) `,
            q` cumulative = as.data.frame(lapply(data[, -1], cumsum)) `,
            q` colnames(cumulative) = c("Benzene ", "Butene ", "Ethene ", "Ethane ", "Propene ", "Propane ", "Isoprene ", "Ethylbenzene ", "2-Methylpropane ", "2-Methylbutane ", "2-Methypropene ", "m-Xylene ", "Butane ", "Pentane ", "Hexane ", "Heptane ", "Octane ", "o-Xylene ", "p-Xylene ", "Toluene ") `,
            q` cumulative$Time = Time `,
            q` cumulative.data = gather(cumulative, VOC, TOPP, -Time) `,
            q` cumulative.data$VOC = factor(cumulative.data$VOC, levels = c("Ethane ", "Propane ", "Butane ", "2-Methylpropane ", "Pentane ", "2-Methylbutane ", "Hexane ", "Heptane ", "Octane ", "Ethene ", "Propene ", "Butene ", "2-Methylpropene ", "Isoprene ", "Benzene ", "Toluene ", "m-Xylene ", "o-Xylene ", "p-Xylene ", "Ethylbenzene ")) `,
            q` plotting(daily.data, daily.title, daily.filename) `,
            q` plotting(cumulative.data, cumulative.title, cumulative.filename) `,
    );
    #my $p = $R->run(q` print(cumulative.data) `);
    #print $p, "\n";
}

$R->stop();

sub get_TOPPs {
    my ($mecca, $kpp, $mechanism) = @_;
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
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production{$species}{$parent} += $rate(1:$NTIME-2);
            } else {
                my $reaction_string = $kpp->reaction_string($reaction);
                $production{$species}{$reaction_string} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $consumption{$species}{$parent} += $rate(1:$NTIME-2);
            } else {
                my $reaction_string = $kpp->reaction_string($reaction);
                $consumption{$species}{$reaction_string} += $rate(1:$NTIME-2);
            }
        }
    }

    remove_common_processes($production{"HO2x"}, $consumption{"HO2x"});
    my $total_ho2x_production = zeroes(PDL::float, $NTIME-2);
    $total_ho2x_production += $production{'HO2x'}{$_} foreach (keys %{$production{'HO2x'}});

    if ($mechanism =~ /CRI/) {
        foreach my $reaction (keys %{$production{'HO2x'}}) {
            $production{"Ox_$mechanism"}{$reaction} += $production{"Ox_$mechanism"}{'HO2 + NO = OH + NO2'} * $production{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption{"Ox_$mechanism"}{$reaction} += $consumption{"Ox_$mechanism"}{'HO2 + O3 = OH'} * $consumption{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption{"Ox_$mechanism"}{$reaction} += $consumption{"Ox_$mechanism"}{'HO2 + NO3 = OH + NO2'} * $consumption{'HO2x'}{$reaction} / $total_ho2x_production;
        }
        delete $production{"Ox_$mechanism"}{'HO2 + NO = OH + NO2'};
        delete $consumption{"Ox_$mechanism"}{'HO2 + O3 = OH'};
        delete $consumption{"Ox_$mechanism"}{'HO2 + NO3 = OH + NO2'};
    } else {
        foreach my $reaction (keys %{$production{'HO2x'}}) {
            $production{"Ox_$mechanism"}{$reaction} += $production{"Ox_$mechanism"}{'HO2 + NO = NO2 + OH'} * $production{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption{"Ox_$mechanism"}{$reaction} += $consumption{"Ox_$mechanism"}{'HO2 + O3 = OH'} * $consumption{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption{"Ox_$mechanism"}{$reaction} += $consumption{"Ox_$mechanism"}{'HO2 + NO3 = NO2 + OH'} * $consumption{'HO2x'}{$reaction} / $total_ho2x_production;
        }
        delete $production{"Ox_$mechanism"}{'HO2 + NO = NO2 + OH'};
        delete $consumption{"Ox_$mechanism"}{'HO2 + O3 = OH'};
        delete $consumption{"Ox_$mechanism"}{'HO2 + NO3 = NO2 + OH'};
    }
    remove_common_processes($production{"Ox_$mechanism"}, $consumption{"Ox_$mechanism"});

    my %emissions;
    if ($mechanism eq "CB05") {
        foreach my $VOC (sort keys %{$production{"Ox_$mechanism"}}) {
            if ($VOC eq "CH4") {
                my $source = $mecca->balance($VOC);
                $emissions{$VOC} = $source->sum * $dt;
            } else {
                my $primary = cb05_allocations($VOC);
                foreach my $species (@$primary) {
                    my $name = "${species}_$VOC"; 
                    my $emission_reaction = $kpp->producing_from($name, "UNITY");
                    next if (@$emission_reaction == 0);
                    my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
                    my $emission_rate = $mecca->rate($reaction_number);
                    $emission_rate = $emission_rate(1:$NTIME-2);
                    if (defined $emission_rate) {
                        if ($VOC =~ /C3H8/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 1.5; #C3H8 => 1.5 PAR
                        } elsif ($VOC =~ /NC4H10/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 4; #NC4H10 => 4 PAR
                        } elsif ($VOC =~ /IC4H10/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 4; #IC4H10 => 4 PAR
                        } elsif ($VOC =~ /NC5H12/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 5; #NC5H12 => 5 PAR
                        } elsif ($VOC =~ /IC5H12/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 5; #IC5H12 => 5 PAR
                        } elsif ($VOC =~ /NC6H14/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 6; #NC6H14 => 6 PAR
                        } elsif ($VOC =~ /NC7H16/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 7; #NC7H16 => 7 PAR
                        } elsif ($VOC =~ /NC8H18/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 8; #NC8H18 => 8 PAR
                        } elsif ($VOC =~ /BUT1ENE/ and $VOC eq 'PAR') {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 2; #BUT1ENE => 2 PAR
                        } elsif ($VOC =~ /MEPROPENE/ and $VOC eq 'PAR') {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 3; #MEPROPENE => 3 PAR
                        } else {
                            $emissions{$VOC} += $emission_rate->sum * $dt;
                        }
                    }
                }
            }
        }
    } elsif ($mechanism eq "CBM-IV") {
        foreach my $VOC (sort keys %{$production{"Ox_$mechanism"}}) {
            if ($VOC eq "CH4") {
                my $source = $mecca->balance($VOC);
                $emissions{$VOC} = $source->sum * $dt;
            } else {
                my $primary = cbm4_allocations($VOC);
                foreach my $species (@$primary) {
                    my $name = "${species}_$VOC"; 
                    my $emission_reaction = $kpp->producing_from($name, "UNITY");
                    next if (@$emission_reaction == 0);
                    my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
                    my $emission_rate = $mecca->rate($reaction_number);
                    $emission_rate = $emission_rate(1:$NTIME-2);
                    if (defined $emission_rate) {
                        if ($VOC =~ /C3H8/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 1.5; #C3H8 => 1.5 PAR
                        } elsif ($VOC =~ /C2H6/) { 
                            $emissions{$VOC} += $emission_rate->sum * $dt / 0.4; #C2H6 => 0.4 PAR
                        } elsif ($VOC =~ /NC4H10/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 4; #NC4H10 => 4 PAR
                        } elsif ($VOC =~ /IC4H10/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 4; #IC4H10 => 4 PAR
                        } elsif ($VOC =~ /NC5H12/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 5; #NC5H12 => 5 PAR
                        } elsif ($VOC =~ /IC5H12/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 5; #IC5H12 => 5 PAR
                        } elsif ($VOC =~ /NC6H14/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 6; #NC6H14 => 6 PAR
                        } elsif ($VOC =~ /NC7H16/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 7; #NC7H16 => 7 PAR
                        } elsif ($VOC =~ /NC8H18/) {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 8; #NC8H18 => 8 PAR
                        } elsif ($VOC =~ /BUT1ENE/ and $VOC eq 'PAR') {
                            $emissions{$VOC} += $emission_rate->sum * $dt / 2; #BUT1ENE => 2 PAR
                        } else {
                            $emissions{$VOC} += $emission_rate->sum * $dt;
                        }
                    }
                }
            }
        }
    } else {
        foreach my $VOC (sort keys %{$production{"Ox_$mechanism"}}) {
            if ($VOC eq "CH4") {
                my $source = $mecca->balance($VOC);
                $emissions{$VOC} = $source->sum * $dt;
            } else {
                my $emission_reaction = $kpp->producing_from($VOC, "UNITY");
                next if (@$emission_reaction == 0);
                my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
                my $emission_rate = $mecca->rate($reaction_number);
                $emission_rate = $emission_rate(1:$NTIME-2);
                $emissions{$VOC} = $emission_rate->sum * $dt;
            }
        }
    }

    my %TOPP;
    foreach my $VOC (sort keys %{$production{"Ox_$mechanism"}}) {
        next if ( $VOC =~ /=/);
        my $rate = $production{"Ox_$mechanism"}{$VOC}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $production = $rate->sumover * $dt;
        $TOPP{$VOC} = $production(0:13:2) / $emissions{$VOC};
    }

    foreach my $VOC (sort keys %TOPP) {
        if ($mechanism eq "RADM2") { 
            if ($VOC eq "HC3") {
                $TOPP{"C3H8"} = $TOPP{$VOC} * 3 / 2.9;
                $TOPP{"NC4H10"} = $TOPP{$VOC} * 4 / 2.9;
                $TOPP{"IC4H10"} = $TOPP{$VOC} * 4 / 2.9;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC5") {
                $TOPP{"NC5H12"} = $TOPP{$VOC} * 5 / 4.8;
                $TOPP{"IC5H12"} = $TOPP{$VOC} * 5 / 4.8;
                $TOPP{"NC6H14"} = $TOPP{$VOC} * 6 / 4.8;
                $TOPP{"NC7H16"} = $TOPP{$VOC} * 7 / 4.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC8") {
                $TOPP{"NC8H18"} = $TOPP{$VOC} * 8 / 7.9;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLT") { 
                $TOPP{"C3H6"} = $TOPP{$VOC} * 3 / 3.8;
                $TOPP{"BUT1ENE"} = $TOPP{$VOC} * 4 / 3.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLI") {
                $TOPP{"MEPROPENE"} = $TOPP{$VOC} * 4 / 4.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "TOL") {
                $TOPP{"BENZENE"} = $TOPP{$VOC} * 6 / 7.1;
                $TOPP{"TOLUENE"} = $TOPP{$VOC} * 7 / 7.1;
                $TOPP{"EBENZ"} = $TOPP{$VOC} * 8 / 7.1;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "XYL") {
                $TOPP{"MXYL"} = $TOPP{$VOC} * 8 / 8.9;
                $TOPP{"OXYL"} = $TOPP{$VOC} * 8 / 8.9;
                $TOPP{"PXYL"} = $TOPP{$VOC} * 8 / 8.9;
                delete $TOPP{$VOC};
            }
        } elsif ($mechanism eq "RACM") {
            if ($VOC eq "HC3") {
                $TOPP{"C3H8"} = $TOPP{$VOC} * 3 / 2.9;
                $TOPP{"NC4H10"} = $TOPP{$VOC} * 4 / 2.9;
                $TOPP{"IC4H10"} = $TOPP{$VOC} * 4 / 2.9;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC5") {
                $TOPP{"NC5H12"} = $TOPP{$VOC} * 5 / 4.8;
                $TOPP{"IC5H12"} = $TOPP{$VOC} * 5 / 4.8;
                $TOPP{"NC6H14"} = $TOPP{$VOC} * 6 / 4.8;
                $TOPP{"NC7H16"} = $TOPP{$VOC} * 7 / 4.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC8") {
                $TOPP{"NC8H18"} = $TOPP{$VOC} * 8 / 7.9;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLT") { 
                $TOPP{"C3H6"} = $TOPP{$VOC} * 3 / 3.8;
                $TOPP{"BUT1ENE"} = $TOPP{$VOC} * 4 / 3.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLI") {
                $TOPP{"MEPROPENE"} = $TOPP{$VOC} * 4 / 5;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "TOL") {
                $TOPP{"BENZENE"} = $TOPP{$VOC} * 6 / 7.1;
                $TOPP{"TOLUENE"} = $TOPP{$VOC} * 7 / 7.1;
                $TOPP{"EBENZ"} = $TOPP{$VOC} * 8 / 7.1;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "XYL") {
                $TOPP{"MXYL"} = $TOPP{$VOC} * 8 / 8.9;
                $TOPP{"OXYL"} = $TOPP{$VOC} * 8 / 8.9;
                $TOPP{"PXYL"} = $TOPP{$VOC} * 8 / 8.9;
                delete $TOPP{$VOC};
            }
        } elsif ($mechanism eq "RACM2") {
            if ($VOC eq "HC3") {
                $TOPP{"C3H8"} = $TOPP{$VOC} * 3 / 3.6;
                $TOPP{"NC4H10"} = $TOPP{$VOC} * 4 / 3.6;
                $TOPP{"IC4H10"} = $TOPP{$VOC} * 4 / 3.6;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC5") {
                $TOPP{"NC5H12"} = $TOPP{$VOC} * 5 / 5.6;
                $TOPP{"IC5H12"} = $TOPP{$VOC} * 5 / 5.6;
                $TOPP{"NC6H14"} = $TOPP{$VOC} * 6 / 5.6;
                $TOPP{"NC7H16"} = $TOPP{$VOC} * 7 / 5.6;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "HC8") {
                $TOPP{"NC8H18"} = $TOPP{$VOC} * 8 / 7.9;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLT") { 
                $TOPP{"C3H6"} = $TOPP{$VOC} * 3 / 3.8;
                $TOPP{"BUT1ENE"} = $TOPP{$VOC} * 4 / 3.8;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "OLI") {
                $TOPP{"MEPROPENE"} = $TOPP{$VOC} * 4 / 5;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "TOL") {
                $TOPP{"TOLUENE"} = $TOPP{$VOC} * 7 / 7.1;
                $TOPP{"EBENZ"} = $TOPP{$VOC} * 8 / 7.1;
                delete $TOPP{$VOC};
            }
        } elsif ($mechanism eq "MOZART-4") {
            if ($VOC eq "BIGALK") {
                $TOPP{"NC4H10"} = $TOPP{$VOC} * 4 / 5;
                $TOPP{"IC4H10"} = $TOPP{$VOC} * 4 / 5;
                $TOPP{"NC5H12"} = $TOPP{$VOC} * 5 / 5;
                $TOPP{"IC5H12"} = $TOPP{$VOC} * 5 / 5;
                $TOPP{"NC6H14"} = $TOPP{$VOC} * 6 / 5;
                $TOPP{"NC7H16"} = $TOPP{$VOC} * 7 / 5;
                $TOPP{"NC8H18"} = $TOPP{$VOC} * 8 / 5;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "BIGENE") {
                $TOPP{"BUT1ENE"} = $TOPP{$VOC} * 4 / 4;
                $TOPP{"MEPROPENE"} = $TOPP{$VOC} * 4 / 4;
                delete $TOPP{$VOC};
            } elsif ($VOC eq "TOLUENE") {
                $TOPP{"BENZENE"} = $TOPP{$VOC} * 6 / 7;
                $TOPP{"EBENZ"} = $TOPP{$VOC} * 8 / 7;
                $TOPP{"MXYL"} = $TOPP{$VOC} * 8 / 7;
                $TOPP{"OXYL"} = $TOPP{$VOC} * 8 / 7;
                $TOPP{"PXYL"} = $TOPP{$VOC} * 8 / 7;
            }
        }
    }
    return \%TOPP;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
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

sub cb05_allocations { #each parent VOC expressed as CB05 allocation
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

sub cbm4_allocations { #each parent VOC expressed as CBM-IV allocation
    my ($parent) = @_;
    my @allocations = ();
    
    if ($parent eq "C2H6" or $parent eq "C3H8" or $parent eq "NC4H10" or $parent eq "IC4H10" or $parent eq "NC5H12" or $parent eq "IC5H12" or $parent eq "NC6H14" or $parent eq "NC7H16" or $parent eq "NC8H18" or $parent eq "BENZENE"){
        @allocations = qw( PAR );
    } elsif ($parent eq "CH4" ){
        @allocations = qw( CH4 );
    } elsif ($parent eq "C3H6" or $parent eq "BUT1ENE") {
        @allocations = qw( OLE PAR );
    } elsif ($parent eq "MEPROPENE") {
        @allocations = qw( PAR HCHO ALD2 );
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

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane ';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane ';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane ';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane ';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane ';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane ';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane ';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane ";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane ";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene ';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene ';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene ";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene ';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene ";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene ";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene ';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene ";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene ';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene ";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene ";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
