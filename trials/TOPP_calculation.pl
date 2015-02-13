#! /usr/bin/env perl
# Calculate daily TOPP of all VOCs in each mechanism for constant emission runs
# Version 0: Jane Coates 5/12/2014

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
my $index = 0;
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
    $index++;
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
                                                         plot = plot + theme(legend.key = element_blank());
                                                         plot = plot + theme(panel.border = element_rect(colour = "black"));
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
        $R->set('voc', $VOC);
        $R->set('topp', [ map { $_ } $data{$mechanism}{$VOC}->dog ]);
        $R->run(q` data[voc] = topp `);
    }
    $R->run(q` daily.data = gather(data, VOC, TOPP, -Time) `,
            q` cumulative = as.data.frame(lapply(data[, -1], cumsum)) `,
            q` cumulative$Time = Time `,
            q` cumulative.data = gather(cumulative, VOC, TOPP, -Time) `,
            q` plotting(daily.data, daily.title, daily.filename) `,
            q` plotting(cumulative.data, cumulative.title, cumulative.filename) `,
    );
    #my $p = $R->run(q` print(cumulative) `);
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
        my $rate = $production{"Ox_$mechanism"}{$VOC}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $production = $rate->sumover * $dt;
        $production = $production(0:13:2);
        if (defined $emissions{$VOC}) {
            $TOPP{$VOC} = $production / $emissions{$VOC};
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
