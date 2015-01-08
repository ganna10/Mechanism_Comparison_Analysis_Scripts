#!/usr/bin/env perl 
# First day total Ox production budget plots from all tagged mechanisms, lumped species are de-lumped into constituent VOC compared to MCM v3.2
# Version 0: Jane Coates 8/1/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base_dir = "/local/home/coates/MECCA/";
my $mecca = MECCA->new("$base_dir/CB05_tagged/boxmodel"); 
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my (%data, %families, %weights);
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 ); 
#my @mechanisms = qw( MCMv3.2 RACM2 );

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base_dir/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base_dir/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base_dir/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]; 
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(dplyr) `,
);

$R->set('VOC', [sort keys %{$data{"MCMv3.2"}}]);
$R->run(q` data = data.frame(VOC) `);
foreach my $mechanism (keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` ox = c() `);
    foreach my $VOC (sort keys %{$data{"MCMv3.2"}}) {
        $R->set('value', $data{$mechanism}{$VOC}->at(0));
        $R->run(q` ox = c(ox, value) `);
    }
    $R->run(q` data[mechanism] = ox `);
}
$R->run(q` data = gather(data, Mechanism, Ox, -VOC, -MCMv3.2) `);
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c( "Inorganic" = "blue", "CO" = "#000000", "Methane" = "#ff0000", "Ethane" = "#696537", "Propane" = "#f9c600", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb", "Ethene" = "#86b650", "Propene" = "#6c254f", "Butene" = "#ee6738", "2-Methylpropene" = "#58691b", "Isoprene" = "#8ed6d5", "Benzene" = "#f3aa7f", "Toluene" = "#c65d6c", "m-Xylene" = "#888a87", "o-Xylene" = "#0e5c28", "p-Xylene" = "#b569b3", "Ethylbenzene" = "#2c9def" ) `,
    #q` data$VOC = factor(data$VOC, levels = c("CO", "Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene , "Ethylbenzene")) `,
        q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `,
);

$R->run(q` plot = ggplot(data = data, aes(x = MCMv3.2, y = Ox, colour = VOC, group = VOC)) `,
        q` plot = plot + geom_point(shape = 19) `,
        q` plot = plot + facet_wrap( ~ Mechanism, ncol = 2) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + xlab("MCM v3.2 Ox Production") `,
        q` plot = plot + ylab("Ox Production") `,
        q` plot = plot + geom_abline(intercept = 0, slope = 1) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + scale_colour_manual(values = my.colours, guide = guide_legend(nrow = 3)) `,
);

$R->run(q` CairoPDF(file = "first_day_Ox_production.pdf", width = 8, height = 11.3) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
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
            print "No family for $species\n";
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
                my $string = $kpp->reaction_string($reaction);
                $string = "CO " if ($string =~ /CO \+ OH|OH \+ CO/);
                $production{$species}{$string} += $rate(1:$NTIME-2);
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
                my $string = $kpp->reaction_string($reaction);
                $string = "CO " if ($string =~ /CO \+ OH|OH \+ CO/);
                $consumption{$species}{$string} += $rate(1:$NTIME-2);
            }
        }
    }
    remove_common_processes($production{"HO2x"}, $consumption{"HO2x"});
    my $total_ho2x_production = zeroes(PDL::float, $NTIME-2);
    $total_ho2x_production += $production{"HO2x"}{$_} foreach (keys %{$production{"HO2x"}});
    
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

    foreach my $process (keys %{$production{"Ox_$mechanism"}}) {
        my $reshape = $production{"Ox_$mechanism"}{$process}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{"Ox_$mechanism"}{$process} = $integrate;
    }

    #de-lump lumped VOC to MCM species
    if ($mechanism eq "RADM2" or $mechanism eq "RACM") {
        $production{"Ox_$mechanism"}{"C3H8"} = 0.628 * $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"NC4H10"} = 0.243 * $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"IC4H10"} = 0.129 * $production{"Ox_$mechanism"}{"HC3"};
        delete $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"NC5H12"} = 0.264 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"IC5H12"} = 0.615 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"NC6H14"} = 0.086 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"NC7H16"} = 0.035 * $production{"Ox_$mechanism"}{"HC5"};
        delete $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"C3H6"} = 0.875 * $production{"Ox_$mechanism"}{"OLT"}; 
        $production{"Ox_$mechanism"}{"BUT1ENE"} = 0.125 * $production{"Ox_$mechanism"}{"OLT"}; 
        delete $production{"Ox_$mechanism"}{"OLT"};
        $production{"Ox_$mechanism"}{"BENZENE"} = 0.232 * $production{"Ox_$mechanism"}{"TOL"};
        $production{"Ox_$mechanism"}{"TOLUENE"} = 0.667 * $production{"Ox_$mechanism"}{"TOL"};
        $production{"Ox_$mechanism"}{"EBENZ"} = 0.101 * $production{"Ox_$mechanism"}{"TOL"};
        delete $production{"Ox_$mechanism"}{"TOL"};
        $production{"Ox_$mechanism"}{"MXYL"} = 0.5 * $production{"Ox_$mechanism"}{"XYL"};
        $production{"Ox_$mechanism"}{"OXYL"} = 0.244 * $production{"Ox_$mechanism"}{"XYL"};
        $production{"Ox_$mechanism"}{"PXYL"} = 0.256 * $production{"Ox_$mechanism"}{"XYL"};
        delete $production{"Ox_$mechanism"}{"XYL"};
    } elsif ($mechanism eq "RACM2") {
        $production{"Ox_$mechanism"}{"C3H8"} = 0.628 * $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"NC4H10"} = 0.243 * $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"IC4H10"} = 0.129 * $production{"Ox_$mechanism"}{"HC3"};
        delete $production{"Ox_$mechanism"}{"HC3"};
        $production{"Ox_$mechanism"}{"NC5H12"} = 0.264 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"IC5H12"} = 0.615 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"NC6H14"} = 0.086 * $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"NC7H16"} = 0.035 * $production{"Ox_$mechanism"}{"HC5"};
        delete $production{"Ox_$mechanism"}{"HC5"};
        $production{"Ox_$mechanism"}{"C3H6"} = 0.875 * $production{"Ox_$mechanism"}{"OLT"}; 
        $production{"Ox_$mechanism"}{"BUT1ENE"} = 0.125 * $production{"Ox_$mechanism"}{"OLT"}; 
        delete $production{"Ox_$mechanism"}{"OLT"};
        $production{"Ox_$mechanism"}{"BENZENE"} = 0.232 * $production{"Ox_$mechanism"}{"TOL"};
        $production{"Ox_$mechanism"}{"TOLUENE"} = 0.667 * $production{"Ox_$mechanism"}{"TOL"};
        $production{"Ox_$mechanism"}{"EBENZ"} = 0.101 * $production{"Ox_$mechanism"}{"TOL"};
        delete $production{"Ox_$mechanism"}{"TOL"};
    } elsif ($mechanism =~ /MOZ/) {
        $production{"Ox_$mechanism"}{"BENZENE"} = 0.166 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"TOLUENE_MOZART"} = 0.478 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"MXYL"} = 0.142 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"OXYL"} = 0.069 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"PXYL"} = 0.073 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"EBENZ"} = 0.073 * $production{"Ox_$mechanism"}{"TOLUENE"};
        $production{"Ox_$mechanism"}{"TOLUENE"} = $production{"Ox_$mechanism"}{"TOLUENE_MOZART"};
        delete $production{"Ox_$mechanism"}{"TOLUENE_MOZART"}; 
        $production{"Ox_$mechanism"}{"NC4H10"} = 0.285 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"IC4H10"} = 0.151 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"NC5H12"} = 0.146 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"IC5H12"} = 0.340 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"NC6H14"} = 0.048 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"NC7H16"} = 0.020 * $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"NC8H18"} = 0.010 * $production{"Ox_$mechanism"}{"BIGALK"};
        delete $production{"Ox_$mechanism"}{"BIGALK"};
        $production{"Ox_$mechanism"}{"BUT1ENE"} = 0.333 * $production{"Ox_$mechanism"}{"BIGENE"};
        $production{"Ox_$mechanism"}{"MEPROPENE"} = 0.667 * $production{"Ox_$mechanism"}{"BIGENE"};
        delete $production{"Ox_$mechanism"}{"BIGENE"};
    }

    foreach my $process (sort keys %{$production{"Ox_$mechanism"}}) {
        if ($process =~ / = /) {
            $production{"Ox_$mechanism"}{"Inorganic"} += $production{"Ox_$mechanism"}{$process};
            delete $production{"Ox_$mechanism"}{$process};
        } else {
            my $name = get_chemical_name($process);
            $production{"Ox_$mechanism"}{$name} = $production{"Ox_$mechanism"}{$process};
            delete $production{"Ox_$mechanism"}{$process};
        }
    }
    return $production{"Ox_$mechanism"};
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

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'CO ') {
        $chemical_species = 'CO'
    } elsif ($VOC eq 'CH4') {
        $chemical_species = 'Methane';
    } elsif ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene";
    } elsif ($VOC eq 'Others') {
        $chemical_species = 'Others';
    } else {
        print "No chemical species found for $VOC\n";
    }
}
