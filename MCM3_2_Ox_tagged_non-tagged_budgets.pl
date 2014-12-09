#!/usr/bin/perl 
# calculate Ox production from tagged and non-tagged model runs
# Version 0: Jane Coates 27/09/2014
# Version 1: Jane Coates 9/12/2014 re-factoring code for constant emissions runs

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_no_tagging/boxmodel"); 
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @runs = qw( no_tagging tagged );
my $index = 0; 
my (%families, %weights, %data);

foreach my $run (@runs) {
    my $boxmodel = "$base/MCMv3.2_$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqnfile = "$base/MCMv3.2_$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base/MCMv3.2_$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{$run} =  [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]; 
    $weights{$run} = { NO3 => 2, N2O5 => 3};
    $data{$run} = get_data($mecca, $kpp, $run);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(tidyr) `);
$R->run(q` library(Cairo) `);

$R->set('Time', [("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")]);
$R->run(q` data = data.frame() `);
foreach my $run (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$run}}) {
        foreach my $item (sort keys %$ref) {
            $R->set('process', $item);
            $R->set('rate', [ map { $_ } $ref->{$item}->dog ]);
            $R->run(q` pre[process] = rate `);
        }
    }
    if ($run =~ /no/) {
        $R->set('run', "Not Tagged");
    } else {
        $R->set('run', "Tagged");
    }
    $R->run(q` pre$Run = rep(run, length(Time)) `,
            q` pre = gather(pre, Process, Rate, -Time, -Run) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
        q` my.colours = c( "Production Others" = "#696537", 
                           "C2H5O = CH3CHO + HO2" = "#f9c600", 
                           "C2H5O2 + NO = C2H5O + NO2" = "#76afca", 
                           "HCHO + hv = CO + HO2 + HO2" = "#dc3522", 
                           "CH3CO3 + NO = CH3O2 + NO2" = "#8c6238", 
                           "HCHO + OH = CO + HO2" = "#9bb08f", 
                           "CH3O = HCHO + HO2" = "#8b1537", 
                           "CH3O2 + NO = CH3O + NO2" = "#e7e85e", 
                           "CO + OH = HO2" = "#2c9def", 
                           "O3 + OH = HO2" = "#603912", 
                           "IC3H7O2 + NO = IC3H7O + NO2" = "#b569b3",
                           "NC7H16" = "#f9c600", 
                           "EBENZ" = "#76afca", 
                           "BENZENE" = "#dc3522", 
                           "OXYL" = "#8c6238", 
                           "PXYL" = "#9bb08f", 
                           "NC6H14" = "#8b1537", 
                           "IC4H10" = "#e7e85e", 
                           "C3H6" = "#0352cb", 
                           "C2H6" = "#86b650", 
                           "MXYL" = "#6c254f", 
                           "C5H8" = "#ee6738", 
                           "C2H4" = "#58691b", 
                           "NC5H12" = "#8ed6d5", 
                           "C3H8" = "#f3aa7f", 
                           "TOLUENE" = "#c65d6c", 
                           "NC4H10" = "#888a87", 
                           "IC5H12" = "#0e5c28", 
                           "CH4" = "#b569b3") `,
        q` my.names = c("NC7H16" = "Heptane", 
                        "EBENZ" = "Ethylbenzene", 
                        "BENZENE" = "Benzene", 
                        "OXYL" = "o-Xylene", 
                        "PXYL" = "p-Xylene", 
                        "NC6H14" = "Hexane", 
                        "IC4H10" = "2-Methylpropane", 
                        "C3H6" = "Propene", 
                        "C2H6" = "Ethane", 
                        "MXYL" = "m-Xylene", 
                        "C5H8" = "Isoprene", 
                        "C2H4" = "Ethene", 
                        "NC5H12" = "Pentane", 
                        "C3H8" = "Propane", 
                        "TOLUENE" = "Toluene", 
                        "NC4H10" = "Butane", 
                        "IC5H12" = "2-Methylbutane", 
                        "CH4" = "Methane", 
                        "CO + OH = HO2" = "CO",
                        "Production Others" = "Production Others", 
                        "C2H5O = CH3CHO + HO2" = "C2H5O = CH3CHO + HO2", 
                        "C2H5O2 + NO = C2H5O + NO2" = "C2H5O2 + NO = C2H5O + NO2", 
                        "HCHO + hv = CO + HO2 + HO2" = "HCHO + hv = CO + HO2 + HO2", 
                        "CH3CO3 + NO = CH3O2 + NO2" = "CH3CO3 + NO = CH3O2 + NO2", 
                        "HCHO + OH = CO + HO2" = "HCHO + OH = CO + HO2", 
                        "CH3O = HCHO + HO2" = "CH3O = HCHO + HO2", 
                        "CH3O2 + NO = CH3O + NO2" = "CH3O2 + NO = CH3O + NO2") `,
);

$R->run(q` plot = ggplot(data, aes( x = Time, y = Rate, fill = Process )) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Run ) `,
        q` plot = plot + theme_bw() `,
);

$R->run(q` CairoPDF(file = "MCMv3.2_tagged_non_tagged_Ox_budget.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop(); 

sub get_data {
    my ($mecca, $kpp, $Ox) = @_;
    $families{'HO2x'} = ( [ qw( HO2 HO2NO2 ) ] );
    my @species = ($Ox, "HO2x");

    my (%production, %consumption); 
    foreach my $species (@species) { #get all production and consumption rates
        my ($producers, $producer_yields, $consumers, $consumer_yields);
        if (exists $families{$species}) { #get family reaction numbers and yields
            $kpp->family({ 
                    name    => $species,
                    members => $families{$species},
                    weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);  
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);  
        } else { #get reaction numbers and yields
            print "No family for $species found\n";
        } 
        die "No producers found for $species\n" if (@$producers == 0);
        die "No consumers found for $species\n" if (@$consumers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string = $kpp->reaction_string($reaction);
            if (defined $parent) { # for tagged reactions
                $string = $parent; # in order to merge all production rates from all parent species reactions into 1 pdl
            }
            $production{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string = $kpp->reaction_string($reaction);
            if (defined $parent) { # for tagged reactions
                $string = $parent; # in order to merge all consumption rates from all parent species reactions into 1 pdl
            }
            $consumption{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    } 
    remove_common_processes($production{'HO2x'}, $consumption{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production{'HO2x'}{$_} for (keys %{ $production{'HO2x'} });

    foreach my $reaction( keys %{ $production{'HO2x'} }) {
        $production{$Ox}{$reaction} += $production{$Ox}{'HO2 + NO = NO2 + OH'} * $production{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption{$Ox}{$reaction} += $consumption{$Ox}{'HO2 + O3 = OH'} * $consumption{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption{$Ox}{$reaction} += $consumption{$Ox}{'HO2 + NO3 = NO2 + OH'} * $consumption{'HO2x'}{$reaction} / $ho2x_total_production;
    }
    delete $production{$Ox}{'HO2 + NO = NO2 + OH'};
    delete $consumption{$Ox}{'HO2 + O3 = OH'};
    delete $consumption{$Ox}{'HO2 + NO3 = NO2 + OH'}; 
    remove_common_processes($production{$Ox}, $consumption{$Ox});

    my $prod_others_max = 8e8;
    my $sort_function = sub { $_[0]->sum };
    foreach my $item (keys %{$production{$Ox}}) {
        if ($production{$Ox}{$item}->sum < $prod_others_max) { #get production others
            $production{$Ox}{'Production Others'} += $production{$Ox}{$item};
            delete $production{$Ox}{$item};
        }
    }

    foreach my $item (keys %{$production{$Ox}}) {
        my $reshape = $production{$Ox}{$item}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production{$Ox}{$item} = $integrate;
    }

    my @sorted_prod = sort { &$sort_function($production{$Ox}{$b}) <=> &$sort_function($production{$Ox}{$a}) } keys %{$production{$Ox}}; 
    my @sorted_plot_data;
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production{$Ox}{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production{$Ox}{'Production Others'} } if (defined $production{$Ox}{'Production Others'}); #add Production Others to the beginning 
    return \@sorted_plot_data;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, "<:encoding(utf-8)", $file or die $!; 
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
