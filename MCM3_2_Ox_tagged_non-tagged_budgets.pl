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
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(tidyr) `);
$R->run(q` library(Cairo) `);
$R->run(q` library(grid) `);
$R->run(q` library(gridExtra) `);

$R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
        q` my.colours = c( "Production Others" = "#696537", 
                           "C2H5O" = "#f9c600", 
                           "C2H5O2 + NO" = "#76afca", 
                           "HCHO + hv" = "#dc3522", 
                           "CH3CO3 + NO2" = "#8c6238", 
                           "HCHO + OH" = "#9bb08f", 
                           "CH3O" = "#e7e85e", 
                           "CH3O2 + NO" = "#8b1537", 
                           "CO + OH" = "#f9c500", 
                           "O3 + OH" = "#603912", 
                           "CH3CO3 + NO" = "#b569b3",
                           "NC7H16" = "#f9c600", 
                           "EBENZ" = "#76afca", 
                           "BENZENE" = "#dc3522", 
                           "OXYL" = "#f3aa7f", 
                           "PXYL" = "#be2448", 
                           "NC6H14" = "#86b650", 
                           "IC4H10" = "#e7e85e", 
                           "C3H6" = "#2b9eb3", 
                           "C2H6" = "#0352cb", 
                           "MXYL" = "#6db875", 
                           "C5H8" = "#0c3f78", 
                           "C2H4" = "#cc6329", 
                           "NC5H12" = "#b569b3", 
                           "C3H8" = "#0e5c28", 
                           "TOLUENE" = "#8c1531", 
                           "NC4H10" = "#ef6638", 
                           "IC5H12" = "#4c9383", 
                           "CH4" = "#6c254f") `,
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
                        "CO + OH" = "CO",
                        "Production Others" = "Production Others" ) `, 
);
$R->run(q` plotting = function (data, title) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Process)) ;
                                                plot = plot + geom_bar(stat = "identity") ;
                                                plot = plot + theme_bw() ;
                                                plot = plot + ggtitle(title) ;
                                                plot = plot + theme(plot.title = element_text(face = "bold")) ;
                                                plot = plot + ylab("Reaction Rate (molecules cm-3 s-1)") ;
                                                plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.7)) ;
                                                plot = plot + theme(axis.title = element_text(face = "bold")) ;
                                                plot = plot + theme(axis.title.x = element_blank()) ;
                                                plot = plot + theme(panel.grid = element_blank()) ;
                                                plot = plot + theme(legend.title = element_blank()) ;
                                                plot = plot + theme(legend.position = c(1.031, 1.031), legend.justification = c(1.031, 1.031)) ;
                                                plot = plot + scale_fill_manual(limits = rev(levels(data$Process)), labels = my.names, values = my.colours) ;
                                                plot = plot + scale_y_continuous(limits=c(0, 1.5e9), breaks=seq(0, 1.5e9, 2e8), expand = c(0, 1e7));
                                                plot = plot + theme(legend.key = element_blank()) ;
                                                plot = plot + theme(panel.border = element_rect(colour = "black")) ;
                                                plot = plot + theme(plot.margin = unit(c(0, 0, 0, -0.04), "cm")) ;
                                                return(plot) } `,
);

$R->set('Time', [("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")]);
$R->run(q` plots = list() `);
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
        $R->set('run', "(b) Not Tagged");
    } else {
        $R->set('run', "(a) Tagged");
    }
    $R->run(q` pre = gather(pre, Process, Rate, -Time) `,
            q` plot = plotting(pre, run) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(plot) `);
#print $p, "\n";
$R->run(q` CairoPDF(file = "MCMv3_2_tagged_non_tagged_Ox_budget.pdf", width = 8.5, height = 5.8) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[2]] ,
                                                    plots[[1]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1) `,
        q` print(multiplot) `,
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
            if (defined $parent) { # for tagged reactions
                $production{$species}{$parent} += $rate(1:$NTIME-2);
            } else {
                my $string = $kpp->reaction_string($reaction);
                my ($reactants, $products) = split / = /, $string;
                $production{$species}{$reactants} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
            }
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            if (defined $parent) { # for tagged reactions
                $consumption{$species}{$parent} += $rate(1:$NTIME-2);
            } else {
                my $string = $kpp->reaction_string($reaction);
                my ($reactants, $products) = split / = /, $string;
                $consumption{$species}{$reactants} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
            }
        }
    } 
    remove_common_processes($production{'HO2x'}, $consumption{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production{'HO2x'}{$_} for (keys %{ $production{'HO2x'} });

    foreach my $reaction( keys %{ $production{'HO2x'} }) {
        $production{$Ox}{$reaction} += $production{$Ox}{'HO2 + NO'} * $production{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption{$Ox}{$reaction} += $consumption{$Ox}{'HO2 + O3'} * $consumption{'HO2x'}{$reaction} / $ho2x_total_production;
        $consumption{$Ox}{$reaction} += $consumption{$Ox}{'HO2 + NO3'} * $consumption{'HO2x'}{$reaction} / $ho2x_total_production;
    }
    delete $production{$Ox}{'HO2 + NO'};
    delete $consumption{$Ox}{'HO2 + O3'};
    delete $consumption{$Ox}{'HO2 + NO3'}; 
    remove_common_processes($production{$Ox}, $consumption{$Ox});

    my $prod_others_max = 1.5e8;
    foreach my $item (keys %{$production{$Ox}}) {
        if ($production{$Ox}{$item}->sum < $prod_others_max) {
            $production{$Ox}{'Production Others'} += $production{$Ox}{$item};
            delete $production{$Ox}{$item};
        }
    }

    foreach my $item (keys %{$production{$Ox}}) {
        my $reshape = $production{$Ox}{$item}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production{$Ox}{$item} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
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
