#!/usr/bin/env perl 
# Ox production budget plots from all tagged mechanisms
# Version 0: Jane Coates 28/08/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $mcm_3_2_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mcm_3_2_mecca = MECCA->new($mcm_3_2_run); 
my $NTIME = $mcm_3_2_mecca->time->nelem;
my $times = $mcm_3_2_mecca->time;
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog;
my @time_blocks;
foreach my $time (@time_axis) {
    if ($time <= 12) {
        push @time_blocks, "Day 1";
    } elsif ($time > 12 and $time <= 24) {
        push @time_blocks, "Night 1";
    } elsif ($time > 24 and $time <= 36) {
        push @time_blocks, "Day 2";
    } elsif ($time > 36 and $time <= 48) {
        push @time_blocks, "Night 2";
    } elsif ($time > 48 and $time <= 60) {
        push @time_blocks, "Day 3",
    } elsif ($time > 60 and $time <= 72) {
        push @time_blocks, "Night 3";
    } elsif ($time > 72 and $time <= 84) {
        push @time_blocks, "Day 4";
    } elsif ($time > 84 and $time <= 96) {
        push @time_blocks, "Night 4";
    } elsif ($time > 96 and $time <= 108) {
        push @time_blocks, "Day 5";
    } elsif ($time > 108 and $time <= 120) {
        push @time_blocks, "Night 5";
    } elsif ($time > 120 and $time <= 132) {
        push @time_blocks, "Day 6";
    } elsif ($time > 132 and $time <= 144) {
        push @time_blocks, "Night 6";
    } elsif ($time > 144 and $time <= 156) {
        push @time_blocks, "Day 7";
    } else {
        push @time_blocks, "Night 7";
    }
}

my %families = ( 'HO2x' => [ qw( HO2 HO2NO2 ) ] );
my (%weights, %production_rates, %sorted_plot_data); 
my $base = "/work/users/jco/MECCA";
my @model_runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging ); 
my @run_name = qw( mcm_3_2 mcm_3_1 cri mozart radm2 racm racm2 cbm4 cb05 );
my $run_name_index = 0;

foreach my $run (@model_runs) {
    print "$run and $run_name[$run_name_index]\n";
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{$run_name[$run_name_index]} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]; 
    $weights{$run_name[$run_name_index]} = { NO3 => 2, N2O5 => 3 };
    ($production_rates{$run_name[$run_name_index]}) = get_rates($mecca, $kpp, $run_name[$run_name_index]);
    ($sorted_plot_data{$run_name[$run_name_index]}) = get_sorted_data($production_rates{$run_name[$run_name_index]});
    $run_name_index++;
}

my ($reaction_plot_data) = plot(\@time_blocks, \%sorted_plot_data);

sub plot { #create dataframe and then create plot
    my ($time, $plot_data) = @_;
    
    my $R = Statistics::R->new();
    $R->set('time', [@$time]);
    $R->run(q` mcm3.2.data = data.frame(time) `,
            q` mcm3.1.data = data.frame(time) `,
            q` cri.data = data.frame(time) `,
            q` mozart.data = data.frame(time) `,
            q` radm2.data = data.frame(time) `,
            q` racm.data = data.frame(time) `,
            q` racm2.data = data.frame(time) `,
            q` cbm4.data = data.frame(time) `,
            q` cb05.data = data.frame(time) `,
    );
    sub get_chemical_name {
        my ($VOC) = @_;
        my $chemical_species;
        if ($VOC eq 'CO + OH = HO2' or $VOC eq 'OH + CO = HO2') {
            $chemical_species = 'CO ';
        } elsif ($VOC eq 'CH4') {
            $chemical_species = 'Methane ';
        } elsif ($VOC eq 'C2H6' or $VOC eq 'ETH') {
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
        } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
            $chemical_species = 'Ethene ';
        } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
            $chemical_species = 'Propene ';
        } elsif ($VOC eq 'OLI') {
            $chemical_species = '2-Methylpropene ';
        } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
            $chemical_species = "Isoprene ";
        } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
            $chemical_species = 'Toluene ';
        } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
            $chemical_species = "m-Xylene ";
        } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
            $chemical_species = 'o-Xylene ';
        } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
            $chemical_species = "p-Xylene ";
        } elsif ($VOC eq 'Production Others') {
            $chemical_species = 'Production Others';
        } else {
            print "No chemical species found for $VOC\n";
        }
    }
    foreach my $run (sort keys %$plot_data) {
        for my $item (0..$#{$plot_data->{$run}}) {
            foreach my $VOC (sort keys %{$plot_data->{$run}[$item]}) {
                my $chemical = get_chemical_name($VOC);
                $R->set('name', $chemical);
                $R->set('rate', [@{$plot_data->{$run}[$item]{$VOC}}]);
                if ($run =~ /mcm_3_2/) {
                    $R->run(q` mcm3.2.data[name] = rate `);
                } elsif ($run =~ /mcm_3_1/) {
                    $R->run(q` mcm3.1.data[name] = rate `);
                } elsif ($run =~ /cri/) {
                    $R->run(q` cri.data[name] = rate `);
                } elsif ($run =~ /mozart/) {
                    $R->run(q` mozart.data[name] = rate `);
                } elsif ($run =~ /radm2/) {
                    $R->run(q` radm2.data[name] = rate `);
                } elsif ($run =~ /racm2/) {
                    $R->run(q` racm2.data[name] = rate `);
                } elsif ($run =~ /racm/) {
                    $R->run(q` racm.data[name] = rate `);
                } elsif ($run =~ /cbm4/) {
                    $R->run(q` cbm4.data[name] = rate `);
                } elsif ($run =~ /cb05/) {
                    $R->run(q` cb05.data[name] = rate `);
                } 
            }
        }
    }

    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(Cairo) `);
    
    $R->run(q` arrange.data = function(data, mechanism) {  data = ddply(data, .(time), colwise(sum)) ;
                                                           data = data[1:7,] ;
                                                           data = melt(data) ;
                                                           colnames(data) = c("Time", "VOC", "Rate") ;
                                                           data["Mechanism"] = rep(mechanism, length(data$Time)) ;
                                                           return(data) }`
    );

    $R->run(q` mcm3.2.data = arrange.data(mcm3.2.data, "(a) MCM v3.2") `,
            q` mcm3.1.data = arrange.data(mcm3.1.data, "(b) MCM v3.1") `,
            q` cri.data = arrange.data(cri.data, "(c) CRI v2") `,
            q` mozart.data = arrange.data(mozart.data, "(g) MOZART-4") `,
            q` radm2.data = arrange.data(radm2.data, "(d) RADM2") `,
            q` racm.data = arrange.data(racm.data, "(e) RACM") `,
            q` racm2.data = arrange.data(racm2.data, "(f) RACM2") `,
            q` cbm4.data = arrange.data(cbm4.data, "(h) CBM-IV") `,
            q` cb05.data = arrange.data(cb05.data, "(i) CB05") `,
            q` data = rbind(mcm3.2.data, mcm3.1.data, cri.data, mozart.data, radm2.data, racm.data, racm2.data, cbm4.data, cb05.data) `,
            q` VOC.levels = c(  "Methane ", 
                                "CO ", 
                                "Ethane ", 
                                "Propane ", 
                                "2-Methylpropane ", 
                                "Butane ", 
                                "Pentane ", 
                                "2-Methylbutane ", 
                                "Hexane ", 
                                "Ethene ", 
                                "Propene ", 
                                "2-Methylpropene ",
                                "Isoprene ",
                                "Toluene ",
                                "m-Xylene ",
                                "o-Xylene ",
                                "p-Xylene ",
                                "Production Others" ) `,

            #set colours and legend names
            q` my.colours = c(  "Production Others" = "#696537", 
                                "Heptane " = "#f9c600", 
                                "Ethybenzene " = "#76afca", 
                                "Benzene " = "#dc3522",
                                "o-Xylene " = "#8c6238", 
                                "p-Xylene " = "#9bb08f", 
                                "Hexane " = "#8b1537", 
                                "2-Methylpropane " = "#e7e85e", 
                                "Propene " = "#0352cb",
                                "Ethane " = "#86b650",
                                "m-Xylene " = "#6c254f",
                                "Isoprene " = "#ee6738",
                                "Ethene " = "#58691b",
                                "Pentane " = "#8ed6d5",
                                "Propane " = "#f3aa7f",
                                "Toluene " = "#c65d6c",
                                "Butane " = "#888a87", 
                                "2-Methylbutane " = "#0e5c28", 
                                "2-Methylpropene " = "red",
                                "Methane " = "#b569b3", 
                                "CO " = "#2c9def" ) `,
            
            q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis

            #plot 
            q` CairoPDF(file = "Ox_tagged_budget_overlay_all_mechanisms.pdf", width = 141, height = 200) `,
            q` plot = ggplot(data, aes(x = Time, y = Rate, fill = VOC)) `,
            q` plot = plot + geom_bar(stat = "identity") `,
            q` plot = plot + facet_wrap( ~ Mechanism, nrow = 3) `,
            q` plot = plot + theme_bw() `,
            q` plot = plot + theme(panel.margin = unit(5, 'lines')) `,
            q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
            q` plot = plot + theme(strip.background = element_blank()) `,
            q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
            q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
            q` plot = plot + theme(axis.title.x = element_text(size = 120)) `,
            q` plot = plot + ylab(expression(bold(paste(O[x], " Production Rate (molecules ", cm^- 3, s^- 1, ")")))) `,
            q` plot = plot + xlab("\n") `,
            q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
            q` plot = plot + theme(axis.title.y = element_text(size = 200)) `,
            q` plot = plot + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) `,
            q` plot = plot + theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) `,
            q` plot = plot + theme(legend.position = "bottom", legend.title = element_blank()) `,
            q` plot = plot + scale_y_continuous(limits=c(0, 1.8e9), breaks=seq(0, 1.8e9, 2e8), label = scientific_10)`,
            q` plot = plot + theme(legend.key = element_blank()) `,
            q` plot = plot + theme(legend.text = element_text(size = 140)) `,
            q` plot = plot + scale_fill_manual( limits = VOC.levels, values = my.colours, guide = guide_legend(nrow = 2)) `,
            q` print(plot) `,
            q` dev.off() `,
    );

    $R->stop(); 
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

sub get_rates {
    my ($mecca, $kpp, $run) = @_;

    my @lookups = ($run, 'HO2x');
    my (%production_reaction_rates, %species_production_rates, %consumption_reaction_rates, %species_consumption_rates);
    foreach my $species (@lookups) { #get all production and consumption rates
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
            print "No family found for $species\n";
        }

        #check that species reactions are found
        die "No producers found for $species\n" if (@$producers == 0);
        die "No consumers found for $species\n" if (@$consumers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                $string = $kpp->reaction_string($reaction);
                $string =~ s/_$parent//g; #removing tag from reaction strings
                $species_production_rates{$species}{$parent}{$string} += $rate;
                $string = $parent; # in order to merge all production rates from all parent species reactions into 1 pdl
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $production_reaction_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                $string = $kpp->reaction_string($reaction);
                $string =~ s/_$parent//g; #removing tag from reaction strings
                $species_consumption_rates{$species}{$parent}{$string} += $rate;
                $string = $parent; # in order to merge all consumption rates from all parent species reactions into 1 pdl
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $consumption_reaction_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    } 
    remove_common_processes($production_reaction_rates{'HO2x'}, $consumption_reaction_rates{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production_reaction_rates{'HO2x'}{$_} for (keys %{ $production_reaction_rates{'HO2x'} });

    if ($run =~ /cri/) {
        foreach my $reaction( keys %{ $production_reaction_rates{'HO2x'} }) {
            $production_reaction_rates{$run}{$reaction} += $production_reaction_rates{$run}{'HO2 + NO = OH + NO2'} * $production_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
            $consumption_reaction_rates{$run}{$reaction} += $consumption_reaction_rates{$run}{'HO2 + O3 = OH'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
            $consumption_reaction_rates{$run}{$reaction} += $consumption_reaction_rates{$run}{'HO2 + NO3 = OH + NO2'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
        }
        delete $production_reaction_rates{$run}{'HO2 + NO = OH + NO2'};
        delete $consumption_reaction_rates{$run}{'HO2 + O3 = OH'};
        delete $consumption_reaction_rates{$run}{'HO2 + NO3 = OH + NO2'}; 
    } else { 
        foreach my $reaction( keys %{ $production_reaction_rates{'HO2x'} }) {
            $production_reaction_rates{$run}{$reaction} += $production_reaction_rates{$run}{'HO2 + NO = NO2 + OH'} * $production_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
            $consumption_reaction_rates{$run}{$reaction} += $consumption_reaction_rates{$run}{'HO2 + O3 = OH'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
            $consumption_reaction_rates{$run}{$reaction} += $consumption_reaction_rates{$run}{'HO2 + NO3 = NO2 + OH'} * $consumption_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
        }
        delete $production_reaction_rates{$run}{'HO2 + NO = NO2 + OH'};
        delete $consumption_reaction_rates{$run}{'HO2 + O3 = OH'};
        delete $consumption_reaction_rates{$run}{'HO2 + NO3 = NO2 + OH'}; 
    }
    remove_common_processes($production_reaction_rates{$run}, $consumption_reaction_rates{$run});
    return $production_reaction_rates{$run};
}

sub get_sorted_data {
    my ($data) = @_;
    my $prod_others_max = 1e8;

    my $sort_function = sub { $_[0]->sum };
    foreach my $item (keys %$data) { 
        if ($data->{$item}->sum < $prod_others_max) { #get production others
            $data->{'Production Others'} += $data->{$item};
            delete $data->{$item};
        }
    }

    my @sorted_data = sort { &$sort_function($data->{$b}) <=> &$sort_function($data->{$a}) } keys %$data;

    my @final_sorted_data;
    foreach (@sorted_data) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others' or $_ eq 'CH4' or $_ eq "CO + OH = HO2" or $_ eq "OH + CO = HO2");
        push @final_sorted_data, { $_ => $data->{$_} };
    }

    push @final_sorted_data, { 'Production Others' => $data->{'Production Others'} } if (defined $data->{'Production Others'}); 
    unshift @final_sorted_data, { 'CO + OH = HO2' => $data->{'CO + OH = HO2'} } if (defined $data->{'CO + OH = HO2'});
    unshift @final_sorted_data, { 'OH + CO = HO2' => $data->{'OH + CO = HO2'} } if (defined $data->{'OH + CO = HO2'});
    unshift @final_sorted_data, { 'CH4' => $data->{'CH4'} } if (defined $data->{'CH4'});

    my @plot_data;
    foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    } 
    return \@plot_data;
}
