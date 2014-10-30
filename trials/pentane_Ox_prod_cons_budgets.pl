#! /usr/bin/env perl
# pentane Ox production and consumption budget for each mechanism
# Version 0: Jane Coates 24/10/2014
# Version 1: Jane Coates 29/10/2014 changing to integrating using pdl not R

use strict;
use diagnostics;
use KPP;
use MECCA;
use Statistics::R;
use PDL;
use PDL::NiceSlice;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05" );
#my @runs = qw( CBM4_tagging CB05_tagging );
#my @mechanisms = ( "(h) CBM-IV", "(i) CB05" );
my $index = 0;

my (%families, %weights, %plot_data);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $run);
        $plot_data{$mechanisms[$index]} = get_data($kpp, $mecca, $mechanisms[$index], $mech_species);
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(plyr) `,
        q` library(grid) `,
        q` library(dplyr) `,
);

my @days = ("Day1", "Day2", "Day3", "Day4", "Day5", "Day6", "DAy7");
$R->set('time', [@days]);
$R->run(q` data = data.frame() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` pre = data.frame(time) `);
    foreach my $reaction (sort keys %{$plot_data{$run}}) {
        $R->set('reaction', $reaction);
        $R->set('rate', [map { $_ } $plot_data{$run}{$reaction}->dog]);
        $R->run(q` pre[reaction] = rate `);
    }
    $R->set('mechanism', $run);
    $R->run(q` pre$Mechanism = rep(mechanism, length(time)) `,
            q` pre$VOC = rep("Pentane", length(time)) `,
            q` pre = melt(pre, id.vars = c("time", "Mechanism", "VOC"), variable.name = "Reaction", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(data = subset(data, Rate > 0), stat = "identity") `, 
        q` plot = plot + geom_bar(data = subset(data, Rate < 0), stat = "identity") `, 
        q` plot = plot + facet_grid(time ~ VOC) `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05"))) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text.x = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.text.y = element_text(size = 200, face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 200)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 160)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 160)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(10, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 170)) `,
);
my $p = $R->run(q` print(data) `);
print $p, "\n";

$R->run(q` CairoPDF(file = "pentane_Ox_intermediates.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC) = @_;
    $families{"HO2x_${mechanism}"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_${mechanism}", "HO2x_${mechanism}");

    my ($producers, $producer_yields, %production_reaction_rates, $consumers, $consumer_yields, %consumption_reaction_rates);
    foreach my $species (@loop) {
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
        
        my $prod_max = 5e9;
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string;
            if ($reactants =~ /XO2/ and $species =~ /RACM|CB|RADM/) {
                my $operator = "XO2_" . $VOC;
                my $op_producers = $kpp->producing($operator);
                my $op_producer_yields = $kpp->effect_on($operator, $op_producers); 
                die "No producers found for $operator\n" if (@$op_producers == 0);

                for (0..$#$op_producers) { #get rates for all producing reactions
                    my $reaction = $op_producers->[$_];
                    my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
                    my $reaction_number = $kpp->reaction_number($reaction);
                    my $rate = $op_producer_yields->[$_] * $mecca->rate($reaction_number); 
                    next if ($rate->sum == 0); # do not include reactions that do not occur 
                    my $reaction_string = $kpp->reaction_string($reaction);
                    $reaction_string =~ s/_(.*?)\b//;
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
                } 
            } else {
                $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_$VOC\b//g;
            my($reactants, $products) = split / = /, $reaction_string; 
            $consumption_reaction_rates{$reactants} += $rate(1:$NTIME-2); 
        }
    
        remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
        foreach my $reaction (sort keys %production_reaction_rates) {
            if ($production_reaction_rates{$reaction}->sum < $prod_max) {
                $production_reaction_rates{"Production Others"} += $production_reaction_rates{$reaction};
                delete $production_reaction_rates{$reaction};
            }
        }

        foreach my $reaction (sort keys %consumption_reaction_rates) {
            if ($consumption_reaction_rates{$reaction}->sum > -$prod_max) {
                $consumption_reaction_rates{"Consumption Others"} += $consumption_reaction_rates{$reaction};
                delete $consumption_reaction_rates{$reaction};
            }
        }
    }

    my $dt = $mecca->dt->at(0); #model time step
    my $parent_emissions;
    if ($mechanism =~ /CB/) {
        my $name = "PAR_NC5H12";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
    } else {
        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    $production_reaction_rates{$_} = $production_reaction_rates{$_} * $dt / $parent_emissions foreach (sort keys %production_reaction_rates);
    $consumption_reaction_rates{$_} = $consumption_reaction_rates{$_} * $dt / $parent_emissions foreach (sort keys %consumption_reaction_rates);
    
    my %plot_data;
    my $n_per_day = 43200 / $dt;
    my $n_days = int ($NTIME / $n_per_day);
    foreach (keys %production_reaction_rates) {
        my $reshape = $production_reaction_rates{$_}->copy->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $plot_data{$_} = $integrate;
    }
    foreach (keys %consumption_reaction_rates) {
        my $reshape = $consumption_reaction_rates{$_}->copy->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $plot_data{$_} = $integrate;
    }
    return \%plot_data;
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

sub get_model_name {
    my ($VOC, $mechanism) = @_;
    my $species;
    if ($VOC eq "Pentane") {
        if ($mechanism  =~ /MCM|CRI|CB/) {
            $species = "NC5H12";
        } elsif ($mechanism =~ /MOZART/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RADM|RACM/) {
            $species = "HC5";
        } else {
            print "No mechanism for $VOC\n";
        }
    } elsif ($VOC eq "Toluene") {
        if ($mechanism =~ /MCM|CRI|MOZART|CB/) {
            $species = "TOLUENE";
        } elsif ($mechanism =~ /RADM|RACM/) {
            $species = "TOL";
        } else {
            print "No mechanism for $VOC\n";
        }
    } else {
        print "No species found for $VOC\n";
    }
}
