#! /usr/bin/env perl
# OPE analysis for each mechanism
# Version 0: Jane Coates 29/10/2014

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

#my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
#my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my @runs = qw( MOZART_tagging CB05_tagging );
my @mechanisms = ( "MOZART-4", "CB05" );
my $index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    my $carbon_file = "$base/$run/carbons.txt";
    my $carbons = get_carbons($run, $carbon_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $run);
        $plot_data{$mechanisms[$index]}{$VOC} = get_data($kpp, $mecca, $mechanisms[$index], $mech_species, $carbons);
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(plyr) `,
        q` library(grid) `,
        q` library(gridExtra) `,
        q` library(dplyr) `,
        q` library(scales) `,
);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('Time', [@days]);
$R->run(q` data = data.frame() `);
foreach my $run (sort keys %plot_data) {
    foreach my $VOC (sort keys %{$plot_data{$run}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $carbon (sort keys %{$plot_data{$run}{$VOC}} ) {
            $R->set('C.number', $carbon);
            $R->set('OPE', [map { $_ } $plot_data{$run}{$VOC}{$carbon}->dog]);
            $R->run(q` pre[C.number] = OPE `);
        }
        $R->set('Mechanism', $run);
        $R->set('VOC', $VOC);
        $R->run( q` pre$Mechanism = rep(Mechanism, length(pre$Time)) `,
                q` pre$VOC = rep(VOC, length(pre$Time)) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` scientific_10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = OPE, fill = C.number)) `,
        q` plot = plot + geom_bar(stat = "identity") `, 
        q` plot = plot + facet_grid(Time ~ VOC) `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05"))) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + ylab(expression(bold(paste("\nNormalised Ox Production Efficiency (molecules ", (VOC)^-1, cm^3, "s)")))) `,
        #q` plot = plot + scale_y_continuous(labels = scientific_10) `,
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
        #q` plot = plot + theme(legend.position = c(0.10, 0.03)) `,
        #q` plot = plot + theme(legend.justification = c(0.10, 0.03)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(10, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 170)) `,
        #q` plot = plot + guides(fill = guide_legend(direction = "horizontal", nrow = 2)) `,
        #q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

#$R->run(q` CairoPDF(file = "OPE_analysis.pdf", width = 141, height = 200) `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC, $carbons) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_$VOC"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_${mechanism}_$VOC", "HO2x_${mechanism}_$VOC");

    my ($producers, $producer_yields, %production_reaction_rates, $consumers, $consumer_yields, %consumption_reaction_rates, %carbon_production, %carbon_consumption);
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
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants) = $kpp->reactants($reaction);
            $reactants = join ' + ', @$reactants;
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
                    $production_reaction_rates{$reaction} += $rate(1:$NTIME-2);
                } 
            } else {
                $production_reaction_rates{$reaction} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $consumption_reaction_rates{$reaction} += $rate(1:$NTIME-2); 
        }
    
        remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    }

    foreach my $reaction (sort keys %production_reaction_rates) { 
        my ($reactants) = $kpp->reactants($reaction); 
        foreach (@$reactants) {
            if ($_ =~ /_/) {
                my ($lookup, $tag) = split /_/, $_;
                if (defined $carbons{$lookup}) {
                    $carbon_production{"C$carbons{$lookup}"} += $production_reaction_rates{$reaction};
                } elsif ($lookup =~ /\bCO\b/) {
                    $carbon_production{"C1"} += $production_reaction_rates{$reaction};
                } else {
                    print "Nothing found for $lookup\n";
                }
            } elsif ($_ =~ /HC5\b/) {
                my $lookup = "HC5";
                $carbon_production{"C$carbons{$lookup}"} += $production_reaction_rates{$reaction};
            }
        }
    }
    
    foreach my $reaction (sort keys %consumption_reaction_rates) {
        my ($reactants) = $kpp->reactants($reaction);
        foreach (@$reactants) {
            if ($_ =~ /_/) {
                next if ($_ =~ /XO2_/ and $mechanism =~ /RADM2|RACM|CB/);
                my ($lookup, $tag) = split /_/, $_;
                if (defined $carbons{$lookup}) {
                    $carbon_consumption{"C$carbons{$lookup}"} += -$consumption_reaction_rates{$reaction};
                } else {
                    print "Nothing found for $lookup in $mechanism\n";
                }
            }
        }
    }

    my $dt = $mecca->dt->at(0); #model time step

    my $parent_emissions;
    if ($mechanism =~ /CB/ and $VOC eq "NC5H12") {
        my $name = "PAR_NC5H12";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
    } elsif ($mechanism =~ /CB/ and $VOC eq "TOLUENE") {
        my $name = "TOL_TOLUENE";
        my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt ; 
    } else {
        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $production_reaction_rates{$_} /= $parent_emissions foreach (sort keys %production_reaction_rates);
    $consumption_reaction_rates{$_} /= $parent_emissions foreach (sort keys %consumption_reaction_rates);

    my $n_per_day = 43200 / $dt;
    my $n_days = int ($NTIME / $n_per_day);
    my %daily_OPEs;
    foreach my $carbon (sort keys %carbon_production) {
        next unless (exists $carbon_consumption{$carbon} and exists $carbon_production{$carbon});

        my $reshaped_prod = $carbon_production{$carbon}->copy->reshape($n_per_day, $n_days);
        my $integ_prod = $reshaped_prod->sumover;
        my $reshaped_cons = $carbon_consumption{$carbon}->copy->reshape($n_per_day, $n_days);
        my $integ_cons = $reshaped_cons->sumover;

        my $original_carbon_consumption = $carbon_consumption{$carbon};
        $carbon_consumption{$carbon}->where($original_carbon_consumption == 0) += 1; #if consumption rate is 0 then just need production rate
        my $OPE = $carbon_production{$carbon} / $carbon_consumption{$carbon} ;
        my $reshaped_OPE = $OPE->copy->reshape($n_per_day, $n_days);
        my $integrated_OPE = $reshaped_OPE->sumover;
        $integrated_OPE = $integrated_OPE(0:13:2); #choose day time period
        $daily_OPEs{$carbon} = $integrated_OPE;
    }
    return \%daily_OPEs;
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

sub get_carbons {
    my ($run, $file) = @_;
    
    my $carbons;
    if ($run =~ /MCM/) {
        $carbons = mcm_n_carbon($file);
    } elsif ($run =~ /CRI|RADM|RACM|CB/) {
        $carbons = carbons_others($file);
    } elsif ($run eq "MOZART_tagging") {
        $carbons = mozart_n_carbon($file);
    } else {
        print "$run doesn't match\n";
    }
    return $carbons;
}

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub mozart_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my $words = join ',', (<$in>);
    close $in;
    my ($string) = $words =~ /Solution(.*?)End\sSolution/s;
    $string =~ s/^\s+,//;
    $string =~ s/,\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/RO2(.*?)->//;
    $string =~ s/ROOH(.*?)->//;
    my @species = split ',', $string;
    my %carbons;
    foreach my $species (@species) {
        $species =~ s/^\s+|\s+$//g;
        my $C_number = 0;
        if ($species !~ /->/ and $species !~ /(C[0-9])/) {
            $C_number ++ while ($species =~ m/C/g);
            $carbons{$species} = $C_number;
        } elsif ($species !~ /->/ and $species =~ /(C[0-9])/) { 
            my ($c_nr) = $species =~ /(C[0-9]+)/s;
            $c_nr =~ s/C//; 
            $C_number = $c_nr;
            $carbons{$species} = $C_number;
        } else {
            my ($mech, $molecule) = split ' -> ', $species;
            $mech =~ s/^\s+|\s+$//g;
            if ($molecule =~ /(C[0-9]+)/) { 
                my ($c_nr) = $molecule =~ /(C[0-9]+)/s;
                $c_nr =~ s/C//; 
                $C_number = $c_nr;
                $carbons{$mech} = $C_number;
            } else {
                $C_number ++ while ($molecule =~ m/C/g);
                $carbons{$mech} = $C_number;
            }
        }
    } 
    return \%carbons;
}

sub carbons_others { #get C-number from file names that have species and C# separated by space
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Cannot open file $file: $!";
    my (@lines) = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}
