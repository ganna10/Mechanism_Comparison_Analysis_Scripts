#! /usr/bin/env perl
# Plot Ox production of carbon numbers as percentage of total Ox production for pentane and toluene for all mechanisms.
# Faceting by day, x-axis is percentage and y-axis is mechanism
# Version 0: Jane Coates 21/9/2014
# Version 1: Jane Coates 1/10/2014 Rounding up carbon numbers to nearest integer and updating colour scheme
# Version 2: Jane Coates 11/11/2014 adapting to local computation
# Version 3: Jane Coates 9/12/2014 re-factoring code for constant emissions

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

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2",  "CBM-IV", "CB05" );
my (%n_carbon, %families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $ro2file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanism"} = get_carbons($mechanism, $carbon_file);
    my @parents = qw( Pentane Toluene );
    foreach my $NMVOC (@parents) {
        my $parent = get_mechanism_species($NMVOC, $mechanism);
        ($data{$mechanism}{$NMVOC}) = get_data($kpp, $mecca, $mechanism, $n_carbon{"Ox_$mechanism"}, $parent);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(scales) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);

foreach my $run (sort keys %data) {
    foreach my $VOC (sort keys %{$data{$run}}) {
        $R->run(q` pre.sort = data.frame(Time) `);
        foreach my $ref (@{$data{$run}{$VOC}}) {
            foreach my $carbons (sort keys %$ref) {
                $R->set('carbons', $carbons);
                $R->set('rate', [map { $_ } $ref->{$carbons}->dog]);
                $R->run(q` pre.sort[carbons] = rate `);
            }
        }
        $R->set('mechanism', $run);
        $R->set('Parent', $VOC); 
        $R->run(q` if("C2.4" %in% colnames(pre.sort)) { pre.sort$C2 = pre.sort$C2 + pre.sort$C2.4 ; pre.sort$C2.4 = NULL }`,
                q` if("C2.9" %in% colnames(pre.sort)) { pre.sort$C3 = pre.sort$C3 + pre.sort$C2.9 ; pre.sort$C2.9 = NULL }`,
                q` if("C3.5" %in% colnames(pre.sort)) { pre.sort$C4 = pre.sort$C4 + pre.sort$C3.5 ; pre.sort$C3.5 = NULL }`,
                q` if("C3.6" %in% colnames(pre.sort)) { pre.sort$C4 = pre.sort$C4 + pre.sort$C3.6 ; pre.sort$C3.6 = NULL }`,
                q` if("C3.9" %in% colnames(pre.sort)) { pre.sort$C4 = pre.sort$C4 + pre.sort$C3.9 ; pre.sort$C3.9 = NULL }`,
                q` if("C4.2" %in% colnames(pre.sort)) { pre.sort$C4 = pre.sort$C4 + pre.sort$C4.2 ; pre.sort$C4.2 = NULL }`,
                q` if("C4.5" %in% colnames(pre.sort)) { pre.sort$C5 = pre.sort$C5 + pre.sort$C4.5 ; pre.sort$C4.5 = NULL }`,
                q` if("C4.8" %in% colnames(pre.sort)) { pre.sort$C5 = pre.sort$C4.8 ; pre.sort$C4.8 = NULL }`,
                q` if("C5.6" %in% colnames(pre.sort)) { pre.sort$C6 = pre.sort$C5.6 ; pre.sort$C5.6 = NULL }`,
                q` if("C6.6" %in% colnames(pre.sort)) { pre.sort$C7 = pre.sort$C6.6 ; pre.sort$C6.6 = NULL }`,
                q` if("C7.1" %in% colnames(pre.sort)) { pre.sort$C7 = pre.sort$C7 + pre.sort$C7.1 ; pre.sort$C7.1 = NULL }`,
                q` if("C7.75" %in% colnames(pre.sort)) { pre.sort$C8 = pre.sort$C7.75 ; pre.sort$C7.75 = NULL }`,
                q` pre.sort$Mechanism = rep(mechanism, length(time)) `,
                q` pre.sort$VOC = rep(Parent, length(time)) `,
                q` pre.sort = gather(pre.sort, C.number, Rate, -Time, -Mechanism, -VOC) `,
                q` pre.sort = pre.sort %>% group_by(Time) %>% do(data.frame(C.number = .$C.number, Mechanism = .$Mechanism, VOC = .$VOC, Percent = .$Rate / sum(.$Rate))) `,
                q` data = rbind(data, pre.sort) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
$R->run(q` my.names = c("C8" = "C8 ", "C7" = "C7 ", "C6" = "C6 ", "C5" = "C5 ", "C4" = "C4 ", "C3" = "C3 ", "C2" = "C2 ", "C1" = "C1 ") `);
$R->run(q` data$C.number = factor(data$C.number, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `);

$R->run(q` plot = ggplot(data, aes(y = Percent, x = Mechanism, fill = C.number)) `,
        q` plot = plot + geom_bar(stat = "identity", width = 0.5) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + facet_grid( Time ~ VOC ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05"))) `,
        q` plot = plot + scale_y_continuous(labels = percent_format(), expand = c(0, 0.009)) `,
        q` plot = plot + ylab("Carbon Number of Degradation Products Percent Contribution to Daily Ox Production") `,
        q` plot = plot + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) `,
        q` plot = plot + theme(legend.margin = unit(0, "lines")) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(panel.margin.x = unit(0.6, "cm")) `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text.y = element_text(face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "carbon_percent_total_Ox_production.pdf", width = 8.7, height = 10) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, $producer_yields, %production_rates);
    my @families = ("Ox_${mechanism}_$VOC", "HO2x_${mechanism}_$VOC");
    foreach my $family (@families) {
        if (exists $families{$family}) { 
            $kpp->family({ 
                    name    => $family,
                    members => $families{$family},
                    weights => $weights{$family},
            });
            $producers = $kpp->producing($family);
            $producer_yields = $kpp->effect_on($family, $producers);  
        } else {
            print "No family found for $family\n";
        } 
        die "No producers found for $family\n" if (@$producers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $family =~ /RADM2|RACM|RACM2|CBM-IV|CB05/);
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) { 
                        $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } elsif ($lookup =~ /CO/) {
                        $production_rates{'C1'} += $rate(1:$NTIME-2);
                    } else {
                        print "nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                }
            }
        }
    }
    
    if ($mechanism =~ /RADM2|RACM|RACM2|CBM-IV|CB05/) {#operator allocation for those mechanisms that use it: RADM2, RACM, RACM2, CBM4, CB05 -> XO2
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
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) {
                        $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } else {
                        print "$mechanism => nothing found for $lookup\n";
                    }
                }
            } 
        } 
    }

    #get parent species emissions for each mechanism
#    my $dt = $mecca->dt->at(0); #model time step
#    my $parent_emissions;
#    if ($mechanism =~ /CBM-IV/ or $mechanism =~ /CB05/) {
#        if ($VOC =~ /NC5H12/) {
#            my $name = "PAR_NC5H12";
#            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
#            $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
#        } elsif ($VOC =~ /TOLUENE/) {
#            my $name = "TOL_TOLUENE";
#            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
#            $parent_emissions += $parent_source->sum * $dt ; 
#        } else {
#            print "No emissions data for $VOC\n";
#        }
#    } else {
#        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
#        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
#    }
#    
#    #normalise by dividing reaction rate of intermediate (molecules (Product) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
#    $production_rates{$_} /= $parent_emissions foreach (sort keys %production_rates);
    # no more normalisation by emissions as constant emissions are used
    foreach my $C (keys %production_rates) {
        my $reshape = $production_rates{$C}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_rates{$C} = $integrate;
    }
    my @prod_sorted_data = sort { $a cmp $b } keys %production_rates; 
    my @final_sorted_data;
    push @final_sorted_data, { $_ => $production_rates{$_} } foreach (@prod_sorted_data); 
    return \@final_sorted_data;
}

sub get_carbons {
    my ($run, $file) = @_;
    my $carbons;
    if ($run =~ /MCM/) {
        $carbons = mcm_n_carbon($file);
    } elsif ($run =~ /MOZART/) {
        $carbons = mozart_n_carbon($file);
    } elsif ($run =~ /CRI|RADM2|RACM|CB/) {
        $carbons = carbons_others($file);
    } else {
        print "$run doesn't match\n";
    }
    return $carbons;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
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

sub get_mechanism_species {
    my ($NMVOC, $run) = @_;

    my $mechanism_species;
    if ($NMVOC eq "Pentane") {
        if ($run =~ /MCM|CRI|CB/) {
            $mechanism_species = "NC5H12";
        } elsif ($run =~ /MOZART/) {
            $mechanism_species = "BIGALK";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "HC5";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } elsif ($NMVOC eq "Toluene") {
        if ($run =~ /MCM|CRI|MOZART|CB/) {
            $mechanism_species = "TOLUENE";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "TOL";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } else {
        print "No $NMVOC data\n";
    }
    return $mechanism_species;
}
