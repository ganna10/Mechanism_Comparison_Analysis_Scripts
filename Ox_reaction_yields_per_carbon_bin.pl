#! /usr/bin/env perl
# Plot Ox production efficiency of carbon numbers as percentage of total Ox production for pentane and toluene for all mechanisms.
# Faceting by day, x-axis is percentage and y-axis is mechanism
# Version 0: Jane Coates 2/10/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2",  "CBM-IV", "CB05" );
#my @runs = qw( RACM_tagging CB05_tagging ) ;
#my @mechanisms = qw( RACM CB05 );
my $array_index = 0;

my (%n_carbon, %families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $ro2file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanisms[$array_index]"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$array_index]"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanisms[$array_index]"} = get_carbons($run);
    my @parents = qw( Pentane Toluene );
    foreach my $NMVOC (@parents) {
        my $parent = get_mechanism_species($NMVOC, $run);
        ($plot_data{$mechanisms[$array_index]}{$NMVOC}) = get_data($kpp, $mecca, $mechanisms[$array_index], $n_carbon{"Ox_$mechanisms[$array_index]"}, $parent);
    }
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(scales) `,
        q` library(grid) `,
);

$R->run(q` data = data.frame() `);
foreach my $run (sort keys %plot_data) {
    foreach my $VOC (sort keys %{$plot_data{$run}}) {
        $R->set('mechanism', $run);
        $R->set('Parent', $VOC);
        $R->run(q` pre = data.frame(carbon.number = numeric(0), Ox.yield = numeric(0)) `);
        foreach my $ref (@{$plot_data{$run}{$VOC}}) {
            foreach my $carbons (sort keys %$ref) {
                $R->set('carbons', $carbons);
                $R->set('yield', $ref->{$carbons});
                $R->run(q` pre = rbind(data.frame(carbon.number = carbons, Ox.yield = yield), pre) `);
            }
        }

        $R->run(q` pre$Mechanism = rep(mechanism, length(pre$Ox.yield)) `,
                q` pre$VOC = rep(Parent, length(pre$Ox.yield)) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
$R->run(q` my.names = c("C8" = "C8 ", "C7" = "C7 ", "C6" = "C6 ", "C5" = "C5 ", "C4" = "C4 ", "C3" = "C3 ", "C2" = "C2 ", "C1" = "C1 ") `);
$R->run(q` data$carbon.number = factor(data$carbon.number, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `);
$R->run(q` data$VOC = factor(data$VOC, labels = c("Pentane\n", "Toluene\n")) `);

$R->run(q` plot = ggplot(data, aes(y = Ox.yield, x = Mechanism)) `,
        q` plot = plot + geom_bar(stat = "identity", width = 0.5) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + facet_grid( carbon.number ~ VOC, scales = "free_x" ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05"))) `,
        #q` plot = plot + scale_y_continuous(labels = percent_format()) `,
        #q` plot = plot + ylab("\nCarbon Number of Degradation Products Percent Contribution to Daily Ox Production\n") `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text.y = element_text(size = 160, face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.text.x = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 160, face = "bold")) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 140)) `,
        #q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "Ox_production_efficiency_by_C_number.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, $producer_yields, %yields, %uniques, %bin_size);
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

        #check that species reactions are found
        die "No producers found for $family\n" if (@$producers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $yield =  $producer_yields->[$_];
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $family =~ /RADM2|RACM|RACM2|CBM-IV|CB05/);
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) { 
                        $yields{"C$carbons{$lookup}"} += $yield;
                        $uniques{$lookup} += 1;
                    } elsif ($lookup =~ /CO/) {
                        $yields{'C1'} += $yield;
                        $uniques{"CO"} += 1;
                    } else {
                        print "nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $yields{"C$carbons{$lookup}"} += $yield;
                    $uniques{"HC5"} += 1;
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
            my ($reactants) = $kpp->reactants($reaction);
            my $yield = $op_producer_yields->[$_];
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) {
                        $yields{"C$carbons{$lookup}"} += $yield;
                        $uniques{$lookup} += 1;
                    } else {
                        print "$mechanism => nothing found for $lookup\n";
                    }
                }
            } 
        } 
    }

    print "$mechanism and $VOC\n";
    foreach my $species (sort keys %uniques) {
        $bin_size{"C$carbons{$species}"} += 1;
    }
    foreach my $cs (sort keys %bin_size) {
        print "\t$cs => $bin_size{$cs}\n";
    }

    #normalise by number of Ox producing reactions per size bin
    #$yields{$_} /= $bin_size{$_} foreach (sort keys %yields);

    my @all_Cs = qw( C1 C2 C2.4 C2.9 C3 C3.5 C3.6 C3.9 C4 C4.2 C4.5 C4.8 C5 C5.6 C6 C6.6 C7 C7.1 C7.75 C8 );
    my @present_Cs = keys %yields;
    foreach my $C (@all_Cs) {
        next if ($C ~~ @present_Cs);
        $yields{$C} = 0;
    }

    $yields{'C2'} += $yields{'C2.4'};
    delete $yields{'C2.4'};
    $yields{'C3'} += $yields{'C2.9'};
    delete $yields{'C2.9'};
    $yields{'C4'} += $yields{'C3.5'} + $yields{'C3.6'} + $yields{'C3.9'} + $yields{'C4.2'};
    delete $yields{'C3.5'};
    delete $yields{'C3.6'};
    delete $yields{'C3.9'};
    delete $yields{'C4.2'};
    $yields{'C5'} += $yields{'C4.5'} + $yields{'C4.8'};
    delete $yields{'C4.5'};
    delete $yields{'C4.8'};
    $yields{'C6'} += $yields{'C5.6'};
    delete $yields{'C5.6'};
    $yields{'C7'} += $yields{'C6.6'} + $yields{'C7.1'};
    delete $yields{'C6.6'};
    delete $yields{'C7.1'};
    $yields{'C8'} += $yields{'C7.75'};
    delete $yields{'C7.75'};
    
    my @prod_sorted_data = sort { $a cmp $b } keys %yields; 
    my @final_sorted_data;
    push @final_sorted_data, { $_ => $yields{$_} } foreach (@prod_sorted_data);

    return \@final_sorted_data;
}

sub get_carbons {
    my ($run) = @_;
    my $mechanism_base = "/work/users/jco/Mechanisms";
    my ($file, $carbons);
    if ($run eq "MCM_3.1_tagged_3.2rates") {
        $file = "$mechanism_base/MCM/MCM-3.1/species_smiles.txt";
        $carbons = mcm_n_carbon($file);
    } elsif ($run eq "MCM_3.2_tagged") {
        $file = "$mechanism_base/MCM/MCM-3.2/smiles.out";
        $carbons = mcm_n_carbon($file);
    } elsif ($run eq "CRI_tagging") {
        $file = "$mechanism_base/CRI/CRI_v2_full/carbons.txt";
        $carbons = carbons_others($file);
    } elsif ($run eq "MOZART_tagging") {
        $file = "$mechanism_base/MOZART/MOZART/chem_mech.in";
        $carbons = mozart_n_carbon($file);
    } elsif ($run eq "RADM2_tagged") {
        $file = "$mechanism_base/RADM2/carbon_numbers.txt";
        $carbons = carbons_others($file);
    } elsif ($run eq "RACM_tagging") {
        $file = "$mechanism_base/RACM/carbon_numbers.txt";
        $carbons = carbons_others($file);
    } elsif ($run eq "RACM2_tagged") {
        $file = "$mechanism_base/RACM2/carbon_numbers.txt";
        $carbons = carbons_others($file);
    } elsif ($run eq "CBM4_tagging") {
        $file = "$mechanism_base/CBM-IV/carbons.txt";
        $carbons = carbons_others($file);
    } elsif ($run eq "CB05_tagging") {
        $file = "$mechanism_base/CB05/carbons.txt";
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
