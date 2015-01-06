#! /usr/bin/env perl
# Plot Ox production rates of the degradation products for octane
# Version 0: Jane Coates 5/1/2015

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
#my @mechanisms = qw( RADM2 CB05 );
my (%families, %weights, %n_carbon, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanism"} = get_carbons($mechanism, $carbon_file);
    my $parent = get_mechanism_species("Octane", $mechanism);
    ($data{$mechanism}) = get_data($kpp, $mecca, $mechanism, $n_carbon{"Ox_$mechanism"}, $parent);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(scales) `,
);

$R->set('time', [ ("Day 1", "Day 2") ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (keys %data) {
    $R->run(q` pre = data.frame(Time = time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $C (sort keys %$ref) {
            $R->set('c', $C);
            $R->set('rate', [ map { $_ } $ref->{$C}->dog ]);
            $R->run(q` pre[c] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` if("C2.4" %in% colnames(pre)) { pre$C2 = pre$C2 + pre$C2.4 ; pre$C2.4 = NULL }`,
            q` if("C2.9" %in% colnames(pre)) { pre$C3 = pre$C3 + pre$C2.9 ; pre$C2.9 = NULL }`,
            q` if("C3.5" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.5 ; pre$C3.5 = NULL }`,
            q` if("C3.6" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.6 ; pre$C3.6 = NULL }`,
            q` if("C3.9" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.9 ; pre$C3.9 = NULL }`,
            q` if("C4.2" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C4.2 ; pre$C4.2 = NULL }`,
            q` if("C4.5" %in% colnames(pre)) { pre$C5 = pre$C5 + pre$C4.5 ; pre$C4.5 = NULL }`,
            q` if("C4.8" %in% colnames(pre)) { pre$C5 = pre$C4.8 ; pre$C4.8 = NULL }`,
            q` if("C5.6" %in% colnames(pre)) { pre$C6 = pre$C5.6 ; pre$C5.6 = NULL }`,
            q` if("C6.6" %in% colnames(pre)) { pre$C7 = pre$C6.6 ; pre$C6.6 = NULL }`,
            q` if("C7.1" %in% colnames(pre)) { pre$C7 = pre$C7 + pre$C7.1 ; pre$C7.1 = NULL }`,
            q` if("C7.75" %in% colnames(pre)) { pre$C8 = pre$C7.75 ; pre$C7.75 = NULL }`,
            q` if("C7.9" %in% colnames(pre)) { pre$C8 = pre$C7.9 ; pre$C7.9 = NULL }`,
            q` pre$Mechanism = rep(mechanism, length(time)) `,
            q` pre = gather(pre, C, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `);
$R->run(q` data$C = factor(data$C, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `);
#my $p = $R->run(q` print(as.data.frame(data)) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = C)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "dodge") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        #q` plot = plot + scale_y_continuous(labels = percent_format()) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "octane_day2_Ox_production_rates_C_number.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, $producer_yields, %production_rates, $consumers, $consumer_yields, %consumption_rates);
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
            $consumers = $kpp->consuming($family);
            $consumer_yields = $kpp->effect_on($family, $consumers);  
        } else {
            print "No family found for $family\n";
        } 
        die "No producers found for $family\n" if (@$producers == 0);
        die "No consumers found for $family\n" if (@$consumers == 0);
        
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

        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $family =~ /RADM2|RACM|RACM2|CBM-IV|CB05/);
                    my ($lookup, $rest) = split '_', $_;
                    if (defined $carbons{$lookup}) { 
                        $consumption_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } elsif ($lookup =~ /CO/) {
                        $consumption_rates{'C1'} += $rate(1:$NTIME-2);
                    } else {
                        print "nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $consumption_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
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
    remove_common_processes(\%production_rates, \%consumption_rates);

    foreach my $C (keys %production_rates) {
        if ($mechanism =~ /MOZ/) {
            $production_rates{$C} *= 0.01;
        } elsif ($mechanism =~ /CB/) {
            $production_rates{$C} /= 8;
        }
        my $reshape = $production_rates{$C}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $production_rates{$C} = $integrate(0:3:2);
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
    if ($NMVOC eq "Octane") {
        if ($run =~ /MCM|CRI|CB/) {
            $mechanism_species = "NC8H18";
        } elsif ($run =~ /MOZART/) {
            $mechanism_species = "BIGALK";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "HC8";
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
