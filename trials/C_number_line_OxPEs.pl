#! /usr/bin/env perl
# Calculate OxPE of degradation products C number in each mechanism by normalising over total Ox consumption
# Version 0: Jane Coates 14/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT ;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#my @runs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "MCM v3.2", "MCM v3.1", "CRI v2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my $index = 0;
my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $carbons_file = "$base/$run/carbons.txt";
    my $carbons = get_carbons($run, $carbons_file);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Toluene );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $run);
        ($plot_data{$mechanisms[$index]}{$VOC}) = get_data($mecca, $kpp, $mechanisms[$index], $mech_species, $carbons);
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %plot_data) {
    foreach my $VOC (sort keys %{$plot_data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $carbon (sort keys %{$plot_data{$mechanism}{$VOC}}) {
            $R->set('carbon', $carbon);
            $R->set('oxpe', [map { $_ } $plot_data{$mechanism}{$VOC}{$carbon}->dog]);
            $R->run(q` pre[carbon] = oxpe `);
        }
        $R->set('mechanism', $mechanism);
        $R->set('voc', $VOC);
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
                q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre$VOC = rep(voc, length(Time)) `,
                q` pre = melt(pre, id.vars = c("Time", "Mechanism", "VOC"), variable.name = "C.number", value.name = "OxPE") `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` data$C.number = factor(data$C.number, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `);
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRI v2" = "#ef6638", "MCM v3.1" = "#000000", "MCM v3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line(size = 14) `,
        q` plot = plot + geom_point(size = 35) `,
        q` plot = plot + facet_wrap( ~ C.number) `,
        q` plot = plot + ylab("Ox Production Efficiency\n") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text.x = element_text(size = 200, face = "bold")) `,
        #q` plot = plot + theme(strip.text.y = element_text(size = 200, face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        #q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 140)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
        q` plot = plot + theme(axis.ticks.length = unit(2, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
);

$R->run(q` CairoPDF(file = "OxPEs_line_by_C_number.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $VOC, $carbons) = @_;
    my %carbons = %$carbons;
    $families{"HO2x_$VOC"} = [ qw( HO2 HO2NO2 ) ];
    $families{"Ox_${mechanism}_$VOC"} = $families{"Ox_$mechanism"};
    my @loop = ("Ox_${mechanism}_$VOC", "HO2x_$VOC");
    my (%production, $consumption);

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

        foreach (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($r_number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $reactants = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $mechanism =~ /RADM2|RACM|CB/);
                    my ($lookup, $tag) = split /_/, $_;
                    if (defined $carbons{$lookup}) {
                        $production{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } elsif ($lookup =~ /\bCO\b/) {
                        $production{"C1"} += $rate(1:$NTIME-2);
                    } else {
                        print "Nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $production{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                }
            }
        }

        if ($species =~ /Ox/ and $mechanism =~ /RADM2|RACM|CB/) {#operator allocation for those mechanisms that use it: RADM2, RACM, RACM2, CBM4, CB05 -> XO2
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
                            $production{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                        } else {
                            print "$mechanism => nothing found for $lookup\n";
                        }
                    }
                } 
            } 
        }

        foreach (0..$#$consumers) {
            next if ($species eq "HO2x");
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            $consumption += $rate(1:$NTIME-2);
        }
    }

    $production{$_} /= -$consumption foreach (keys %production); #normalise each processes' Ox production by total consumption

    foreach my $carbon (keys %production) {
        my $reshape = $production{$carbon}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$carbon} = $integrate;
    }
    return \%production;
}

sub get_no2_reservoirs {
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