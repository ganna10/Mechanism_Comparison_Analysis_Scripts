#! /usr/bin/env perl
# OxPE contributing reaction during NC5H12 degradation in each mechanism by normalising over total Ox consumption
# Version 0: Jane Coates 31/12/2014

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
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT ;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#my @mechanisms = qw( CB05 );
my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $mechanism);
        ($data{$mechanism}{$VOC}) = get_data($mecca, $kpp, $mechanism, $mech_species);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $ref (@{$data{$mechanism}{$VOC}}) {
            foreach my $reaction (sort keys %$ref) {
                $R->set('reaction', $reaction);
                $R->set('oxpe', [map { $_ } $ref->{$reaction}->dog]);
                $R->run(q` pre[reaction] = oxpe `);
            }
        }
        $R->set('mechanism', $mechanism);
        $R->set('voc', $VOC);
        $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre$VOC = rep(voc, length(Time)) `,
                q` pre = gather(pre, Reaction, OxPE, -Time, -Mechanism, -VOC) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + ylab("Ox Production Efficiency") `,
        q` plot = plot + theme_bw() `,
        #q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) `,
        q` plot = plot + theme(legend.margin = unit(0, "lines")) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
);

$R->run(q` CairoPDF(file = "pentane_OxPE_by_reaction.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $VOC) = @_;
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
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $products) = split / = /, $reaction_string;
            if ($species =~ /Ox/ and $reactants =~ /XO2/ and $mechanism =~ /RADM2|RACM|CB/) {
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
                    $reaction_string =~ s/_(.*?)\b//g;
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $production{$reactants} += $rate(1:$NTIME-2);
                } 
            } else {
                $production{$reactants} += $rate(1:$NTIME-2);
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
    my $others = 2;
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others) {
            $production{"Others"} += $production{$reaction};
            delete $production{$reaction};
        }
    }

    foreach my $reaction (keys %production) {
        if ($VOC =~ /TOL/) {
            if ($mechanism =~ /MOZ/) {
                $production{$reaction} *= 0.232;
            } elsif ($mechanism =~ /RADM2|RACM\b/) {
                $production{$reaction} *= 0.667;
            } elsif ($mechanism =~ /RACM2/) {
                $production{$reaction} *= 0.868;
            }
        } else { #pentane
            if ($mechanism =~ /CB/) {
                $production{$reaction} /= 5;
            } elsif ($mechanism =~ /MOZ/) {
                $production{$reaction} *= 0.146;
            } elsif ($mechanism =~ /RA/) {
                $production{$reaction} *= 0.264;
            }
        }
        my $reshape = $production{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$reaction} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Others') ;
        push @final_sorted_data, { $_ => $production{$_} };
    } 
    push @final_sorted_data, { 'Others' => $production{'Others'} } if (defined $production{'Others'}); 
    return \@final_sorted_data;
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
    } elsif ($run =~ /MOZ/) {
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
