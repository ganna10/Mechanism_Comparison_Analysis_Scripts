#! /usr/bin/env perl
# Show reactions contributing to pentane C1 OxPE in each mechanism
# Version 0: Jane Coates 17/11/2014
# Version 1: Jane Coates 5/1/2015

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

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my (%families, %weights, %plot_data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $spc_file = "$base/${mechanism}_tagged/gas.spc";
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $carbons_file = "$base/${mechanism}_tagged/carbons.txt";
    my $carbons = get_carbons($mechanism, $carbons_file);
    my $RO2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $mechanism);
        ($plot_data{$mechanism}{$VOC}) = get_data($mecca, $kpp, $mechanism, $mech_species, $carbons, $spc_file);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "CH3O2 + NO" = "#6c254f", "MO2 + NO" = "#6c254f", "MEO2 + NO" = "#6c254d",
                            "HCHO + OH" = "#f9c500", "CH2O + OH" = "#f9c500", "FORM + OH" = "#f9c500",
                            "CO + OH" = "#0e5c28",
                            "HCHO + hv" = "#ef6638", "CH2O + hv" = "#ef6638", "FORM + hv" = "#ef6638",
                            "OH + PAR" = "#6db875",
                            "HOCH2OO" = "#0352cb") `);

$R->run(q` plotting = function (data) { plot = ggplot(data, aes(x = Time, y = Rate, fill = OxPE)) ;
                                        plot = plot + geom_bar(stat = "identity") ;
                                        plot = plot + ggtitle(data$Mechanism) ;
                                        plot = plot + theme_bw() ;
                                        plot = plot + scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) ;
                                        plot = plot + theme(plot.title = element_text(size = 200, face = "bold")) ;
                                        plot = plot + theme(axis.title.x = element_blank()) ;
                                        plot = plot + theme(axis.title.y = element_blank()) ;
                                        plot = plot + theme(axis.text.x = element_text(size = 160, angle = 45, vjust = 0.5)) ;
                                        plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                        plot = plot + theme(legend.title = element_blank()) ;
                                        plot = plot + theme(legend.key = element_blank()) ;
                                        plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                        plot = plot + theme(legend.text = element_text(size = 120)) ;
                                        plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                        plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                        plot = plot + theme(axis.ticks.length = unit(2, "cm")) ;
                                        plot = plot + theme(axis.ticks.margin = unit(1, "cm")) ;
                                        plot = plot + theme(panel.grid.major = element_blank()) ;
                                        plot = plot + theme(panel.grid.minor = element_blank()) ;
                                        plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$OxPE))) ;
                                        return(plot) } `);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plots = c() `);
foreach my $mechanism (sort keys %plot_data) {
    foreach my $VOC (sort keys %{$plot_data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $ref (@{$plot_data{$mechanism}{$VOC}}) {
            foreach my $OxPE (sort keys %$ref) {
                next if ($OxPE eq "CH3O2" or $OxPE eq "MO2 + NO3");
                $R->set('oxpe', $OxPE);
                $R->set('rate', [ map { $_ } $ref->{$OxPE}->dog ]);
                $R->run(q` pre[oxpe] = rate `);
            }
        }
        $R->set('mechanism', $mechanism);
        $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "OxPE", value.name = "Rate") `,
                q` plot = plotting(pre) `,
                q` plots = c(plots, list(plot)) `,
        );
    }
}
#my $p = $R->run(q` print(pre) `);
#print "$p\n";

$R->run(q` CairoPDF(file = "OxPEs_C1_pentane.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                                plots[[2]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                                plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                plots[[6]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                plots[[7]],
                                                plots[[8]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                plots[[9]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                                                nrow = 3),
                                    nrow = 1, ncol = 1,
                                    sub = textGrob("\n", gp = gpar(fontsize = 50)),
                                    left = textGrob("Ox Production Efficiency\n", rot = 90, gp = gpar(fontsize = 180, fontface = "bold"), vjust = 0.8) ) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $VOC, $carbons, $species_file) = @_;
    my %carbons = %$carbons;
    my @ho2_reservoirs = get_ho2_reservoirs($kpp, $species_file) unless ($mechanism =~ /MCM/);
    if (@ho2_reservoirs == 0) {
        $families{"HO2x_$VOC"} = [ qw( HO2 HO2NO2 ) ];
    } else {
        $families{"HO2x_$VOC"} = [ qw( HO2 HO2NO2 ), @ho2_reservoirs ];
    }
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
                    if (defined $carbons{$lookup} and $carbons{$lookup} == 1) {
                        my ($reaction_string) = $kpp->reaction_string($reaction);
                        $reaction_string =~ s/_(.*?)\b//g;
                        my ($label, $rest) = split / = /, $reaction_string;
                        if ($label eq "ROR" and $mechanism =~ /CB/) {
                            $label = "OH + PAR";
                        } elsif ($label eq "CH3O" and $mechanism =~ /MCM/) {
                            $label = "CH3O2 + NO";
                        }
                        $production{$label} += $rate(1:$NTIME-2);
                    } elsif ($lookup =~ /\bCO\b/) {
                        my ($reaction_string) = $kpp->reaction_string($reaction);
                        $reaction_string =~ s/_(.*?)\b//g;
                        my ($label, $rest) = split / = /, $reaction_string;
                        $production{$label} += $rate(1:$NTIME-2);
                    } elsif ( $carbons{$lookup} != 1) {
                        next;
                    } else {
                        print "Nothing found for $lookup\n";
                    }
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
                        if (defined $carbons{$lookup} and $carbons{$lookup} == 1) {
                            my ($reaction_string) = $kpp->reaction_string($reaction);
                            $reaction_string =~ s/_(.*?)\b//g;
                            my ($label, $rest) = split / = /, $reaction_string;
                            if ($label eq "ROR" and $mechanism =~ /CB/) {
                                $label = "OH + PAR";
                            } elsif ($label eq "CH3O" and $mechanism =~ /MCM/) {
                                $label = "CH3O2 + NO";
                            }
                            $production{$label} += $rate(1:$NTIME-2);
                        } elsif ($carbons{$lookup} != 1) {
                            next;
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

    my $others = 1e-1;
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others) {
            $production{"Production Others"} += $production{$reaction};
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
    my @prod_sorted_data = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    
    my @final_sorted_data;
    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $production{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $production{'Production Others'} } if (defined $production{'Production Others'}); 

    return \@final_sorted_data;
}

sub get_ho2_reservoirs {
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @species;
    for (<$in>) {
        push @species, split /\s+/, $_; 
    }
    close $in;
    my %ho2_reservoirs;
    foreach my $species (@species) {
        my ($reactions) = $kpp->reacting_with($species, 'HO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            next unless ($products->[0] =~ /HOCH2OO|HCO3/);
            $ho2_reservoirs{$products->[0]} += 1;
        }   
    }   
    return keys %ho2_reservoirs;
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
    } elsif ($run =~ /MOZART/) {
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
