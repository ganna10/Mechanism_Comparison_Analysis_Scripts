#! /usr/bin/env perl
# pentane Ox net Ox production budget for each mechanism
# Version 0: Jane Coates 23/10/2014

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
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
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
    my $carbon_file = "$base/$run/carbons.txt";
    my $carbons = get_carbons($run, $carbon_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 HO2NO2 NO3 N2O5 O1D O ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    my @VOCs = qw( Pentane );
    foreach my $VOC (@VOCs) {
        my $mech_species = get_model_name($VOC, $run);
        $plot_data{$mechanisms[$index]} = get_data($kpp, $mecca, $mechanisms[$index], $mech_species, $carbons);
    }
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(grid) `,
        q` library(gridExtra) `,
        q` library(dplyr) `,
        q` library(scales) `,
);

$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
$R->run(q` my.names = c("C8" = "C8 ", "C7" = "C7 ", "C6" = "C6 ", "C5" = "C5 ", "C4" = "C4 ", "C3" = "C3 ", "C2" = "C2 ", "C1" = "C1 ") `);
$R->run(q` scientific_10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `);

$R->run(q` plotting = function (data, mechanism) {  plot = ggplot(data, aes(x = Time, y = net.rate, fill = C.number)) ;
                                                    plot = plot + geom_bar(data = subset(data, net.rate > 0), stat = "identity") ; 
                                                    plot = plot + geom_bar(data = subset(data, net.rate < 0), stat = "identity") ; 
                                                    plot = plot + scale_y_continuous(limits = c(0, 3.1), breaks = seq(0, 3.1, 1), labels = scientific_10) ;
                                                    plot = plot + theme_bw() ;
                                                    plot = plot + ggtitle(mechanism) ;
                                                    plot = plot + theme(plot.title = element_text(size = 140, face = "bold")) ;
                                                    plot = plot + theme(axis.title.x = element_blank()) ;
                                                    plot = plot + theme(axis.title.y = element_blank()) ;
                                                    plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                                    plot = plot + theme(axis.text.x = element_text(size = 140)) ;
                                                    plot = plot + theme(panel.grid.major = element_blank()) ;
                                                    plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                    plot = plot + theme(legend.position = "bottom") ;
                                                    plot = plot + theme(legend.key = element_blank()) ;
                                                    plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                                    plot = plot + theme(legend.title = element_blank()) ;
                                                    plot = plot + theme(legend.text = element_text(size = 140)) ;
                                                    plot = plot + scale_fill_manual(values = my.colours, labels = my.names) ;
                                                    return(plot) } `,
);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('Time', [@days]);
$R->run(q` Plots = list() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $carbon (sort keys %{$plot_data{$run}} ) {
        $R->set('C.number', $carbon);
        $R->set('net.rate', [map { $_ } $plot_data{$run}{$carbon}->dog]);
        $R->run(q` pre[C.number] = net.rate `);
    }
    $R->set('Mechanism', $run);
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
            q` pre = melt(pre, id.vars = c("Time"), variable.name = "C.number", value.name = "net.rate") `,
            q` pre$C.number = factor(pre$C.number, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `,
            q` pre$Mechanism = rep(Mechanism, length(pre$Time)) `,
            q` plot = plotting(pre, Mechanism) `,
            q` Plots = c(Plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(pre) `);
#print "$p\n"; 

$R->run(q` CairoPDF(file = "pentane_net_Ox_daytime_budget.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(arrangeGrob(Plots[[1]] ,
                                                Plots[[2]] ,
                                                Plots[[3]] ,
                                                Plots[[4]] ,
                                                Plots[[5]] ,
                                                Plots[[6]] ,
                                                Plots[[7]] ,
                                                Plots[[8]] ,
                                                Plots[[9]] ,
                                                nrow = 3),
                                    nrow = 1, ncol = 1,
                                    left = textGrob(expression(bold(paste("\nNet Ox Production Budget (molecules (Ox) ", cm^-3, s^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $VOC, $carbons) = @_;
    my %carbons = %$carbons;
    $families{"HO2x_${mechanism}"} = [ qw( HO2 HO2NO2 )];
    my @loop = ("Ox_${mechanism}", "HO2x_${mechanism}");

    my ($producers, $producer_yields, %production_rates, $consumers, $consumer_yields, %consumption_rates);
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
        
        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number);
            next if ($rate->sum == 0);
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $species =~ /RADM2|RACM|CB/);
                    my ($lookup, $tag) = split /_/, $_;
                    if (defined $carbons{$lookup}) {
                        $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } elsif ($lookup =~ /\bCO\b/) {
                        $production_rates{"C1"} += $rate(1:$NTIME-2);
                    } else {
                        print "Nothing found for $lookup\n";
                    }
                } elsif ($_ =~ /HC5\b/) {
                    my $lookup = "HC5";
                    $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
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
                            $production_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                        } else {
                            print "$mechanism => nothing found for $lookup\n";
                        }
                    }
                } 
            } 
        }
        
        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number);
            next if ($rate->sum == 0);
            my ($reactants) = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $species =~ /RADM2|RACM|CB/);
                    my ($lookup, $tag) = split /_/, $_;
                    if (defined $carbons{$lookup}) {
                        $consumption_rates{"C$carbons{$lookup}"} += $rate(1:$NTIME-2);
                    } else {
                        print "Nothing found for $lookup\n";
                    }
                }
            }
        } 
    }

    my $dt = $mecca->dt->at(0);
    my $n_per_day = 43200 / $dt;
    my $n_days = int ($NTIME / $n_per_day);
    my %net_rates;
    foreach my $carbon (sort keys %production_rates) {
        my $reshape = $production_rates{$carbon}->copy->reshape($n_per_day, $n_days);
        my $integrated = $reshape->sumover;
        $integrated = $integrated(0:13:2); #choose only day period
        #print "production: $carbon => $integrated\n";
        $production_rates{$carbon} = $integrated;
    }
    foreach my $carbon (sort keys %consumption_rates) {
        my $reshape = $consumption_rates{$carbon}->copy->reshape($n_per_day, $n_days);
        my $integrated = $reshape->sumover;
        $integrated = $integrated(0:13:2); #choose only day period
        #print "consumption: $carbon => $integrated\n";
        $consumption_rates{$carbon} = $integrated;
        $net_rates{$carbon} = $production_rates{$carbon} + $consumption_rates{$carbon};
    }

    my $parent_emissions;
    if ($mechanism =~ /CB/) {
        if ($VOC =~ /NC5H12/) {
            my $name = "PAR_NC5H12";
            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
            $parent_emissions += $parent_source->sum * $dt / 5; #NC5H12 => 5 PAR
        } elsif ($VOC =~ /TOLUENE/) {
            my $name = "TOL_TOLUENE";
            my $parent_source = $mecca->balance($name); #in molecules (VOC)/cm3/s
            $parent_emissions += $parent_source->sum * $dt ; 
        } else {
            print "No emissions data for $VOC\n";
        }
    } else {
        my $parent_source = $mecca->balance($VOC); #in molecules (VOC)/cm3/s
        $parent_emissions += $parent_source->sum * $dt; #in molecules (VOC)/cm3
    }
    
    $net_rates{$_} = $net_rates{$_} * $dt / $parent_emissions foreach (sort keys %net_rates);#normalise by dividing reaction rate of intermediate (molecules (Product) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    #print "Net: $_ => $net_rates{$_}\n" foreach sort keys %net_rates;
    return \%net_rates;
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
