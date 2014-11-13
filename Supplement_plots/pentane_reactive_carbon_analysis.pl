#! /usr/bin/env perl
# analysis of rate of reactive carbon loss during pentane degradation in each mechanism
# Version 0: Jane Coates 23/9/2014
# Version 1: Jane Coates 13/11/2014 refining plot for inclusion in supplement

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog; 
my @time_blocks;
foreach my $time (@time_axis) {#map to day and night
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

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2",  "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
my $array_index = 0;
my $NMVOC = "Pentane";

my (%n_carbon, %families, %weights, %plot_data, %legend);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $carbon_file = "$base/$run/carbons.txt";
    my $ro2file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanisms[$array_index]"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$array_index]"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanisms[$array_index]"} = get_carbons($run, $carbon_file);
    my $parent = get_mechanism_species($NMVOC, $run);
    ($plot_data{$mechanisms[$array_index]}, $legend{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index], $n_carbon{"Ox_$mechanisms[$array_index]"}, $parent);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "Consumption Others" = "#cc6329",
                            "C2O3 + NO" = "#6c254f", "CH3CO3 + NO" = "#6c254f", "ACO3 + NO" = "#6c254f",
                            "ROR" = "#1b695b", 
                            "C2H5CO3 + NO" = "#0c3f78",
                            "CO2C4CO3 + NO" = "#8ed6d2",
                            "HOC2H4CO3 + NO" = "#f9c500", 
                            "HOCH2CO3 + NO" = "#4c9383",
                            "CO2C3CO3 + NO" = "#86b650",
                            "CH3CO3" = "#623812",
                            "ALKO2 + NO" = "#77aecc", 
                            "NO + RCO3" = "#c9a415",
                            "ETHP + NO" = "#ae4901",
                            "HC3P + NO" = "#b569b3", 
                            "HC5 + OH" = "#0352cb",
                            "HC5P + NO" = "#0e5c28",
                            "MEKP + NO" = "#58691b",
                            "KETP + NO" = "#ef6638", 
                            "OH + ONIT" = "#58691b") `,
);

$R->run(q` plotting = function (data, mechanism, legend) {  plot = ggplot(data, aes(x = Time, y = Carbon.loss.rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, data$Carbon.loss.rate < 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, data$Carbon.loss.rate > 0), stat = "identity") ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank());
                                                            plot = plot + theme(axis.ticks.margin = unit(1, "cm")) ;
                                                            plot = plot + theme(axis.ticks.length = unit(2, "cm")) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.01)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.01)) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 150, angle = 45, vjust = 0.5)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                                            plot = plot + theme(plot.title = element_text(size = 200, face = "bold")) ;
                                                            plot = plot + theme(legend.text = element_text(size = 140)) ;
                                                            plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                                            plot = plot + scale_y_continuous(limits = c(-2.4e8, 2e7), breaks = seq(-2.4e8, 2e7, 4e7)) ;
                                                            plot = plot + scale_fill_manual(values = my.colours, limits = legend) ;
                                                            return(plot) } `);

$R->set('Time', [@time_blocks]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('mechanism', $run);
    $R->set('legend', [@{$legend{$run}}]);
    $R->run(q` data = ddply(data, .(Time), colwise(sum)) `,
            q` data = data[1:7,] `,
            q` data = melt(data, id.vars = c("Time"), variable.name = "Reaction", value.name = "Carbon.loss.rate") `,
            q` reaction.levels = rev(levels(factor(data$Reaction))) `,
            q` data$Reaction = ordered(data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(data, mechanism, legend) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` CairoPDF(file = "pentane_reactive_carbon_loss.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(   arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[2]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[6]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[7]], 
                                                    plots[[8]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[9]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 200), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, %carbon_loss_rate, %carbon_gain_rate);
    my @families = ("Ox_${mechanism}_$VOC", "HO2x_${mechanism}_$VOC");
    foreach my $family (@families) {
        if (exists $families{$family}) { 
            $kpp->family({ 
                    name    => $family,
                    members => $families{$family},
                    weights => $weights{$family},
            });
            $producers = $kpp->producing($family);
        } else {
            print "No family found for $family\n";
        }

        die "No producers found for $family\n" if (@$producers == 0);
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            next unless (defined $parent and $parent eq $VOC);
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            next if ($reaction_string eq "CO + OH = HO2"); 
            next if (exists $carbon_loss_rate{$reaction_string});
            my ($net_carbon) = get_total_C($reaction_string, $carbons, $kpp);
            next if ($net_carbon == 0);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $net_carbon * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($reactants, $products) = split / = /, $reaction_string;
            if ($rate->sum > 0) {
                $carbon_gain_rate{$reactants} += $rate(1:$NTIME-2);
            } else {
                $carbon_loss_rate{$reactants} += $rate(1:$NTIME-2);
            }
        }
    }

    my $others_max = 3e6;
    foreach my $reaction (keys %carbon_gain_rate) {
        if ($carbon_gain_rate{$reaction}->sum < $others_max) {
            $carbon_gain_rate{"Production Others"} += $carbon_gain_rate{$reaction};
            delete $carbon_gain_rate{$reaction};
        }
    }

    foreach my $reaction (keys %carbon_loss_rate) {
        if ($carbon_loss_rate{$reaction}->sum > -$others_max) {
            $carbon_loss_rate{"Consumption Others"} += $carbon_loss_rate{$reaction};
            delete $carbon_loss_rate{$reaction};
        }
    }

    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($carbon_gain_rate{$b}) <=> &$sort_function($carbon_gain_rate{$a}) } keys %carbon_gain_rate;
    my @cons_sorted_data = reverse sort { &$sort_function($carbon_loss_rate{$b}) <=> &$sort_function($carbon_loss_rate{$a}) } keys %carbon_loss_rate;
    
    my @final_sorted_data;
    foreach (@cons_sorted_data) { 
        next if ($_ eq 'Consumption Others') ;
        push @final_sorted_data, { $_ => $carbon_loss_rate{$_} };
    } 
    push @final_sorted_data, {'Consumption Others' => $carbon_loss_rate{'Consumption Others'}} if (defined $carbon_loss_rate{'Consumption Others'}) ;

    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $carbon_gain_rate{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $carbon_gain_rate{'Production Others'} } if (defined $carbon_gain_rate{'Production Others'}); 
    my (@plot_data, @legend_pos, @legend_neg, @legend);
    foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            if ($ref->{$item}->sum > 0) {
                push @legend_pos, $item;
            } else {
                push @legend_neg, $item;
            }
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    } 
    push @legend, reverse @legend_pos;
    push @legend, @legend_neg; 
    return (\@final_sorted_data, \@legend);
}

sub get_total_C {
    my ($reaction_string, $carbons, $kpp) = @_;
    my ($reactant_c, $product_c, @reactants, @products);

    my @inorganic = qw( hv OH HO2 O3 NO NO2 NO3 H2O HNO3 H2 PAROP O CO2 XO2 XO2N );
    my ($reactants, $products) = split / = /, $reaction_string;
    push @reactants, split / \+ /, $reactants;
    push @products, split / \+ /, $products;
    
    foreach my $reactant (@reactants) {
        next if ($reactant ~~ @inorganic);
        $reactant_c += get_species_carbon($reactant, $carbons);
    }
    
    return 0 unless (defined $reactant_c);
    foreach my $product (@products) {
        my ($yield, $item);
        if ($product =~ /^[0-9]|^\.[0-9]/) {
            ($yield, $item) = split / /, $product;
            next if ($item ~~ @inorganic);
            $product_c += $yield * get_species_carbon($item, $carbons);
        } else {
            next if ($product ~~ @inorganic);
            $product_c += get_species_carbon($product, $carbons);
        } 
    }
    $product_c = 0 unless (defined $product_c);
    return $product_c - $reactant_c;
}

sub get_species_carbon {
    my ($species, $carbons) = @_;
    my %carbons = %$carbons;
    my $carbon;
    if (exists $carbons{$species}) {
        $carbon = $carbons{$species};
    } else {
        print "No C found for species: $species\n";
    }
    return $carbon;
}

sub get_carbons {
    my ($run, $file) = @_;
    my $carbons;
    if ($run =~ /MCM/) {
        $carbons = mcm_n_carbon($file);
    } elsif ($run eq "MOZART_tagging") {
        $carbons = mozart_n_carbon($file);
    } elsif ($run =~ /CRI|RADM2|RACM|CB/) {
        $carbons = carbons_others($file);
    } else {
        print "$run doesn't match\n";
    }
    return $carbons;
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
