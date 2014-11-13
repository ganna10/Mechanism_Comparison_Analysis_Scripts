#! /usr/bin/env perl
# analysis of rate of reactive carbon loss during pentane degradation in each mechanism
# Version 0: Jane Coates 23/9/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
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
    my $ro2file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanisms[$array_index]"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$array_index]"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanisms[$array_index]"} = get_carbons($run);
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
                            "Consumption Others" = "#1c3e3d",
                            "C2O3 + NO = MEO2 + NO2" = "#8c1531", "C2O3 + NO = HCHO + HO2 + NO2 + XO2" = "#8c1531", "CH3CO3 + NO = CH3O2 + NO2" = "#8c1531", "CH3CO3 + NO = CH3O2 + CO2 + NO2" = "#8c1531", "ACO3 + NO = MO2 + NO2" = "#8c1531",
                            "ROR = HO2" = "#6c254f",
                            "ROR = 0.6 ALD2 + 0.5 ALDX + 0.94 HO2\n+ 2.1 + PAROP + 0.02 ROR + 0.96 XO2 + 0.04 XO2N" = "#8ed6d2", "ROR = 1.1 ALD2 + 0.94 HO2 + 2.1 PAROP\n+ 0.02 ROR + 0.96 XO2 + 0.04 XO2N" = "#8ed6d2",
                            "C2H5CO3 + NO = C2H5O2 + NO2" = "#f3aa7f",
                            "CO2C4CO3 + NO = MEKAO2 + NO2" = "#0e5c28",
                            "HOC2H4CO3 + NO = HOCH2CH2O2 + NO2" = "#f9c500", 
                            "CO2C3CO3 + NO = CH3COCH2O2 + NO2" = "#898989",
                            "CH3CO3 = .3 CH2O + .1 CH3COOH + .9 CH3O2\n+ .9 CO2 + .5 HO2" = "#a67c52",
                            "NO + RCO3 = ETHP + NO2" = "#86b650",
                            "ALKO2 + NO = .1 CH2O + .4 CH3CHO\n+ .25 CH3COCH3 + .9 HO2 + .75 MEK\n+ .9 NO2 + .1 ONIT" = "#77aecc", 
                            "ETHP + NO = ALD + HO2 + NO2" = "#f7c56c",
                            "HC3P + NO = .504 ACD + .165 ACT + .132 ALD\n+ .048 ETHP + .660 HO2 + .042 MEK + .131 MO2\n+ .935 NO2 + .065 ONIT + .089 XO2" = "#dc3522", "HC3P + NO = 0.233 ALD + 0.048 ETHP + 0.063 GLY\n+ 0.047 HCHO + 0.742 HO2 + 0.623 KET + 0.15 MO2\n+ 0.941 NO2 + 0.059 ONIT + 0.048 XO2" = "#dc3522",
                            "HC5 + OH = H2O + 0.75 HC5P\n+ 0.25 HO2 + 0.25 KET" = "#9bb18d",
                            "HC5P + NO = 0.38 ALD + 0.92 HO2\n+ 0.69 KET + 0.92 NO2 + 0.08 ONIT" = "#c9a415", "HC5P + NO = .045 ACD + .217 ACT\n+ .203 ALD + .231 ETHP + .018 HCHO + .272 HKET\n+ .2 HO2 + .039 KET + .033 MEK + .051 MO2\n+ .864 NO2 + .136 ONIT + .235 XO2" = "#c9a415", "HC5P + NO = 0.211 ALD + 0.245 ETHP + 0.021 HCHO\n+ 0.599 HO2 + 0.722 KET + 0.031 MO2\n+ 0.876 NO2 + 0.124 ONIT + 0.334 XO2" = "#c9a415",
                            "KETP + NO = HO2 + MGLY + NO" = "#ef6638", "KETP + NO = 0.23 ACO3 + 0.46 ALD + 0.77 HO2\n+ 0.54 MGLY + NO2 + 0.16 XO2" = "#ef6638",
                            "OH + ONIT = H2O + HC3P + NO2" = "#58691b") `,
);

$R->run(q` plotting = function (data, mechanism, legend) {  plot = ggplot(data, aes(x = Time, y = Carbon.loss.rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(data = subset(data, data$Carbon.loss.rate < 0), stat = "identity") ;
                                                            plot = plot + geom_bar(data = subset(data, data$Carbon.loss.rate > 0), stat = "identity") ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank());
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.0)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.0)) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 80)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 80)) ;
                                                            plot = plot + theme(plot.title = element_text(size = 140, face = "bold")) ;
                                                            plot = plot + theme(legend.text = element_text(size = 67)) ;
                                                            plot = plot + theme(legend.key.size = unit(6.5, "cm")) ;
                                                            plot = plot + scale_y_continuous(limits = c(-2.3e8, 2e7), breaks = seq(-2.3e8, 2e7, 1e7)) ;
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
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
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
        
        my $prod_others_max = 3e6;
        my $max_string_width = 27;
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
            $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
            $reaction_string =~ s/\n$//;
            $reaction_string = string_mapping($reaction_string);
            if ($rate->sum < $prod_others_max and $rate->sum > 0) {
                $carbon_gain_rate{"Production Others"} += $rate(1:$NTIME-2);
            } elsif ($rate->sum > -$prod_others_max and $rate->sum < 0) {
                $carbon_loss_rate{"Consumption Others"} += $rate(1:$NTIME-2);
            } elsif ($rate->sum > 0) {
                $carbon_gain_rate{$reaction_string} += $rate(1:$NTIME-2);
            } else {
                $carbon_loss_rate{$reaction_string} += $rate(1:$NTIME-2);
            }
        }
    }

    sub string_mapping {
        my ($string) = @_;
        if ($string =~ /ROR = 0\.600 ALD2/) { 
            $string = "ROR = 0.6 ALD2 + 0.5 ALDX + 0.94 HO2\n+ 2.1 + PAROP + 0.02 ROR + 0.96 XO2 + 0.04 XO2N";
        } elsif ($string =~ /ROR = 1\.1 ALD2/) {
            $string = "ROR = 1.1 ALD2 + 0.94 HO2 + 2.1 PAROP\n+ 0.02 ROR + 0.96 XO2 + 0.04 XO2N";
        } elsif ($string =~ /HC3P \+ NO = \.504 ACD/) {
            $string = "HC3P + NO = .504 ACD + .165 ACT + .132 ALD\n+ .048 ETHP + .660 HO2 + .042 MEK + .131 MO2\n+ .935 NO2 + .065 ONIT + .089 XO2";
        } elsif ($string =~ /C2O3 \+ NO = HCHO/) {
            $string = "C2O3 + NO = HCHO + HO2 + NO2 + XO2";
        } elsif ($string =~ /CO2C4CO3/) {
            $string = "CO2C4CO3 + NO = MEKAO2 + NO2";
        } elsif ($string =~ /HOC2H4CO3/) {
            $string = "HOC2H4CO3 + NO = HOCH2CH2O2 + NO2";
        } elsif ($string =~ /CO2C3CO3/) {
            $string = "CO2C3CO3 + NO = CH3COCH2O2 + NO2";
        } elsif ($string =~  /CH3CO3 = \.3/) {
            $string = "CH3CO3 = .3 CH2O + .1 CH3COOH + .9 CH3O2\n+ .9 CO2 + .5 HO2";
        } elsif ($string =~ /ALKO2 \+ NO/) {
            $string = "ALKO2 + NO = .1 CH2O + .4 CH3CHO\n+ .25 CH3COCH3 + .9 HO2 + .75 MEK\n+ .9 NO2 + .1 ONIT";
        } elsif ($string =~ /CH3CO3 \+ NO = CH3O2 \+ CO2/) {
            $string = "CH3CO3 + NO = CH3O2 + CO2 + NO2";
        } elsif ($string =~ /HC5P \+ NO = 0\.38 ALD/) {
            $string = "HC5P + NO = 0.38 ALD + 0.92 HO2\n+ 0.69 KET + 0.92 NO2 + 0.08 ONIT";
        } elsif ($string =~ /KETP \+ NO = HO2/) {
            $string = "KETP + NO = HO2 + MGLY + NO";
        } elsif ($string =~ /HC5P \+ NO = \.045 ACD/) {
            $string = "HC5P + NO = .045 ACD + .217 ACT\n+ .203 ALD + .231 ETHP + .018 HCHO + .272 HKET\n+ .2 HO2 + .039 KET + .033 MEK + .051 MO2\n+ .864 NO2 + .136 ONIT + .235 XO2";
        } elsif ($string =~ /KETP \+ NO = 0\.23 ACO3/) {
            $string = "KETP + NO = 0.23 ACO3 + 0.46 ALD + 0.77 HO2\n+ 0.54 MGLY + NO2 + 0.16 XO2";
        } elsif ($string =~ /HC3P \+ NO = 0\.233 ALD/) {
            $string = "HC3P + NO = 0.233 ALD + 0.048 ETHP + 0.063 GLY\n+ 0.047 HCHO + 0.742 HO2 + 0.623 KET + 0.15 MO2\n+ 0.941 NO2 + 0.059 ONIT + 0.048 XO2";
        } elsif ($string =~ /HC5 \+ OH/) {
            $string = "HC5 + OH = H2O + 0.75 HC5P\n+ 0.25 HO2 + 0.25 KET";
        } elsif ($string =~ /HC5P \+ NO = 0\.211 ALD/) {
            $string = "HC5P + NO = 0.211 ALD + 0.245 ETHP + 0.021 HCHO\n+ 0.599 HO2 + 0.722 KET + 0.031 MO2\n+ 0.876 NO2 + 0.124 ONIT + 0.334 XO2";
        } elsif ($string =~ /OH \+ ONIT/) {
            $string = "OH + ONIT = H2O + HC3P + NO2";
        }
        return $string;
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
