#! /usr/bin/env perl
# analysis of rate of reactive carbon loss during pentane degradation in each mechanism
# Version 0: Jane Coates 23/9/2014
# Version 1: Jane Coates 13/11/2014 refining plot for inclusion in supplement
# Version 2: Jane Coates 10/12/2014 refactoring code for constant emission runs
#####need to get script to print out plot#####

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2",  "CBM-IV", "CB05" );
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
my $NMVOC = "Pentane"; 
my (%n_carbon, %families, %weights, %data, %legend);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    my $ro2file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanism"} = get_carbons($mechanism, $carbon_file);
    my $parent = get_mechanism_species($NMVOC, $mechanism);
    ($data{$mechanism}, $legend{$mechanism}) = get_data($kpp, $mecca, $mechanism, $n_carbon{"Ox_$mechanism"}, $parent);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
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

$R->run(q` plotting = function (data, mechanism) {  plot = ggplot(data, aes(x = Time, y = Carbon.loss.rate, fill = Reaction)) ;
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
                                                            plot = plot + scale_fill_manual(values = my.colours) ;
                                                            return(plot) } `);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $run (sort keys %data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('mechanism', $run);
    #$R->set('legend', [@{$legend{$run}}]);
    $R->run(q` data = gather(data, Reaction, Carbon.loss.rate, -Time) `,
            q` reaction.levels = rev(levels(factor(data$Reaction))) `,
            q` data$Reaction = ordered(data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(data, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->stop();

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

    my $others_max = 3e7;
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

    foreach my $reaction (keys %carbon_gain_rate) {
        my $reshape = $carbon_gain_rate{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $carbon_gain_rate{$reaction} = $integrate;
    }

    foreach my $reaction (keys %carbon_loss_rate) {
        my $reshape = $carbon_loss_rate{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $carbon_loss_rate{$reaction} = $integrate;
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
    return \@final_sorted_data;
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
    } elsif ($run =~ /MOZART/) {
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
