#! /usr/bin/env perl
# Show reactions contributing to pentane C5 Ox production on first 2 days in each mechanism (except CBs)
# Version 0: Jane Coates 6/1/2015

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
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT ;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2" );
#my @mechanisms = qw( RACM2 );
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
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %plot_data) {
    foreach my $VOC (sort keys %{$plot_data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $ref (@{$plot_data{$mechanism}{$VOC}}) {
            foreach my $reaction (sort keys %$ref) {
                my $name;
                if ($reaction =~ /NO \+ PEBO2|ALKO2|HC5P|RN16O2/) {
                    $name = "NO + Pentyl Peroxy Radical";
                } else {
                    $name = $reaction;
                }
                $R->set('reaction', $name);
                $R->set('production', [ map { $_ } $ref->{$reaction}->dog ]);
                $R->run(q` pre[reaction] = production `);
            }
        }
        $R->set('mechanism', $mechanism);
        $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre = gather(pre, Reaction, Production, -Time, -Mechanism) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(pre) `);
#print "$p\n";
$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "NO + Pentyl Peroxy Radical" = "#6c254f",
                            "NO + RN15AO2" = "#f9c500", 
                            "KETP + NO" = "#0e5c28", 
                            "NO + PECO2" = "#e7e85e", 
                            "PECO" = "#b569b3", 
                            "HO2C5O" = "#ef6638", 
                            "HO2C5O2 + NO" = "#4c9383" ) `, 
);
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4")) `);
$R->run(q` data$Reaction = factor(data$Reaction, levels = c("NO + Pentyl Peroxy Radical", "HO2C5O", "HO2C5O2 + NO", "NO + PECO2", "PECO", "NO + RN15AO2", "KETP + NO", "Production Others")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Production, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + ylab("Reaction Rates (molecules cm-3 s-1)") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Reaction))) `,
        q` plot = plot + theme(legend.position = c(0.7, 0.17)) `,
);

$R->run(q` CairoPDF(file = "C5_pentane_reaction_rates.pdf", width = 6, height = 9) `,
        q` print(plot) `,
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
    my (%production, %consumption);
    my @c5s = qw( 5 4.5 4.8 5.6 );
    my $nr = 0;
    my $yield = 0;
    print "$mechanism\n";

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

        foreach my $i (0..$#$producers) {
            my $reaction = $producers->[$i];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$i];
            next if ($rate->sum == 0);
            my ($r_number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $reactants = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $mechanism =~ /RADM2|RACM|CB/);
                    my ($lookup, $tag) = split /_/, $_;
                    if (defined $carbons{$lookup}) {
                        next unless ($carbons{$lookup} ~~ @c5s);
                        $nr++;
                        $yield += $producer_yields->[$i];
                        my ($reaction_string) = $kpp->reaction_string($reaction);
                        $reaction_string =~ s/_(.*?)\b//g;
                        my ($label, $rest) = split / = /, $reaction_string;
                        $production{$label} += $rate(1:$NTIME-2);
                    } else {
                        print "Nothing found for $lookup\n";
                    }
                }
            }
        }

        if ($species =~ /Ox/ and $mechanism =~ /RADM2|RACM|CB/) {#operator allocation for those mechanisms that use it: RADM2, RACM, RACM2
            my $operator = "XO2_" . $VOC;
            my $op_producers = $kpp->producing($operator);
            my $op_producer_yields = $kpp->effect_on($operator, $op_producers); 
            die "No producers found for $operator\n" if (@$op_producers == 0);

            for my $i (0..$#$op_producers) { #get rates for all producing reactions
                my $reaction = $op_producers->[$i];
                my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
                my $reaction_number = $kpp->reaction_number($reaction);
                my $rate = $op_producer_yields->[$i] * $mecca->rate($reaction_number); 
                next if ($rate->sum == 0); # do not include reactions that do not occur 
                my ($reactants) = $kpp->reactants($reaction);
                foreach (@$reactants) {
                    if ($_ =~ /_/) {
                        my ($lookup, $rest) = split '_', $_;
                        if (defined $carbons{$lookup}) {
                            next unless ($carbons{$lookup} ~~ @c5s);
                            $nr++;
                            $yield += $op_producer_yields->[$i];
                            my ($reaction_string) = $kpp->reaction_string($reaction);
                            $reaction_string =~ s/_(.*?)\b//g;
                            my ($label, $rest) = split / = /, $reaction_string;
                            $production{$label} += $rate(1:$NTIME-2);
                        } else {
                            print "$mechanism => nothing found for $lookup production\n";
                        }
                    }
                } 
            } 
        } 

        print "$species\n";
        print "number of C5s => $nr\n" ;
        print "total Ox yield of C5s => $yield\n" ;
        foreach (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($r_number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq $VOC);
            my $reactants = $kpp->reactants($reaction);
            foreach (@$reactants) {
                if ($_ =~ /_/) {
                    next if ($_ =~ /XO2/ and $mechanism =~ /RADM2|RACM|CB/);
                    my ($lookup, $tag) = split /_/, $_;
                    if (defined $carbons{$lookup}) {
                        next unless ($carbons{$lookup} ~~ @c5s);
                        my ($reaction_string) = $kpp->reaction_string($reaction);
                        $reaction_string =~ s/_(.*?)\b//g;
                        my ($label, $rest) = split / = /, $reaction_string;
                        $consumption{$label} += $rate(1:$NTIME-2);
                    } else {
                        print "Nothing found for $lookup\n";
                    }
                }
            }
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others = 6e6;
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others) {
            $production{"Production Others"} += $production{$reaction};
            delete $production{$reaction};
        }
    }

    my $parent;
    if ($mechanism =~ /CB/) {
        $parent = "PAR_NC5H12";
        } elsif ($mechanism =~ /RA/) {
            $parent = "HC5";
        } elsif ($mechanism eq "MOZART-4") {
            $parent = "BIGALK";
    } else {
        $parent = "NC5H12";
    }
    my $emission_reaction = $kpp->producing_from($parent, "UNITY");
    my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
    my $emission_rate = $mecca->rate($reaction_number); 
    $emission_rate = $emission_rate(1:$NTIME-2);
    $emission_rate = $emission_rate->sum * $DT; 
    $emission_rate /= 5 if ($mechanism =~ /CB/);
    
    #normalise by dividing reaction rate of intermediate (molecules (intermediate) /cm3/s) by number density of parent VOC (molecules (VOC) /cm3)
    $production{$_} /= $emission_rate foreach (sort keys %production); 

    foreach my $reaction (keys %production) {
        if ($mechanism =~ /MOZ/) {
            $production{$reaction} *= 0.146;
        } elsif ($mechanism =~ /RA/) {
            $production{$reaction} *= 0.264;
        }
        my $reshape = $production{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:3:2);
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
