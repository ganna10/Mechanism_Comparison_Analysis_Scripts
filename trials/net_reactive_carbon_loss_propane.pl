#! /usr/bin/env perl
# analysis of net rate of reactive carbon loss during propane degradation in each mechanism
# Version 0: Jane Coates 4/1/2015

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
my (%n_carbon, %families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $ro2file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    $n_carbon{"Ox_$mechanism"} = get_carbons($mechanism, $carbon_file);
    my @VOCs = qw(Propane);
    foreach my $NMVOC (@VOCs) {
        my $parent = get_mechanism_species($NMVOC, $mechanism);
        ($data{$mechanism}{$NMVOC}) = get_data($kpp, $mecca, $mechanism, $n_carbon{"Ox_$mechanism"}, $parent);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plot.data = data.frame() `);
foreach my $run (sort keys %data) {
    $R->set('mechanism', $run);
    $R->run(q` data = data.frame(Time) `);
    foreach my $VOC (sort keys %{$data{$run}}) {
        $R->set('voc', $VOC);
        $R->set('rate', [map { $_ } $data{$run}{$VOC}->dog]);
        $R->run(q` data[voc] = rate `);
    }
    $R->run(q` data$Mechanism = rep(mechanism, length(data$Time)) `,
            q` data = gather(data, VOC, Rate, -Mechanism, -Time) `,
            q` plot.data = rbind(plot.data, data) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c(  "CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);

$R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `), #scientific label format for y-axis
$R->run(q` plot = ggplot(data = plot.data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ VOC) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab(expression(bold(paste("Net Carbon Loss Rate (molecules ", cm^-3, s^-1, ")")))) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(legend.justification = c(1, 0)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.position = c(1, 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "net_reactive_carbon_loss_propane.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    $families{"Ox_${mechanism}_${VOC}"} = $families{"Ox_$mechanism"};
    $families{"HO2x_${mechanism}_$VOC"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, %carbon_loss_rate);
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
            $carbon_loss_rate{$reaction_string} += $rate(1:$NTIME-2);
        }
    }

    foreach my $reaction (keys %carbon_loss_rate) {
        if ($mechanism =~ /CB/) {
            $carbon_loss_rate{$reaction} /= 8;
        } elsif ($mechanism =~ /MOZ/) {
            $carbon_loss_rate{$reaction} *= 0.010;
        }
    }

    my $overall_carbon_loss_rate = 0;
    $overall_carbon_loss_rate += $carbon_loss_rate{$_} foreach (keys %carbon_loss_rate);
    $overall_carbon_loss_rate = $overall_carbon_loss_rate->reshape($N_PER_DAY, $N_DAYS);
    $overall_carbon_loss_rate = $overall_carbon_loss_rate->sumover;
    $overall_carbon_loss_rate = $overall_carbon_loss_rate(0:13:2);
    return $overall_carbon_loss_rate;
}

sub get_total_C {
    my ($reaction_string, $carbons, $kpp) = @_;
    my ($reactant_c, $product_c, @reactants, @products);

    my @inorganic = qw( hv OH HO2 O3 NO NO2 NO3 H2O HNO3 H2 PAROP O CO2 XO2 XO2N OHOP );
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
    if ($run =~ /MCMv3\.1|MCMv3\.2/) {
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
    if ($NMVOC eq "Propane") {
        if ($run =~ /MCM|CRI|CB|MOZART/) {
            $mechanism_species = "C3H8";
        } elsif ($run =~ /RA/) {
            $mechanism_species = "HC3";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } else {
        print "No $NMVOC data\n";
    }
    return $mechanism_species;
}
