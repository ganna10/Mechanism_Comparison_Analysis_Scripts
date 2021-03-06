#! /usr/bin/env perl
# analysis of net rate of reactive carbon loss during pentane and toluene degradation in each mechanism, all reactions not just Ox producing ones
# Version 0: Jane Coates 27/2/2015
# Version 1: Jane Coates 2/4/2015 including octane analysis in plot

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
#my @mechanisms = qw( MOZART-4 );
my (%n_carbon, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $n_carbon{$mechanism} = get_carbons($mechanism, $carbon_file);
    my @species = get_species($mechanism);
    foreach my $species (@species) {
        ($data{$mechanism}{$species}) = get_data($kpp, $mecca, $mechanism, $n_carbon{$mechanism}, $species);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(scales) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plot.data = data.frame() `);

my %final_data;
foreach my $mechanism (sort keys %data) {
    foreach my $species (sort keys %{$data{$mechanism}}) {
        if ($species =~ /NC5/) {
            $final_data{$mechanism}{"Pentane"} = $data{$mechanism}{$species};
        } elsif ($species eq "HC5") {
            $final_data{$mechanism}{"Pentane"} = $data{$mechanism}{$species} * 0.264;
        } elsif ($species =~ /NC8/) {
            $final_data{$mechanism}{"Octane"} = $data{$mechanism}{$species};
        } elsif ($species eq "HC8") {
            $final_data{$mechanism}{"Octane"} = $data{$mechanism}{$species};
        } elsif ($species eq "BIGALK") {
            $final_data{$mechanism}{"Pentane"} = $data{$mechanism}{$species} * 0.146;
            $final_data{$mechanism}{"Octane"} = $data{$mechanism}{$species} * 0.01;
        } elsif ($mechanism eq "RACM2" and $species eq "TOL") {
            $final_data{$mechanism}{"Toluene"} = $data{$mechanism}{$species} * 0.868;
        } elsif ($mechanism =~ /RA/ and $species eq "TOL") {
            $final_data{$mechanism}{"Toluene"} = $data{$mechanism}{$species} * 0.667;
        } elsif ($mechanism eq "MOZART-4" and $species =~ "TOL") {
            $final_data{$mechanism}{"Toluene"} = $data{$mechanism}{$species} * 0.478;
        } elsif ($species eq "TOLUENE") {
            $final_data{$mechanism}{"Toluene"} = $data{$mechanism}{$species};
        } 
    }
}

foreach my $run (sort keys %final_data) {
    $R->set('mechanism', $run);
    $R->run(q` data = data.frame(Time) `);
    #print "$run\n";
    foreach my $VOC (sort keys %{$final_data{$run}}) {
        #print "\t$VOC: $final_data{$run}{$VOC}\n";
        $R->set('voc', $VOC);
        $R->set('rate', [map { $_ } $final_data{$run}{$VOC}->dog]);
        $R->run(q` data[voc] = rate `);
    }
    $R->run(q` data$Mechanism = rep(mechanism, length(data$Time)) `,
            q` data = gather(data, VOC, Rate, -Mechanism, -Time) `,
            q` plot.data = rbind(plot.data, data) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` my.colours = c(  "CB05" = "#0352cb", "CBM-IV" = "#ef6638", "CRIv2" = "#b569b3", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot.data$Mechanism = factor(plot.data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `);
$R->run(q` plot.data$VOC = factor(plot.data$VOC, levels = c("Pentane", "Octane", "Toluene")) `);

$R->run(q` plot = ggplot(data = plot.data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + geom_point() `,
        q` plot = plot + facet_wrap( ~ VOC, scales = "free_y") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0.2)) `,
        q` plot = plot + scale_y_continuous(label = scientific_format() ) `,
        q` plot = plot + ylab("Rate of Change of Reactive Carbon (molecules cm-3 s-1)") `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "net_reactive_carbon_loss_pentane_toluene_octane.pdf", width = 8, height = 5.7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    my %carbon_loss_rate;
    my $all_reactions = $kpp->all_reactions(); #get all reactions and loop over them

    for (0..$#$all_reactions) { #get rates for all producing reactions
        my $reaction = $all_reactions->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent =~ $VOC);
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

    my @inorganic = qw( hv OH HO2 O3 NO NO2 NO3 H2O HNO3 H2 PAROP O CO2 XO2 XO2N OHOP UNITY CL SO3 SO2 H2O2 O2 HONO NULL );
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

sub get_species {
    my ($mechanism) = @_;

    my @species;
    if ($mechanism =~ /MCM|CRI|CB/) {
        @species = qw( NC5H12 NC8H18 TOLUENE );
    } elsif ($mechanism =~ /RA/) {
        @species = qw( HC5 HC8 TOL );
    } elsif ($mechanism =~ /MOZ/) {
        @species = qw( BIGALK TOLUENE);
    }
    return @species;
}
