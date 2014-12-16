#! /usr/bin/env perl
# analysis of net rate of reactive carbon loss during Ox degradation in each mechanism, allocated to parent VOCs
# Version 0: Jane Coates 16/12/2014

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
#my @mechanisms = qw( CB05 ) ;
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
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    $n_carbon{"Ox_$mechanism"} = get_carbons($mechanism, $carbon_file);
    $data{$mechanism} = get_data($kpp, $mecca, $mechanism, $n_carbon{"Ox_$mechanism"});
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->run( q` my.colours = c("Carbon Loss Others" = "#696537", 
                           "Heptane" = "#f9c600", 
                           "Ethylbenzene" = "#76afca", 
                           "Benzene" = "#dc3522", 
                           "o-Xylene" = "#8c6238", 
                           "p-Xylene" = "#9bb08f", 
                           "Hexane" = "#8b1537", 
                           "2-Methylpropane" = "#e7e85e", 
                           "Propene" = "#0352cb", 
                           "Ethane" = "#86b650", 
                           "m-Xylene" = "#6c254f", 
                           "Isoprene" = "#ee6738", 
                           "Ethene" = "#58691b", 
                           "Pentane" = "#8ed6d5", 
                           "Propane" = "#f3aa7f", 
                           "Toluene" = "#c65d6c", 
                           "Butane" = "#888a87", 
                           "2-Methylbutane" = "#0e5c28", 
                           "Methane" = "#b569b3") `); 
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    $R->set('mechanism', $mechanism);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $VOC (sort keys %$ref) {
            $R->set('voc', $VOC);
            $R->set('loss.rate', [ map { $_ } $ref->{$VOC}->dog ]);
            $R->run(q` pre[voc] = loss.rate `);
        }
    }
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, VOC, Loss.Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `,
        q` data$VOC = factor(data$VOC, levels = c("Methane", "Ethane", "Propane", "Butane", "Pentane", "2-Methylbutane", "Ethene", "Propene", "Isoprene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Carbon Loss Others")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Loss.Rate, fill = VOC)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab(expression(bold(paste("Reactive Carbon Loss Rate (molecules ", cm^-3, s^-1, ")")))) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `, 
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 1e6)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0.3)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = levels(data$VOC), guide = guide_legend(nrow = 3)) `,
);

$R->run(q` CairoPDF(file = "Net_carbon_loss_VOC_allocated.pdf", width = 7, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $mechanism, $carbons) = @_;
    my %carbons = %$carbons;
    $families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

    my ($producers, %carbon_loss_rate);
    my @families = ("Ox_${mechanism}", "HO2x");
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
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent);
            my $reaction_string = $kpp->reaction_string($reaction);
            next if ($reaction_string eq "CO + OH = HO2"); 
            (my $un_tagged_reaction = $reaction_string) =~ s/_(.*?)\b//g;
            next if (exists $carbon_loss_rate{$reaction_string});
            my ($net_carbon) = get_total_C($un_tagged_reaction, $carbons, $kpp);
            next if ($net_carbon == 0);
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $net_carbon * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            $parent = get_chemical_name($parent);
            $carbon_loss_rate{$parent} += $rate(1:$NTIME-2);
        }
    } 
    
    foreach my $VOC (keys %carbon_loss_rate) {
        my $reshape = $carbon_loss_rate{$VOC}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $carbon_loss_rate{$VOC} = $integrate;
    }

    my $others = -4e7;
    foreach my $VOC (keys %carbon_loss_rate) {
        if ($carbon_loss_rate{$VOC}->sum > $others) {
            $carbon_loss_rate{"Carbon Loss Others"} += $carbon_loss_rate{$VOC};
            delete $carbon_loss_rate{$VOC};
        }
    }
    my $sort_function = sub { $_[0]->sum };

    my @sorted_loss = reverse sort { &$sort_function($carbon_loss_rate{$b}) <=> &$sort_function($carbon_loss_rate{$a}) } keys %carbon_loss_rate;
    my @sorted;
    foreach (@sorted_loss) {
        next if ($_ eq "Carbon Loss Others");
        push @sorted, { $_ => $carbon_loss_rate{$_} };
    }
    push @sorted, {'Carbon Loss Others' => $carbon_loss_rate{"Carbon Loss Others"} };
    return \@sorted;
}

sub get_total_C {
    my ($reaction_string, $carbons, $kpp) = @_;
    my ($reactant_c, $product_c, @reactants, @products);

    my @inorganic = qw( hv O2 H2O2 OH HO2 O3 NO NO2 NO3 H2O HNO3 H2 PAROP O CO2 XO2 XO2N OHOP );
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

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene";
    } elsif ($VOC eq "CH4") {
        $chemical_species = "Methane";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
