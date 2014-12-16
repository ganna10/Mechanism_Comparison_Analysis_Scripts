#! /usr/bin/env perl
# Calculate total OxPE by normalising total Ox production by total Ox Consumption in each mechanism, allocated to parent VOCs
# Version 0: Jane Coates 16/12/2014

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
my $N_PER_DAY = 43200 / $DT;
my $NDAYS = int $NTIME / $N_PER_DAY;


my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
#my @mechanisms = qw( CB05 );
my (%families, %weights, %plot_data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanism"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3};
    ($plot_data{$mechanism}) = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    $R->set('mechanism', $run);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $VOC (sort keys %$ref) {
            $R->set('voc', $VOC);
            $R->set('oxpe', [map { $_ } $ref->{$VOC}->dog]);
            $R->run(q` pre[voc] = oxpe `);
        }
    }
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, VOC, OxPE, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run( q` my.colours = c("Others" = "#696537", 
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
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")) `,
        q` data$VOC = factor(data$VOC, levels = c("Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Ethene", "Propene", "Isoprene", "Toluene", "m-Xylene", "Others")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = OxPE, fill = VOC)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Ox Production Efficiency") `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `, 
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0.4)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0.3)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = levels(data$VOC), guide = guide_legend(nrow = 3)) `,
);

$R->run(q` CairoPDF(file = "Total_OxPE_VOC_allocated.pdf", width = 7, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $run) = @_;
    $families{"HO2x"} = [ qw(HO2 HO2NO2 )];
    my @loop = ("Ox_$run");
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
            my $reaction_string = $kpp->reaction_string($reaction);
            my ($number, $parent) = split /_/, $reaction;
            next unless defined $parent;
            $parent = get_chemical_name($parent);
            $production{$parent} += $rate(1:$NTIME-2);
        }

        foreach (0..$#$consumers) {
            last if ($species eq "HO2x");
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            $consumption += $rate(1:$NTIME-2);
        }
    }
    
    my $others = 1e8;
    foreach my $VOC (keys %production) {
        if ($production{$VOC}->sum < $others) {
            $production{"Others"} += $production{$VOC};
            delete $production{$VOC};
        }
    }
    
    my %OxPE;
    foreach my $VOC (keys %production) {
        my $OxPE = $production{$VOC} / -$consumption;
        my $reshape = $OxPE->reshape($N_PER_DAY, $NDAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $OxPE{$VOC} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_OxPE = sort { &$sort_function($OxPE{$b}) <=> &$sort_function($OxPE{$a}) } keys %OxPE;
    my @sorted;
    foreach (@sorted_OxPE) {
        next if ($_ eq "Others");
        push @sorted, { $_ => $OxPE{$_} };
    }
    push @sorted, { "Others" => $OxPE{"Others"} };
    return \@sorted;
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
