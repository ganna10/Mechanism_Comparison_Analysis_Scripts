#!/usr/bin/env perl 
# Ox production budget plots from all tagged mechanisms, lumped species are de-lumped into constituent VOC. HO2x not assigned to verify if this is why Ox is so low in CBs compared to MCM
# Version 0: Jane Coates 2/2/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base_dir = "/local/home/coates/MECCA/";
my $mecca = MECCA->new("$base_dir/CB05_tagged/boxmodel"); 
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my (%data, %families, %weights);
my @mechanisms = qw( MCMv3.2 CBM-IV CB05 ); 
#my @mechanisms = qw( CB05 );

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base_dir/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base_dir/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $ro2file = "$base_dir/${mechanism}_tagged/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
    $families{"Ox_$mechanism"} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ]; 
    $weights{"Ox_$mechanism"} = { NO3 => 2, N2O5 => 3 };
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(dplyr) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $process (sort keys %$ref) {
            my $voc = get_chemical_name($process);
            $R->set('voc', $voc);
            $R->set('rate', [ map { $_ } $ref->{$process}->dog ]);
            $R->run(q` pre[voc] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, VOC, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "CBM-IV", "CB05")) `,
        q` VOC.levels = c( "HO2", "Methane ", "Propane ", "Butane ", "Pentane ", "2-Methylbutane ", "Isoprene ", "Others" ) `, 
        q` my.colours = c(  "Others" = "#696537", 
                            "HO2" = "#6c254f",
                            "Methane " = "#f9c500", 
                            "Ethane " = "#0352cb", 
                            "Propane " = "#0e5c28", 
                            "2-Methylpropane " = "#e7e85e", 
                            "Butane " = "#ef6638", 
                            "Pentane " = "#b569b3", 
                            "2-Methylbutane " = "#4c9383", 
                            "Hexane " = "#86b650", 
                            "Ethene " = "#cc6329", 
                            "Propene " = "#2b9eb3", 
                            "2-Methylpropene " = "#f7c56c",
                            "Isoprene " = "#0c3f78",
                            "Toluene " = "#8c1531",
                            "m-Xylene " = "#6db875",
                            "o-Xylene " = "#f3aa7f",
                            "p-Xylene " = "#be2448" ) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = VOC )) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Reaction Rate (molecules cm-3 s-1)") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 0.4)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(VOC.levels)) `,
);

$R->run(q` CairoPDF( file = "Ox_production_budgets_by_VOC_de-allocated_CBs_MCM_not_HO2.pdf" , width = 8.5, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    my (%production, %consumption);

    my ($producers, $producer_yields, $consumers, $consumer_yields);
    if (exists $families{"Ox_${mechanism}"}) {
        $kpp->family({
                name    => "Ox_${mechanism}",
                members => $families{"Ox_${mechanism}"},
                weights => $weights{"Ox_${mechanism}"},
        });
        $producers = $kpp->producing("Ox_${mechanism}");
        $producer_yields = $kpp->effect_on("Ox_${mechanism}", $producers);
        $consumers = $kpp->consuming("Ox_${mechanism}");
        $consumer_yields = $kpp->effect_on("Ox_${mechanism}", $consumers);
    } else {
        print "No family for ${mechanism}\n";
    }
    print "No producers found for ${mechanism}\n" if (@$producers == 0);
    print "No consumers found for ${mechanism}\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $production{$parent} += $rate(1:$NTIME-2);
        } else {
            my $string = $kpp->reaction_string($reaction);
            $string = "CO " if ($string =~ /CO \+ OH|OH \+ CO/);
            $production{$string} += $rate(1:$NTIME-2);
        }
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $consumption{$parent} += $rate(1:$NTIME-2);
        } else {
            my $string = $kpp->reaction_string($reaction);
            $string = "CO " if ($string =~ /CO \+ OH|OH \+ CO/);
            $consumption{$string} += $rate(1:$NTIME-2);
        }
    }

    remove_common_processes(\%production, \%consumption);

    foreach my $process (keys %production) {
        my $reshape = $production{$process}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$process} = $integrate;
    }

    my $others = 1.5e8;
    foreach my $process (keys %production) {
        if ($production{$process}->sum < $others) {
            $production{"Others"} += $production{$process};
            delete $production{$process};
        }
    }

    my $sort_function = sub { $_[0]->sum } ;
    my @sorted_prod = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    my @sorted_data;
    foreach (@sorted_prod) {
        next if ($_ eq 'Others' or $_ eq 'Methane ' or $_ eq 'CO ');
        push @sorted_data, { $_ => $production{$_} }
    }
    push @sorted_data, { 'Others' => $production{'Others'} } if (defined $production{'Others'});
    unshift @sorted_data, { 'Methane ' => $production{'Methane '} } if (defined $production{'Methane '});
    return \@sorted_data;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
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

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'CO ') {
        $chemical_species = 'CO ';
    } elsif ($VOC eq 'CH4') {
        $chemical_species = 'Methane ';
    } elsif ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane ';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane ';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane ';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane ';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane ';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane ';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane ';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane ";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane ";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene ';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene ';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene ';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene ";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene ";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene ';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene ";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene ';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene ";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene ";
    } elsif ($VOC eq 'Others') {
        $chemical_species = 'Others';
    } elsif ($VOC =~ /HO2/) {
        $chemical_species = "HO2";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
