#! /usr/bin/env perl
# Ox Consumption budgets attributed to VOCs and inorganic reactions
# Version 0: Jane Coates 13/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT;
my $N_DAYS = int( $NTIME / $N_PER_DAY );

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCM v3.2", "(b) MCM v3.1", "(c) CRI v2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my $index = 0;
my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    $families{"Ox_$mechanisms[$index]"} = [ qw( O3 NO2 O O1D HO2NO2 NO3 N2O5 ), @no2_reservoirs ];
    $weights{"Ox_$mechanisms[$index]"} = { NO3 => 2, N2O5 => 3 };
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, "Ox_$mechanisms[$index]");
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(grid) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);

foreach my $run (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $process (sort keys %$ref) {
            $R->set('process', $process);
            $R->set('rate', [ map { $_ } $ref->{$process}->dog ]);
            $R->run(q` pre[process] = rate `);
        }
    }
    $R->set('mechanism', $run);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "Process", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` my.names = c("TOL" = "TOLUENE", "BEN" = "BENZENE", "HC3" = "C3H8", "XYL" = "MXYL", "HC5" = "NC5H12", "OL2" = "C2H4", "OLI" = "MEPROPENE", "XYM" = "MXYL", "XYP" = "PXYL", "XYO" = "OXYL", "OLT" = "C3H6", "ETE" = "C2H4", "ISO" = "C5H8", "BIGALK"  = "NC5H12") `,
        q` my.colours = c("C2H6" = "#696537", "HC3" = "#f9c600", "C3H8" = "#f9c600", "NC4H10" = "#76afca", "IC4H10" = "#dc3522", "BIGALK" = "#8c6238", "HC5" = "#8c6238", "NC5H12" = "#8c6238", "IC5H12" = "#9bb08f", "NC6H14" = "#8b1537", "NC7H16" = "#ba8b01", "NC8H18" = "#0352cb", "OL2" = "#86b650", "ETE" = "#86b650", "C2H4" = "#86b650", "OLT" = "#6c254f", "C3H6" = "#6c254f", "BUT1ENE" = "#ee6738", "OLI" = "#58691b", "MEPROPENE" = "#58691b", "ISO" = "#8ed6d5", "C5H8" = "#8ed6d5", "BEN" = "#f3aaf7", "BENZENE" = "#f3aa7f", "TOL" = "#c65d6c", "TOLUENE" = "#c65d6c", "XYL" = "#888a87", "XYM" = "#888a87", "MXYL" = "#888a87", "XYO" = "#0e5c28", "OXYL" = "#0e5c28", "XYP" = "#b569b3", "PXYL" = "#b569b3", "EBENZ" = "#2c9def", "CH4" = "#000000", "Consumption Others" = "#898989" ) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Reaction Rate") `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 150)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(axis.ticks.length = unit(2, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(legend.text = element_text(size = 140)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "Ox_consumption_budgets_tagged.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $species) = @_;
    
    my %consumption_rates;
    my ($consumers, $consumer_yields);
    if (exists $families{$species}) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $consumers = $kpp->consuming($species);
        $consumer_yields = $kpp->effect_on($species, $consumers);
    } else {
        print "No family found for $species\n";
    }
    print "No consumers found for $species\n" if (@$consumers == 0);

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $consumption_rates{$parent} += $rate(1:$NTIME-2);
        } else {
            next;
            my $reaction_string =  $kpp->reaction_string($reaction);
            $consumption_rates{"Inorganic"} += $rate(1:$NTIME-2);
        }
    }
    
    my $others_max = -1e7;
    foreach my $process (keys %consumption_rates) {
        if ($consumption_rates{$process}->sum > $others_max) {
            $consumption_rates{"Consumption Others"} += $consumption_rates{$process};
            delete $consumption_rates{$process};
        }
    }

    foreach my $process (keys %consumption_rates) {
        my $reshape = $consumption_rates{$process}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $consumption_rates{$process} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($consumption_rates{$a}) <=> &$sort_function($consumption_rates{$b}) } keys %consumption_rates;
    my @final_sort;
    foreach (@sorted_data) {
        next if ($_ eq "Consumption Others");
        push @final_sort, { $_ => $consumption_rates{$_} };
    }
    push @final_sort, { "Consumption Others" => $consumption_rates{"Consumption Others"} } if (defined $consumption_rates{"Consumption Others"});

    return \@final_sort;
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
