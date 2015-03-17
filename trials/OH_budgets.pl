#! /usr/bin/env perl
# OH production and loss budgets in all mechanism
# Version 0: Jane Coates 2/2/2015
# Version 1: Jane Coates 17/3/2015 including all mechanisms

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

#my @mechanisms = ( "CB05" );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my (%families, %weights, %plot_data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    $families{"HOx"} = [ qw( OH )];
    #$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
    $kpp->family({
            name    => "HOx",
            members => $families{"HOx"},
            weights => $weights{"HOx"},
    });

    my (%production, %consumption);
    #$kpp->family({
    #name    => "HOx",
    #members => $families{"HOx"},
    #weights => $weights{"HOx"},
    #});
    my $producers = $kpp->producing("OH");
    my $producer_yields = $kpp->effect_on("OH", $producers);
    my $consumers = $kpp->consuming("OH");
    my $consumer_yields = $kpp->effect_on("OH", $consumers);
    print "No producers found in $mechanism\n" if (@$producers == 0);
    print "No consumers found in $mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            my $name = get_name($parent);
            $production{$name} += $rate(1:$NTIME-2);
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
#            if ($reaction_string =~ /HO2 \+ NO/) {
#                my $ho2x_producers = $kpp->producing("HO2x");
#                my $ho2x_producer_yields = $kpp->effect_on("HO2x", $ho2x_producers);
#                foreach (0..$#$ho2x_producers) {
#                    my $reaction = $ho2x_producers->[$_];
#                    next unless ($reaction =~ /_/);
#                    my ($number, $parent) = split /_/, $reaction;
#                    my $name = get_name($parent);
#                    my $reaction_number = $kpp->reaction_number($reaction);
#                    my $rate = $mecca->rate($reaction_number) * $ho2x_producer_yields->[$_];
#                    next if ($rate->sum == 0);
            #$production{$name} += $rate(1:$NTIME-2);
                    #}
                    #} else {
                $production{$reaction_string} += $rate(1:$NTIME-2);
                #}
        }
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            my $name = get_name($parent);
            $consumption{$name} += $rate(1:$NTIME-2);
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
            $consumption{$reaction_string} += $rate(1:$NTIME-2);
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others_max = 2e8;
    foreach my $reaction (keys %production) {
        my $reshape = $production{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        if ($production{$reaction}->sum < $others_max) {
            $production{"Production Others"} += $integrate(0:13:2);
            delete $production{$reaction};
        } else {
            $production{$reaction} = $integrate(0:13:2);
        }
    }
    
    foreach my $reaction (keys %consumption) {
        my $reshape = $consumption{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        if ($consumption{$reaction}->sum > -$others_max) {
            $consumption{"Consumption Others"} += $integrate(0:13:2);
            delete $consumption{$reaction};
        } else {
            $consumption{$reaction} = $integrate(0:13:2);
        }
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    my @sorted_cons = reverse sort { &$sort_function($consumption{$b}) <=> &$sort_function($consumption{$a}) } keys %consumption;
    
    my @final_sorted_data;
    foreach (@sorted_cons) { 
        next if ($_ =~ /Others/) ;
        push @final_sorted_data, { $_ => $consumption{$_} };
    } 
    push @final_sorted_data, { 'Consumption Others' => $consumption{'Consumption Others'} } if (defined $consumption{'Consumption Others'}); 
    foreach (@sorted_prod) { 
        next if ($_ =~ /Others/) ;
        push @final_sorted_data, { $_ => $production{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $production{'Production Others'} } if (defined $production{'Production Others'}); 
    $plot_data{$mechanism} = \@final_sorted_data; 
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->set('filename', "${mechanism}_OH_budget.csv");
    $R->run(q` pre$Mechanism = rep(mechanism, rep(length(Time))) `,
            q` write.csv(pre, file = filename, col.names = FALSE, row.names = TRUE) `,
            q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` my.colours = c(  "Production Others" = "#696537", 
                            "Consumption Others" = "#ee6738",
                            "Methane" = "#f9c500",
                            "O1D = OH + OH" = "#6c254f", 
                            "2-Methylbutane" = "#0352cb", 
                            "Butane" = "#e7e85e", 
                            "HONO + hv = NO + OH" = "#0e5c28", 
                            "NO2 + OH = HNO3" = "#2b9eb3", 
                            "CO + OH = HO2" = "#0c3f78", 
                            "NO + OH = HONO" = "#b569b3", 
                            "HO2 + NO = NO2 + OH" = "#000000",
                            "Toluene" = "#898989",
                            "Ethene" = "#77aecc",
                            "Isoprene" = "#86b650",
                            "Propane" = "#cc6329" ) `); 
$R->run(#q` data$Reaction = factor(data$Reaction, levels = c("Production Others", "Ethene", "Isoprene", "Methane", "HONO + hv = NO + OH", "O1D = OH + OH", "NO2 + OH = HNO3", "NO + OH = HONO", "HO2 + HO2 = H2O2", "HO2 + OH = UNITY", "Consumption Others")) `,
    #q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "CBM-IV", "CB05")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(data = filter(data, Rate < 0), stat = "identity") `,
        q` plot = plot + geom_bar(data = filter(data, Rate > 0), stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Reaction Rate (molecules cm-3 s-1)") `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 0.7, hjust = 0.8)) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        #q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "OH_budgets.pdf", width = 8, height = 5.7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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

sub get_name {
    my ($parent) = @_;
    if ($parent eq "CH4") {
        $parent = "Methane";
    } elsif ($parent =~ /C2H6|ETH/) {
        $parent = "Ethane";
    } elsif ($parent =~ /C3H8|HC3/) {
        $parent = "Propane";
    } elsif ($parent eq "NC4H10") {
        $parent = "Butane";
    } elsif ($parent eq "IC4H10") {
        $parent = "2-Methylpropane";
    } elsif ($parent =~ /NC5H12|BIGALK|HC5/) {
        $parent = "Pentane";
    } elsif ($parent eq "IC5H12") {
        $parent = "2-Methylbutane";
    } elsif ($parent eq "NC6H14") {
        $parent = "Hexane";
    } elsif ($parent eq "NC7H16") {
        $parent = "Heptane";
    } elsif ($parent =~ /NC8H18|HC8/) {
        $parent = "Octane";
    } elsif ($parent =~ /C2H4|OL2|ETE/) {
        $parent = "Ethene";
    } elsif ($parent =~ /C3H6|OLT/) {
        $parent = "Propene";
    } elsif ($parent =~ /BUT1ENE|BIGENE/) {
        $parent = "Butene";
    } elsif ($parent =~ /MEPROPENE|OLI/) {
        $parent = "2-Methylpropene";
    } elsif ($parent =~ /C5H8|ISO/) {
        $parent = "Isoprene";
    } elsif ($parent =~ /^BEN/) {
        $parent = "Benzene";
    } elsif ($parent =~ /TOL/) {
        $parent = "Toluene";
    } elsif ($parent =~ /OXYL|XYO/) {
        $parent = "o-Xylene";
    } elsif ($parent =~ /PXYL|XYP/) {
        $parent = "p-Xylene";
    } elsif ($parent eq "EBENZ") {
        $parent = "Ethylbenzene";
    } elsif ($parent =~ /MXYL|XYL|XYM/) {
        $parent = "m-Xylene";
    } elsif ($parent =~ /Others/) {
        $parent = $parent;
    } else {
        print "No chemical name for $parent\n";
    }
    return $parent;
}
