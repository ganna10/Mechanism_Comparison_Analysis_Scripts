#! /usr/bin/env perl
# HOx production and loss budgets in all mechanism
# Version 0: Jane Coates 18/11/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanism = ( "(a) MCM v3.2", "(b) MCM v3.1", "(c) CRI v2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( MOZART_tagging );
#my @mechanism = qw( MOZART-4 );
my $index = 0;
my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    $families{"HOx"} = [ qw( OH HO2 HO2NO2 )];

    my (%production, %consumption);
    $kpp->family({
            name    => "HOx",
            members => $families{"HOx"},
            weights => $weights{"HOx"},
    });
    my $producers = $kpp->producing("HOx");
    my $producer_yields = $kpp->effect_on("HOx", $producers);
    my $consumers = $kpp->consuming("HOx");
    my $consumer_yields = $kpp->effect_on("HOx", $consumers);
    print "No producers found in $run\n" if (@$producers == 0);
    print "No consumers found in $run\n" if (@$consumers == 0);

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
            $production{$reaction_string} += $rate(1:$NTIME-2);
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
            $consumption{$parent} += $rate(1:$NTIME-2);
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
            $consumption{$reaction_string} += $rate(1:$NTIME-2);
        }
    }
    remove_common_processes(\%production, \%consumption);

    my $others_max = 5e7;
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others_max) {
            $production{"Others"} += $production{$reaction};
            delete $production{$reaction};
        }
    }

    foreach my $reaction (keys %production) {
        my $reshape = $production{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production{$reaction} = $integrate;
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($production{$b}) <=> &$sort_function($production{$a}) } keys %production;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Others') ;
        push @final_sorted_data, { $_ => $production{$_} };
    } 
    push @final_sorted_data, { 'Others' => $production{'Others'} } if (defined $production{'Others'}); 
    $plot_data{$mechanism[$index]} = \@final_sorted_data; 
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
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
    $R->run(q` pre$Mechanism = rep(mechanism, rep(length(Time))) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "Reaction", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Reaction Rate\n") `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 180, face = "bold")) `,
        q` plot = plot + theme(axis.ticks.length = unit(2, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 150, angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(7, "cm")) `,
        q` plot = plot + theme(legend.text = element_text(size = 140)) `,
);

$R->run(q` CairoPDF(file = "HOx_budgets.pdf", width = 141, height = 200) `,
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
