#! /usr/bin/env perl
# Compare reactions producing aldehydes in all mechanisms during toluene degradation
# Version 0; Jane Coates 20/11/2014

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

#my @runs = qw( MOZART_tagging );
#my @mechanisms = qw( MOZART-4 );
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2",  "(h) CBM-IV", "(i) CB05" );
my $index = 0;
my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $aldehyde_file = "$base/$run/aldehydes.txt";
    my $aldehydes = get_species($aldehyde_file);
    $families{$mechanisms[$index]} = [ @$aldehydes ];
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(grid) `,
        q` library(scales) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
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
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "Reaction", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c(  "CH3O" = "#b569b3",
                            "Others" = "#623812",
                            "CH3O2 + NO" = "#6c254f", 
                            "C2H5O2 + NO" = "#0c3f78", 
                            "TLBIPERO" = "#ef6638",
                            "NO + RN10O2" = "#2b9eb3",
                            "KETP + NO" = "#f9c500",
                            "NO + TOLP" = "#cc6329",
                            "HC3P + NO" = "#0e5c28",
                            "HOCH2CH2O2 + NO" = "#f7c56c",
                            "NO + TCO3" = "#86b650",
                            "NO + PER1" = "#c9a415",
                            "EPX + O3" = "#8c1531",
                            "NO2 + XOH" = "#77aecc" ) `,
                            q` data$Reaction = factor(data$Reaction, levels = c("TLBIPERO", "CH3O", "CH3O2 + NO", "C2H5O2 + NO", "NO + TOLP", "NO + RN10O2", "KETP + NO", "HC3P + NO", "HOCH2CH2O2 + NO", "NO + TCO3", "NO + PER1", "EPX + O3", "NO2 + XOH", "Others")) `,
);
$R->run(q` scientific_10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab(expression(bold(paste("Reaction Rate (molecules(aldehyde) ", cm^-3, s^-1, ")")))) `,
        q` plot = plot + scale_y_continuous(labels = scientific_10) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 180)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 160, angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140)) `,
        q` plot = plot + theme(axis.ticks.length = unit(2, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 150)) `,
        q` plot = plot + theme(legend.key.size = unit(8, "cm")) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(0.96, 0.10)) `,
        q` plot = plot + theme(legend.justification = c(0.96, 0.10)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Reaction))) `,
);

$R->run(q` CairoPDF(file = "toluene_aldehyde_production_budget.pdf", width = 141, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;

    my ($producers, $producer_yields, $consumers, $consumer_yields, %production, %consumption);
    if (exists $families{$mechanism}) {
        $kpp->family({
                name    => $mechanism,
                members => $families{$mechanism},
                weights => $weights{$mechanism},
        });
        $producers = $kpp->producing($mechanism);
        $producer_yields = $kpp->effect_on($mechanism, $producers);
        $consumers = $kpp->consuming($mechanism);
        $consumer_yields = $kpp->effect_on($mechanism, $consumers);
    } else {
        print "No family for $mechanism\n";
    }
    print "no producers found for $mechanism\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        next unless ($reaction =~ /TOL/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        if ($reactants =~ /MO2 \+ NO|MEO2 \+ NO/) {
            $reactants = "CH3O2 + NO";
        } elsif ($reactants =~ /ETHP \+ NO/) {
            $reactants = "C2H5O2 + NO";
        } elsif ($reactants =~ /RA16O2|TOLO2/) {
            $reactants = "NO + TOLP";
        }
        $production{$reactants} += $rate(1:$NTIME-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        next unless ($reaction =~ /TOL/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        if ($reactants =~ /MO2 \+ NO|MEO2 \+ NO/) {
            $reactants = "CH3O2 + NO";
        } elsif ($reactants =~ /ETHP \+ NO/) {
            $reactants = "C2H5O2 + NO";
        } elsif ($reactants =~ /RA16O2|TOLO2/) {
            $reactants = "NO + TOLP";
        }
        $consumption{$reactants} += $rate(1:$NTIME-2);
    }
    remove_common_processes(\%production, \%consumption);

    my $others = 1e7;
    foreach my $reaction (keys %production) {
        if ($production{$reaction}->sum < $others) {
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
    return \@final_sorted_data;
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

sub get_species {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file for reading : $!";
    chomp(my @all = <$in>);
    close $in;
    return \@all;
} 
