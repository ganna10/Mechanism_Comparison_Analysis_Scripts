#! /usr/bin/env perl
# CH3CHO budget during propane degradation in MCM and RADM2
# Version 0: Jane Coates 13/2/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $N_TIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $N_TIME / $N_PER_DAY;

my @mechanisms = qw( MCMv3.2 RADM2 );
my @species = qw( CH3CHO ALD );
my $index = 0;
my (%families, %weights, %data);
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $spcfile = "$base/${mechanism}_tagged/gas.spc";
    my $all_tagged_species = get_tagged_species($species[$index], $spcfile); 
    $families{$mechanism} = [ @$all_tagged_species ];
    $data{$mechanism} = get_data($mecca, $kpp, $mechanism);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    print "$mechanism\n";
    my $sum = 0;
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $sum += $ref->{$reaction};
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    print "$sum\n";
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` data$Reaction = factor(data$Reaction, levels = rev(c("C2H5O", "OH + PPN", "HC3P + NO", "HC3 + OH", "ETHP + NO", "OH + OP2", "Others"))) `,
        q` my.colours = c( "C2H5O" = "#2b9eb3",
                            "OH + PPN" = "#b569b3",
                            "HC3P + NO" = "#6c254f",  
                            "HC3 + OH" = "#f9c500",
                            "ETHP + NO" = "#0e5c28",
                            "OH + OP2" = "#ef6638",
                            "Others" = "#696537") `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + ylab("Reaction Rate (molecules cm-3 s-1)") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(1, 1)) `,
        q` plot = plot + theme(legend.justification = c(1, 1)) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "C3H8_CH3CHO_budget.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;

    my ($producers, $producer_yields, $consumers, $consumer_yields, %production_rates, %consumption_rates);
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
        print "No family found for $mechanism\n";
    }
    print "No producers for $mechanism\n" if (@$producers == 0);
    print "No consumers for $mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        next unless ($reaction =~ /_C3H8|_HC3/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $string = $kpp->reaction_string($reaction);
        $string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $string;
        $production_rates{$reactants} += $rate(1:$N_TIME-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        next unless ($reaction =~ /_C3H8|_HC3/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $string = $kpp->reaction_string($reaction);
        $string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $string;
        $consumption_rates{$reactants} += $rate(1:$N_TIME-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

    my $others = 5e5;
    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate *= 0.628 if ($mechanism =~ /RA/);
        if ($integrate->sum < $others) {
            $production_rates{"Others"} += $integrate(0:13:2);
            delete $production_rates{$reaction};
        } else {
            $production_rates{$reaction} = $integrate(0:13:2);
        }
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    my @sorted;
    foreach (@sorted_prod) {
        next if ($_ =~ /Others/);
        push @sorted, { $_ => $production_rates{$_} };
    }
    push @sorted, { "Others" => $production_rates{"Others"} } if (exists $production_rates{"Others"});
    return \@sorted;
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

sub get_tagged_species {
    my ($species, $spc_file) = @_;
    my $all_species = read_file($spc_file);
    my @tagged_species;
    foreach my $line (@$all_species) {
        next unless ($line =~ /^$species/);
        $line =~ s/\s=.*$//;
        push @tagged_species, $line;
    }
    return \@tagged_species;
}

sub read_file {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file for reading : $!";
    chomp(my @all = <$in>);
    close $in;
    return \@all;
} 
