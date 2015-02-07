#! /usr/bin/env perl
# compare radical production as prescribed by NO source calculation in CBs and MCM v3.2, budget expressed as which radicals are produced
# Version 0: Jane Coates 2/2/2015

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
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY; 

#my @mechanisms = ( "MCMv3.2", "CBM-IV", "CB05" );
my @mechanisms = ( "CB05" );
my (%families, %weights, %data);
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $radical_file = "$base/${mechanism}_tagged/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{$mechanism} = [ @radicals ];
    ($data{$mechanism}) = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "MGLYOX + hv" = "#f9c500", "MGLY + hv" = "#f9c500", 
                            "OH + PAR" = "#2c9daf", 
                            "C2O3 + NO" = "#e7e85e",
                            "C2H4 + OH" = "#6db875", 
                            "CXO3 + NO" = "#bb8a01",
                            "O1D" = "#6c254f",
                            "HCHO + hv" = "#0e5c28", "FORM + hv" = "#0e5c28",
                            "CH4 + OH" = "#0d3e76",
                            "C2H6 + OH" = "#ef6638",
                            "MEK + hv" = "#86b650" ) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog ]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Reaction, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "CBM-IV", "CB05")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        #q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(reaction.levels)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + ylab("Reaction Rate (molecules cm-3 s-1)") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 0.7, hjust = 0.8)) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
);

$R->run(q` CairoPDF(file = "radical_NOx_production_budgets_CBs_MCM.pdf", width = 8.5, height = 6 ) `,
        q` print(plot) `, 
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    
    $families{'NOx'} = [ qw( NO NO2 NO3 N2O5 ) ];
    $weights{'NOx'} = { N2O5 => 2 };
    $kpp->family({ #NOx family
                    name    => 'NOx',
                    members =>$families{'NOx'},
                    weights => $weights{'NOx'},
    });

    my ($radical_producers, $radical_producer_yields, $NOx_yields, %production_rates);
    if (exists $families{$mechanism}) { #radicals family
        $kpp->family({                                                                                                                           
                        name    => $mechanism,
                        members => $families{$mechanism},
                        weights => $weights{$mechanism},
        }); 
        $radical_producers = $kpp->producing($mechanism);
        $radical_producer_yields = $kpp->effect_on($mechanism, $radical_producers);
        $NOx_yields = $kpp->effect_on('NOx', $radical_producers);
    } else {
        print "No radical family found for $mechanism\n";
    }
    
    die "No producers found for $mechanism\n" if (@$radical_producers == 0);

    for (0..$#$radical_producers) { #get rates for all radical producing reactions
        next if ($NOx_yields->[$_] < 0);
        my $reaction = $radical_producers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $net_radical_yield = $radical_producer_yields->[$_] - $NOx_yields->[$_];
        next if ($net_radical_yield == 0);
        my $rate = $net_radical_yield * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $reaction_string =~ s/_(.*?)\b//g;
        print "$reaction_string : $net_radical_yield\n";
        #$production_rates{$label} += $rate(1:$NTIME-2);
    }

    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_rates{$reaction} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_plot_data;
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    }
    return \@sorted_plot_data;
}

sub get_species {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die $!;
    local $/ = undef;
    my $lines = <$in>;
    close $in;
    my @species = split /\s+/, $lines;
}
