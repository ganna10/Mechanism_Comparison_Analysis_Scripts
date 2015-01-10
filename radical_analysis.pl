#! /usr/bin/env perl
# compare radical production as prescribed by NO source calculation
# Version 0: Jane Coates 3/10/2014
# Version 1: Jane Coates 30/10/2014 Merging OH+PAR and ROR, as they are tied together, plotting aesthetic improvements
# Version 2: Jane Coates 8/12/2014 updating script for constant emissions runs
# Version 3: Jane Coates 28/12/2014 going back to separate plots for each mechanism and not facetting

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

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05" );
my (%families, %weights, %data, %legend);
foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $radical_file = "$base/${mechanism}_tagged/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{$mechanism} = [ @radicals ];
    ($data{$mechanism}, $legend{$mechanism}) = get_data($mecca, $kpp, $mechanism);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "MGLYOX + hv" = "#f9c500", "CH3COCHO + hv" = "#f9c500", "CARB6 + hv" = "#f9c500", "MGLY + hv" = "#f9c500", 
                            "NO + TCO3" = "#603912",
                            "HC5 + OH" = "#c9a415",
                            "EPX + O3" = "#6db875",
                            "DCB1 + O3" = "#ae4901",
                            "BIGALD + hv" = "#8fd5d3",
                            "OH + PAR" = "#bb8a01", 
                            "C2O3 + NO" = "#e7e85e",
                            "C2H4 + OH" = "#9bb08f", 
                            "CXO3 + NO" = "#2c9daf",
                            "DCB + hv" = "#8ed6d2",
                            "O1D" = "#6c254f",
                            "HCHO + hv" = "#0e5c28", "CH2O + hv" = "#0e5c28", "FORM + hv" = "#0e5c28",
                            "DCB + hv" = "#76afca",
                            "KET + hv" = "#1c3e3d",
                            "CH4 + OH" = "#0d3e76",
                            "C2H6 + OH" = "#cc6329",
                            "MEK + hv" = "#86b650" ) `,
        q` reaction.levels = c( "O1D",
                                "HCHO + hv",
                                "MGLY + hv",
                                "NO + TCO3",
                                "DCB + hv",
                                "KET + hv",
                                "CH4 + OH",
                                "MEK + hv",
                                "OH + PAR",
                                "C2O3 + NO",
                                "C2H4 + OH",
                                "CXO3 + NO",
                                "EPX + O3",
                                "DCB + hv",
                                "C2H6 + OH",
                                "BIGALD + hv",
                                "Production Others" ) `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + scale_x_discrete(expand = c(0, 0)) ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 7e8), breaks = seq(0, 7e8, 1e8), expand = c(0, 0)) ;
                                                            plot = plot + scale_fill_manual(values = my.colours, limits = legend) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + theme(plot.title = element_text(size = 22, face = "bold")) ;
                                                            plot = plot + theme(axis.title = element_blank()) ;
                                                            plot = plot + theme(axis.text.x = element_text(face = "bold", size = 20, angle = 45, vjust = 0.7, hjust = 0.8)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 18)) ;
                                                            plot = plot + theme(legend.position = c(1, 1)) ;
                                                            plot = plot + theme(legend.justification = c(1, 1)) ;
                                                            plot = plot + theme(panel.grid = element_blank()) ;
                                                            plot = plot + theme(panel.border = element_rect(colour = "black")) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.text = element_text(size = 14)) ;
                                                            return(plot) } `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` plots = list() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $reaction (sort keys %$ref) {
            next if ($reaction eq "CH3CO3");
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog ]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('legend', [@{$legend{$mechanism}}]);
    $R->set('mechanism', $mechanism);
    $R->run(q` data = gather(data, Reaction, Rate, -Time) `,
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` CairoPDF(file = "radical_NOx_production_budgets.pdf", width = 14.1, height = 20.0) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[9]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                                    plots[[7]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[8]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[6]] + theme(axis.title.y = element_blank()),
                                                    plots[[2]] + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    plots[[1]] + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()),
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob("Reaction Rate (molecules cm-3 s-1)", gp = gpar(fontface = "bold", fontsize = 26), rot = 90, vjust = 0.5) ) `, 
        q` print(multiplot) `, 
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
        my ($reactants, $products) = split / = /, $reaction_string;
        if ($reactants eq "ROR") {
            $reactants = "OH + PAR" ;
        } elsif ($reactants eq "ETH + OH" and $mechanism =~ /CB/) {
            $reactants = "C2H4 + OH";
        } elsif ($reactants eq "ETHA + OH") {
            $reactants = "C2H6 + OH";
        }
        $production_rates{$reactants} += $rate(1:$NTIME-2);
    }

    my $others = 3.3e7;
    foreach my $reaction (keys %production_rates) {
        if ($production_rates{$reaction}->sum < $others) {
            $production_rates{'Production Others'} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
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
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); #add Production Others to the beginning 
    my (@positive, @negative, @legend);
    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            if ($ref->{$item}->sum > 0) {
                push @positive, $item;
            } else {
                push @negative, $item;
            }
        }
    }
    push @legend, reverse @positive;
    push @legend, @negative;
    return (\@sorted_plot_data, \@legend);
}

sub get_species {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die $!;
    local $/ = undef;
    my $lines = <$in>;
    close $in;
    my @species = split /\s+/, $lines;
}
