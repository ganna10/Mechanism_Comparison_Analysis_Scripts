#! /usr/bin/env perl
# compare radical production 
# Version 0: Jane Coates 8/1/2015

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

my @mechanisms = ( "MCMv3.2", "RADM2" );
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
        q` library(dplyr) `,
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
#        q` reaction.levels = c( "O1D",
#                                "HCHO + hv",
#                                "MGLY + hv",
#                                "NO + TCO3",
#                                "DCB + hv",
#                                "KET + hv",
#                                "CH4 + OH",
#                                "MEK + hv",
#                                "OH + PAR",
#                                "C2O3 + NO",
#                                "C2H4 + OH",
#                                "CXO3 + NO",
#                                "EPX + O3",
#                                "DCB + hv",
#                                "C2H6 + OH",
#                                "BIGALD + hv",
#                                "Production Others" ) `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + scale_fill_manual(values = my.colours, limits = legend) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + theme(plot.title = element_text(size = 22, face = "bold")) ;
                                                            plot = plot + theme(axis.title = element_blank()) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.7, hjust = 0.8)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 18)) ;
                                                            plot = plot + theme(plot.margin = unit(c(0, 0, 0, -0.04), "cm")) ;
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
    $R->run(q` HONO = filter(data, Reaction == "HONO + hv") `);
my $p = $R->run(q` print(HONO) `);
print $p, "\n";
}

$R->run(q` CairoPDF(file = "radical_production_budgets.pdf", width = 16.9, height = 26.0) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] ,
                                                    plots[[2]] ,
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm^-3, s^-1, ")"))), gp = gpar(fontsize = 26), rot = 90, vjust = 0.5) ) `, 
        q` print(multiplot) `, 
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism) = @_;
    
    my ($producers, $producer_yields, $consumers, $consumer_yields, %production_rates, %consumption_rates);
    if (exists $families{$mechanism}) { #radicals family
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
        print "No radical family found for $mechanism\n";
    }
    
    die "No producers found for $mechanism\n" if (@$producers == 0);
    die "No consumers found for $mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) { #get rates for all radical producing reactions
        my $reaction = $producers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_]; 
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

    for (0..$#$consumers) { #get rates for all radical consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_]; 
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
        $consumption_rates{$reactants} += $rate(1:$NTIME-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

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
