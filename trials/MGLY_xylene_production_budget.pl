#! /usr/bin/env perl
# Plot MGLY production and consumption budgets for all mechanisms from xylene analysis
# Version 0: Jane Coates 6/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0); 
my $n_per_day = 43200 / $dt;
my $n_days = int ($NTIME / $n_per_day);

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
my @base_name = qw( MGLYOX MGLYOX CARB6 CH3COCHO MGLY MGLY MGLY MGLY MGLY);
my @parent = qw( MXYL MXYL MXYL TOLUENE XYL XYL XYM MXYL MXYL);
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
#my @base_name = qw( MGLY );
my $array_index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $spcfile = "$base/$run/gas.spc";
    my $all_tagged_species = get_tagged_species($base_name[$array_index], $spcfile); 
    $families{$mechanisms[$array_index]} = [ @$all_tagged_species ];
    ($plot_data{$mechanisms[$array_index]}, $legend{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index], $parent[$array_index]);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(dplyr) `,
        q` library(gridExtra) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "Alkoxy Radicals" = "#0e5c28",
                            "KETP + NO" = "#6c254f",
                            "CO2H3CO3 + NO" = "#2b9eb3",
                            "MXYOBPEROH + OH" = "#e7e85e",
                            "C5COO2NO2 + OH" = "#86b650",
                            "NO + XYLP" = "#ae4901", "NO + RA19BO2" = "#ae4901",
                            "OH + RA19NO3" = "#f9c500",
                            "OH + RA19OOH" = "#8ed6d2",
                            "NO + TCO3" = "#0352cb",
                            "NO + PER2" = "#1b695b",
                            "NO + TOLO2" = "#77aecc",
                            "BIGALD + hv" = "#8c1531", 
                            "OH + XYL" = "#0c3f78") `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 2.5e7), breaks = seq(0, 2.5e7, 1e7)) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(plot.title = element_text(size = 140, face = "bold")) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 80, angle = 45, vjust = 0.5)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 80)) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.text = element_text(size = 67)) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(legend.key.size = unit(6.5, "cm")) ;
                                                            plot = plot + scale_fill_manual(limits = legend, values = my.colours) ;
                                                            return(plot) } `,
);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('Time', [@days]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_} $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('legend', [@{$legend{$run}}]);
    $R->set('mechanism', $run);
    if ($run =~ /MCM/) {
        $R->run(q` data = mutate(data, "Alkoxy Radicals" = MXYBIPERO + C3MDIALO + MXYOLO + MXYMUCO) `,
                q` data = select(data, -c(MXYBIPERO, C3MDIALO, MXYOLO, MXYMUCO)) `,
                q` data = data[, c("Time", "Alkoxy Radicals", "CO2H3CO3 + NO", "MXYOBPEROH + OH", "C5COO2NO2 + OH", "Production Others")] `,
                q` legend = rev(c("Alkoxy Radicals", "CO2H3CO3 + NO", "MXYOBPEROH + OH", "C5COO2NO2 + OH", "Production Others")) `,
        );
        #my $p = $R->run(q` print(data) `);
        #print "$p\n";
    }
    $R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Reaction", value.name = "Rate") `,
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `, #add plot to list 
    );
}

$R->run(q` CairoPDF(file = "MGLY_xylene_production_comparison.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[2]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[6]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[7]], 
                                                    plots[[8]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) ,
                                                    plots[[9]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) ,
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 160), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $species, $VOC) = @_;

    my ($consumers, $producers, $consumer_yields, $producer_yields, %production_reaction_rates, %consumption_reaction_rates);
    if (exists $families{$species}) { 
        $kpp->family({ 
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $consumers = $kpp->consuming($species);
        $producers = $kpp->producing($species);
        $consumer_yields = $kpp->effect_on($species, $consumers); 
        $producer_yields = $kpp->effect_on($species, $producers);  
    } else {
        print "No family found for $species\n";
    }

    #check that species reactions are found
    die "No producers found for $species\n" if (@$producers == 0);
    die "No consumers found for $species\n" if (@$consumers == 0);

    my $prod_others_max = 1e6;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($number, $parent) = split /_/, $reaction;
        next unless(defined $parent and $parent eq $VOC);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $production_reaction_rates{$reactants} += $rate(1:$NTIME-2);
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my ($number, $parent) = split /_/, $reaction;
        next unless(defined $parent and $parent eq $VOC);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $consumption_reaction_rates{$reactants} += $rate(1:$NTIME-2);
    } 

    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    foreach my $reaction (keys %production_reaction_rates) {
        if ($production_reaction_rates{$reaction}->sum < $prod_others_max) {
            $production_reaction_rates{"Production Others"} += $production_reaction_rates{$reaction};
            delete $production_reaction_rates{$reaction};
        }
    }

    foreach my $reaction (keys %production_reaction_rates) {
        my $reshape = $production_reaction_rates{$reaction}->copy->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_reaction_rates{$reaction} = $integrate;
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    
    my @final_sorted_data;
    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $production_reaction_rates{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $production_reaction_rates{'Production Others'} } if (defined $production_reaction_rates{'Production Others'}); 

    my (@legend_pos, @legend_neg, @legend);
    foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            if ($ref->{$item}->sum > 0) {
                push @legend_pos, $item;
            } else {
                push @legend_neg, $item;
            }
        }
    } 
    push @legend, reverse @legend_pos;
    push @legend, @legend_neg;
    return (\@final_sorted_data, \@legend);
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
