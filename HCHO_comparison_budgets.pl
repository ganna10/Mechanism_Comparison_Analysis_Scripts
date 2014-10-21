#! /usr/bin/env perl
# Plot HCHO production and comparison budgets for all mechanisms
# Version 0: Jane Coates 17/9/2014
# Version 1: Jane Coates 29/9/2014 changing to day-time production only plots

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog; 
my @time_blocks;
foreach my $time (@time_axis) {#map to day and night
    if ($time <= 12) {
        push @time_blocks, "Day 1";
    } elsif ($time > 12 and $time <= 24) {
        push @time_blocks, "Night 1";
    } elsif ($time > 24 and $time <= 36) {
        push @time_blocks, "Day 2";
    } elsif ($time > 36 and $time <= 48) {
        push @time_blocks, "Night 2";
    } elsif ($time > 48 and $time <= 60) {
        push @time_blocks, "Day 3",
    } elsif ($time > 60 and $time <= 72) {
        push @time_blocks, "Night 3";
    } elsif ($time > 72 and $time <= 84) {
        push @time_blocks, "Day 4";
    } elsif ($time > 84 and $time <= 96) {
        push @time_blocks, "Night 4";
    } elsif ($time > 96 and $time <= 108) {
        push @time_blocks, "Day 5";
    } elsif ($time > 108 and $time <= 120) {
        push @time_blocks, "Night 5";
    } elsif ($time > 120 and $time <= 132) {
        push @time_blocks, "Day 6";
    } elsif ($time > 132 and $time <= 144) {
        push @time_blocks, "Night 6";
    } elsif ($time > 144 and $time <= 156) {
        push @time_blocks, "Day 7";
    } else {
        push @time_blocks, "Night 7";
    }
}

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2",  "(h) CBM-IV", "(i) CB05" );
my @base_name = qw( HCHO HCHO HCHO CH2O HCHO HCHO HCHO HCHO FORM );
#my @runs = qw( RACM_tagging ) ;
#my @mechanisms = qw( RACM );
#my @base_name = qw( HCHO );
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
    ($plot_data{$mechanisms[$array_index]}, $legend{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index]);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(RColorBrewer) `,
        q` library(grid) `,
        q` library(gridExtra) `,
);

$R->run(q` arrange.data = function (data) { data = ddply(data, .(Time), colwise(sum)) ;
                                            data = data[1:7,] ;
                                            data = melt(data, id.vars = c("Time"), variable.name = "Reaction", value.name = "Rate") ;
                                            return(data) } `,
);

            
$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "Consumption Others" = "#1c3e3d",
                            "NO + OLTP = .012 ACD + .06 ACT\n+ .44 ALD + .78 HCHO + .78 HO2\n+ .13 MEK + .97 NO2 + .03 ONIT" = "#603912",
                            "ISOPO2 + NO = .55 CH2O + HO2\n+ .37 HYDRALD + .23 MACR + .32 MVK\n+ .92 NO2 + .08 ONIT" = "#6db875",
                            "ALKO2 + NO = .1 CH2O + .4 CH3CHO\n+ .25 CH3COCH3 + .9 HO2 + .75 MEK\n+ .9 NO2 + .1 ONIT" = "#ae4901",
                            "NO + PO2 = CH2O + CH3CHO\n+ HO2 + NO2" = "#8fd5d3",
                            "CH3COCH3O2 + NO = CH2O\n+ CH3CO3 + NO2" = "#bb8a01",
                            "OH + OLE = ALD2 + HCHO + HO2\n+ PAROP + XO2" = "#4b9483",
                            "ISOP + OH = 0.629 FORM + 0.912 HO2\n+ 0.912 ISPD + 0.991 XO2 + 0.088 XO2N" = "#f8c56c",
                            "OH + OLE = 0.33 ALD2 + 0.62 ALDX\n+ 0.8 FORM + 0.95 HO2 + 0.7 PAROP + 0.8 XO2" = "#5b671f",
                            "ISOPBO = HCHO + HO2 + MVK" = "#9bb08f",
                            "HOCH2OO = CH2O + HO2" = "#0352cb", "HCO3 = FORM + HO2" = "#0352cb",
                            "CH3O = HCHO + HO2" = "#c9a415",
                            "HOCH2CH2O = 2 HCHO + HO2" = "#6c254f", "EO = 2 CH2O + HO2" = "#6c254f",
                            "HOCH2CO3 + NO = HCHO + HO2\n+ NO2" = "#8ed6d2",
                            "CH3COCH2O = CH3CO3 + HCHO" = "#f3aa7f",
                            "NO + RI12O2 = CH3COCH3 + HCHO\n+ HO2 + NO2" = "#0e5c28",
                            "IBUTOLBO = CH3COCH3 + HCHO\n+ HO2" = "#f9c500", 
                            "HYPROPO = CH3CHO + HCHO\n+ HO2" = "#a67c52",
                            "NO + RN9O2 = CH3CHO + HCHO\n+ HO2 + NO2" = "#86b650", "NO + PO2 = CH2O + CH3CHO + HO2 + NO\n2" = "#86b650", "NO + OLTP = ALD + HCHO + HO2 + NO2" = "#86b650", "NO + OLTP = 0.94 ALD + HCHO\n+ HO2 + 0.06 KET + NO2" = "#86b650", "NO + OLTP = .012 ACD + 0.06 ACT + 0\n.44 ALD + 0.78 HCHO + 0.78 HO2 + 0.\n13 MEK + 0.97 NO2 + 0.03 ONIT" = "#86b650",
                            "HOCH2CH2O2 + NO = 2 HCHO\n+ HO2 + NO2" = "#77aecc", 
                            "CH3O2 + NO = HCHO + HO2\n+ NO2" = "#8c1531", "MEO2 + NO = FORM + HO2 + NO2" = "#8c1531", "CH3O2 + NO = CH2O + HO2\n+ NO2" = "#8c1531", "MO2 + NO = HCHO + HO2 + NO2" = "#8c1531", "C2O3 + NO = HCHO + HO2\n+ NO2 + XO2" = "#8c1531",
                            "CH4 + OH = HCHO + HO2 + XO2" = "#898989",
                            "NO + RU14O2 = HCHO + HO2\n+ NO2 + UCARB10" = "#dc3522", 
                            "C2H4 + OH = 0.5 CH2O +\n0.75 EO + 0.25 HO2" = "#9bb18d", "ETH + OH = 0.22 ALD2 + 1.56 HCHO\n+ HO2 + XO2" = "#9bb18d", "ETH + OH = 0.22 ALDX +\n1.56 FORM + HO2 + XO2" = "#9bb18d",
                            "NO + RN8O2 = CH3CO3 + HCHO + NO2" = "#623812", 
                            "ISOP + NO = 0.606 HCHO + 0.847 HO2\n+ 0.446 MACR + 0.847 NO2\n+ 0.384 OLT + 0.153 ONIT" = "#58691b",
                            "HC3P + NO = 0.75 ALD + 0.09 HCHO +\n0.964 HO2 + 0.25 KET + 0.964 NO2 + 0.036 ONIT" = "#4c9383", 
                            "ACTP + NO = ACO3 + HCHO + NO2" = "#cc6638",
                            "NO + OL2P = 0.2 ALD + 1.6 HCHO\n+ HO2 + NO2" = "#ef6638", "ETEP + NO = 0.2 ALD +\n1.6 HCHO + HO2 + NO2" = "#ef6638",
                            "OH + OLE =  0.330 ALD2 + 0.620 ALDX\n + 0.800 FORM + 0.950 HO2 + 0.700 PA\nROP + 0.800 XO2" = "#6d6537", "OH + OLE = ALD2 + HCHO + HO2 + PARO\nP + XO2" = "#6d6537",
                            "ISOP + OH = 0.629 FORM + 0.912 HO2 \n+ 0.912 ISPD + 0.991 XO2 + 0.088 XO\n2N" = "#e7e85e", "ISOP + OH = 0.2 ALD2 + 0.2 C2O3\n+ ETH + HCHO + 0.67 HO2\n+ 0.4 MGLY + XO2 + 0.13 XO2N" = "#e7e85e") `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 4e8), breaks = seq(0, 4e8, 1e8)) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
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

$R->set('Time', [@time_blocks]);
$R->run(q` plots = list() `); #list to fill up with plots from each mechanism
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            next if ($reaction =~ /^CH3CO3 = \.3 CH2O|^CH3O2 = \.7 HCHO/);
            $R->set('reaction', $reaction);
            $R->set('rate', [@{$ref->{$reaction}}]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('legend', [@{$legend{$run}}]);
    $R->set('mechanism', $run);
    $R->run(q` data = arrange.data(data) `,
            q` reaction.levels = levels(factor(data$Reaction)) `,
            q` data$Reaction = ordered(data$Reaction, levels = reaction.levels) `, 
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `, #add plot to list 
    );
}

$R->run(q` CairoPDF(file = "HCHO_budget_comparison.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[2]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[3]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[4]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                                    plots[[5]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[6]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[7]], 
                                                    plots[[8]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    plots[[9]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 3), 
                                       nrow = 1, ncol = 1,
                                       left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm ^-3, s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop;

sub get_data {
    my ($kpp, $mecca, $species) = @_;

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

    my $prod_others_max = 8e6;
    my $max_string_width = 28;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        $reaction_string =~ s/\n$//;
        next if ($reaction_string =~ /^CH3CO3 = \.3 CH2O|^CH3O2 = \.7 CH2O/);
        $reaction_string = string_mapping($reaction_string);
        if ($rate->sum <= $prod_others_max) {
            $production_reaction_rates{"Production Others"} += $rate;
        } else {
            $production_reaction_rates{$reaction_string} += $rate;
        }
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        $reaction_string =~ s/(.{1,$max_string_width})/$1\n/gs;
        $reaction_string =~ s/\n$//;
        $reaction_string = string_mapping($reaction_string);
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 

    sub string_mapping {
        my ($string) = @_;

        if ($string eq "HYPROPO = CH3CHO + HCHO + HO\n2") {
            $string = "HYPROPO = CH3CHO + HCHO\n+ HO2";
        } elsif ($string eq "NO + RN9O2 = CH3CHO + HCHO +\n HO2 + NO2") {
            $string = "NO + RN9O2 = CH3CHO + HCHO\n+ HO2 + NO2";
        } elsif ($string eq "HOCH2CH2O2 + NO = 2 HCHO + H\nO2 + NO2") {
            $string = "HOCH2CH2O2 + NO = 2 HCHO\n+ HO2 + NO2";
        } elsif ($string eq "CH3O2 + NO = HCHO + HO2 + NO\n2") {
            $string = "CH3O2 + NO = HCHO + HO2\n+ NO2";
        } elsif ($string =~ /HC3P \+ NO = 0\.75 ALD/) {
            $string = "HC3P + NO = 0.75 ALD + 0.09 HCHO +\n0.964 HO2 + 0.25 KET + 0.964 NO2 + 0.036 ONIT";
        } elsif ($string =~ /NO \+ OLTP = ALD \+ HCHO \+ HO2/) {
            $string = "NO + OLTP = ALD + HCHO + HO2 + NO2";
        } elsif ($string =~ /NO \+ OL2P/) {
            $string = "NO + OL2P = 0.2 ALD + 1.6 HCHO\n+ HO2 + NO2";
        } elsif ($string =~ /NO \+ OLTP = 0\.94 ALD/) {
            $string = "NO + OLTP = 0.94 ALD + HCHO\n+ HO2 + 0.06 KET + NO2";
        } elsif ($string =~ /ETEP \+ NO = 0\.2 ALD/) {
            $string = "ETEP + NO = 0.2 ALD +\n1.6 HCHO + HO2 + NO2";
        } elsif ($string =~ /C2H4 \+ OH/) {
            $string = "C2H4 + OH = 0.5 CH2O +\n0.75 EO + 0.25 HO2";
        } elsif ($string =~ /CH3O2 \+ NO = CH2O/) {
            $string = "CH3O2 + NO = CH2O + HO2\n+ NO2";
        } elsif ($string =~ /^ISOP \+ OH = 0\.2 ALD2/) {
            $string = "ISOP + OH = 0.2 ALD2 + 0.2 C2O3\n+ ETH + HCHO + 0.67 HO2\n+ 0.4 MGLY + XO2 + 0.13 XO2N";
        } elsif ($string =~ /^C2O3 \+ NO/) {
            $string = "C2O3 + NO = HCHO + HO2\n+ NO2 + XO2";
        } elsif ($string =~ /ETH \+ OH = 0\.22 ALD2/) {
            $string = "ETH + OH = 0.22 ALD2 + 1.56 HCHO\n+ HO2 + XO2";
        } elsif ($string =~ /ETH \+ OH = 0\.220 ALDX/) {
            $string = "ETH + OH = 0.22 ALDX +\n1.56 FORM + HO2 + XO2";
        } elsif ($string =~ /IBUTOLBO/) {
            $string = "IBUTOLBO = CH3COCH3 + HCHO\n+ HO2";
        } elsif ($string =~ /HOCH2CO3 \+ NO/) {
            $string = "HOCH2CO3 + NO = HCHO + HO2\n+ NO2";
        } elsif ($string =~ /NO \+ RI12O2/) {
            $string = "NO + RI12O2 = CH3COCH3 + HCHO\n+ HO2 + NO2";
        } elsif ($string =~ /NO \+ RU14O2/) {
            $string = "NO + RU14O2 = HCHO + HO2\n+ NO2 + UCARB10";
        } elsif ($string =~ /NO \+ RN8O2/) {
            $string = "NO + RN8O2 = CH3CO3 + HCHO + NO2";
        } elsif ($string =~ /ISOP \+ NO = /) {
            $string = "ISOP + NO = 0.606 HCHO + 0.847 HO2\n+ 0.446 MACR + 0.847 NO2\n+ 0.384 OLT + 0.153 ONIT";
        } elsif ($string =~ /ACTP \+ NO/) {
            $string = "ACTP + NO = ACO3 + HCHO + NO2";
        } elsif ($string =~ /NO \+ OLTP/) {
            $string = "NO + OLTP = .012 ACD + .06 ACT\n+ .44 ALD + .78 HCHO + .78 HO2\n+ .13 MEK + .97 NO2 + .03 ONIT";
        } elsif ($string =~ /ISOPO2 \+ NO/) {
            $string = "ISOPO2 + NO = .55 CH2O + HO2\n+ .37 HYDRALD + .23 MACR + .32 MVK\n+ .92 NO2 + .08 ONIT";
        } elsif ($string =~ /ALKO2 \+ NO/) {
            $string = "ALKO2 + NO = .1 CH2O + .4 CH3CHO\n+ .25 CH3COCH3 + .9 HO2 + .75 MEK\n+ .9 NO2 + .1 ONIT";
        } elsif ($string =~ /NO \+ PO2/) {
            $string = "NO + PO2 = CH2O + CH3CHO\n+ HO2 + NO2";
        } elsif ($string =~ /CH3COCH2O2 \+ NO/) {
            $string = "CH3COCH3O2 + NO = CH2O\n+ CH3CO3 + NO2";
        } elsif ($string =~ /OH \+ OLE = ALD2/) {
            $string = "OH + OLE = ALD2 + HCHO + HO2\n+ PAROP + XO2";
        } elsif ($string =~ /ISOP \+ OH/) {
            $string = "ISOP + OH = 0.629 FORM + 0.912 HO2\n+ 0.912 ISPD + 0.991 XO2 + 0.088 XO2N";
        } elsif ($string =~ /OH \+ OLE = 0\.330 ALD2/) {
            $string = "OH + OLE = 0.33 ALD2 + 0.62 ALDX\n+ 0.8 FORM + 0.95 HO2 + 0.7 PAROP + 0.8 XO2";
        }
        return $string;
    }
    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    
    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    
    my @final_sorted_data;
    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $production_reaction_rates{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $production_reaction_rates{'Production Others'} } if (defined $production_reaction_rates{'Production Others'}); 

    my (@plot_data, @legend_pos, @legend_neg, @legend);
    foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
        foreach my $item (keys %$ref) {
            if ($ref->{$item}->sum > 0) {
                push @legend_pos, $item;
            } else {
                push @legend_neg, $item;
            }
            my @rate_array = map { $_ } $ref->{$item}->dog;
            my @new_array = splice @rate_array, 1, 504;
            push @plot_data, { $item => \@new_array };
        }
    } 
    push @legend, reverse @legend_pos;
    push @legend, @legend_neg;
    return (\@plot_data, \@legend);
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
