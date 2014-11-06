#! /usr/bin/env perl
# Allocation MGLY net production to VOC sources in each mechanism
# Version 0: Jane Coates 5/11/2014

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
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT ;
my $N_DAYS = int ( $NTIME / $N_PER_DAY );

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCM v3.2", "(b) MCM v3.1", "(c) CRI v2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
my @base_name = qw( MGLYOX MGLYOX CARB6 CH3COCHO MGLY MGLY MGLY MGLY MGLY );
#my @runs = qw( CBM4_tagging CB05_tagging);
#my @mechanisms = qw( CBM-IV CB05);
#my @base_name = qw( MGLY MGLY);
my $index = 0;

my (%families, %weights, %plot_data);
foreach my $run (@runs) {
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $spc_file = "$base/$run/gas.spc";
    my $all_tagged_species = get_tagged_species($base_name[$index], $spc_file); 
    $families{$mechanisms[$index]} = [ @$all_tagged_species ];
    $plot_data{$mechanisms[$index]} = get_data($kpp, $mecca, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(reshape2) `,
        q` library(grid) `,
        q` library(scales) `,
);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('Time', [@days]);
$R->run(q` data = data.frame() `);

foreach my $mechanism (keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$mechanism}}) {
        foreach my $VOC (keys %$ref) {
            $R->set('voc', $VOC);
            $R->set('rate', [map { $_ } $ref->{$VOC}->dog]) ;
            $R->run(q` pre[voc] = rate `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Mechanism"), variable.name = "VOC", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(levels(data$VOC)) `);
#print $p, "\n";

$R->run(q` scientific_10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `);
$R->run(q` my.colours = c( "Methane" = "#c9a415", "Ethane" = "#ae4901", "Propane" = "#f9c600", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb", "Ethene" = "#86b650", "Propene" = "#6c254f", "Butene" = "#ee6738", "2-Methylpropene" = "#58691b", "Isoprene" = "#8ed6d5", "Benzene" = "#f3aa7f", "Toluene" = "#c65d6c", "m-Xylene" = "#888a87", "o-Xylene" = "#0e5c28", "p-Xylene" = "#b569b3", "Ethylbenzene" = "#2c9def", "Others" = "#696537" ) `,
    q` data$VOC = factor(data$VOC, levels = c("Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Ethylbenzene", "Others")) `,
);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = VOC)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
        q` plot = plot + xlab("\n") `,
        q` plot = plot + ylab(expression(bold(paste("Methyl Glyoxal Production Rates Allocated to Parent VOC (molecules (MGLY) ", cm^-3, s^-1,")")))) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_y_continuous(label = scientific_10) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 100)) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 180, face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 160)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 150)) `,
        q` plot = plot + theme(axis.ticks.length = unit(2.5, "cm")) `,
        q` plot = plot + theme(axis.ticks.margin = unit(1, "cm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(10, "cm")) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + scale_fill_manual(values = my.colours, guide = guide_legend(label.theme = element_text(size = 140, angle = 0), label.vjust = 0.5)) `,
);

$R->run(q` CairoPDF(file = "MGLY_VOC_allocated_production_rates.pdf", width = 165, height = 200) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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

    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my ($number, $parent) = split /_/, $reaction;
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        if (defined $parent) {
            my $chemical_name = get_name($parent);
            $production_reaction_rates{$chemical_name} += $rate(1:$NTIME-2);
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
            $production_reaction_rates{$reaction_string} += $rate(1:$NTIME-2);
        }
    }

    for (0..$#$consumers) { #get rates for all consuming reactions
        my $reaction = $consumers->[$_];
        my ($number, $parent) = split /_/, $reaction;
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        if (defined $parent) {
            my $chemical_name = get_name($parent);
            $consumption_reaction_rates{$chemical_name} += $rate(1:$NTIME-2);
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
            $consumption_reaction_rates{$reaction_string} += $rate(1:$NTIME-2);
        }
    } 

    #remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    my $others_max = 2e6;
    foreach my $VOC (keys %production_reaction_rates) {
        if ($production_reaction_rates{$VOC}->sum < $others_max) {
            $production_reaction_rates{"Others"} += $production_reaction_rates{$VOC};
            delete $production_reaction_rates{$VOC};
        }
    }

    foreach my $VOC (keys %production_reaction_rates) {
        my $reshape = $production_reaction_rates{$VOC}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_reaction_rates{$VOC} = $integrate;
    }
    
    my $sort_function = sub { $_[0]->sum };
    my @prod_sorted_data = sort { &$sort_function($production_reaction_rates{$b}) <=> &$sort_function($production_reaction_rates{$a}) } keys %production_reaction_rates;
    
    my @final_sorted_data;
    foreach (@prod_sorted_data) { 
        next if ($_ eq 'Others') ;
        push @final_sorted_data, { $_ => $production_reaction_rates{$_} };
    } 
    push @final_sorted_data, { 'Others' => $production_reaction_rates{'Others'} } if (defined $production_reaction_rates{'Others'}); 

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
    } elsif ($parent =~ /MXYL|XYL|XYM/) {
        $parent = "m-Xylene";
    } elsif ($parent =~ /OXYL|XYO/) {
        $parent = "o-Xylene";
    } elsif ($parent =~ /PXYL|XYP/) {
        $parent = "p-Xylene";
    } elsif ($parent eq "EBENZ") {
        $parent = "Ethylbenzene";
    } else {
        print "No chemical name for $parent\n";
    }
    return $parent;
}
