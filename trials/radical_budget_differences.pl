#!/usr/bin/perl
# Plotting differences between radical process contributions from each mechanism relative to MCM v3.2
# Version 0 : Jane Coates 16/9/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

#Create x-axis for plot in hours
my $run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $mecca = MECCA->new($run); 
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

my (%families, %weights, %production_rates, %consumption_rates);
my $base = "/work/users/jco/MECCA";
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
#my @runs = qw( MCM_3.2_tagged RACM_tagging ) ;
#my @mechanisms = qw( MCMv3.2 RACM );
my $array_index = 0;

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $eqnfile = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqnfile); 
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{$mechanisms[$array_index]} = [ @radicals ];
    ($production_rates{$mechanisms[$array_index]}, $consumption_rates{$mechanisms[$array_index]}) = get_data($kpp, $mecca, $mechanisms[$array_index]);
    $array_index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(plyr) `,
        q` library(grid) `,
        q` library(gridExtra) `,
        q` library(RColorBrewer) `,
);

$R->set('time', [@time_blocks]);
$R->run(q` mcm3.2.data = data.frame(time) `,
        q` mcm3.1.data = data.frame(time) `,
        q` cri.data = data.frame(time) `,
        q` mozart.data = data.frame(time) `,
        q` radm2.data = data.frame(time) `,
        q` racm.data = data.frame(time) `,
        q` racm2.data = data.frame(time) `,
        q` cbm4.data = data.frame(time) `,
        q` cb05.data = data.frame(time) `,
);

foreach my $mechanism (sort keys %production_rates) {
    $R->set('mechanism', $mechanism);
    foreach my $reaction (sort keys %{$production_rates{$mechanism}}) {
        my @rates = map { $_ } $production_rates{$mechanism}{$reaction}->dog;
        $R->set('reaction', $reaction);
        $R->set('rate', [@rates]);
        if ($mechanism eq "MCMv3.2") {
            $R->run(q` mcm3.2.data[reaction] = rate `);
        } elsif ($mechanism eq "MCMv3.1") {
            $R->run(q` mcm3.1.data[reaction] = rate `);
        } elsif ($mechanism eq "CRIv2") {
            $R->run(q` cri.data[reaction] = rate `);
        } elsif ($mechanism eq "MOZART-4") {
            $R->run(q` mozart.data[reaction] = rate `);
        } elsif ($mechanism eq "RADM2") {
            $R->run(q` radm2.data[reaction] = rate `);
        } elsif ($mechanism eq "RACM") {
            $R->run(q` racm.data[reaction] = rate `);
        } elsif ($mechanism eq "RACM2") {
            $R->run(q` racm2.data[reaction] = rate `);
        } elsif ($mechanism eq "CBM-IV") {
            $R->run(q` cbm4.data[reaction] = rate `);
        } elsif ($mechanism eq "CB05") {
            $R->run(q` cb05.data[reaction] = rate `);
        }
    }
}
foreach my $mechanism (sort keys %consumption_rates) {
    $R->set('mechanism', $mechanism);
    foreach my $reaction (sort keys %{$consumption_rates{$mechanism}}) {
        my @rates = map { $_ } $consumption_rates{$mechanism}{$reaction}->dog;
        $R->set('reaction', $reaction);
        $R->set('rate', [@rates]);
        if ($mechanism eq "MCMv3.2") {
            $R->run(q` mcm3.2.data[reaction] = rate `);
        } elsif ($mechanism eq "MCMv3.1") {
            $R->run(q` mcm3.1.data[reaction] = rate `);
        } elsif ($mechanism eq "CRIv2") {
            $R->run(q` cri.data[reaction] = rate `);
        } elsif ($mechanism eq "MOZART-4") {
            $R->run(q` mozart.data[reaction] = rate `);
        } elsif ($mechanism eq "RADM2") {
            $R->run(q` radm2.data[reaction] = rate `);
        } elsif ($mechanism eq "RACM") {
            $R->run(q` racm.data[reaction] = rate `);
        } elsif ($mechanism eq "RACM2") {
            $R->run(q` racm2.data[reaction] = rate `);
        } elsif ($mechanism eq "CBM-IV") {
            $R->run(q` cbm4.data[reaction] = rate `);
        } elsif ($mechanism eq "CB05") {
            $R->run(q` cb05.data[reaction] = rate `);
        }
    }
}

$R->run(q` mcm3.2.data = ddply(mcm3.2.data, .(time), colwise(sum)) `,
        q` mcm3.1.data = ddply(mcm3.1.data, .(time), colwise(sum)) `,
        q` cri.data = ddply(cri.data, .(time), colwise(sum)) `,
        q` mozart.data = ddply(mozart.data, .(time), colwise(sum)) `,
        q` radm2.data = ddply(radm2.data, .(time), colwise(sum)) `,
        q` racm.data = ddply(racm.data, .(time), colwise(sum)) `,
        q` racm2.data = ddply(racm2.data, .(time), colwise(sum)) `,
        q` cbm4.data = ddply(cbm4.data, .(time), colwise(sum)) `,
        q` cb05.data = ddply(cb05.data, .(time), colwise(sum)) `,
        q` mcm3.2.data = mcm3.2.data[1:7,] `,
        q` mcm3.1.data = mcm3.1.data[1:7,] `,
        q` cri.data = cri.data[1:7,] `,
        q` mozart.data = mozart.data[1:7,] `,
        q` radm2.data = radm2.data[1:7,] `,
        q` racm.data = racm.data[1:7,] `,
        q` racm2.data = racm2.data[1:7,] `,
        q` cbm4.data = cbm4.data[1:7,] `,
        q` cb05.data = cb05.data[1:7,] `,
);

$R->run(q` days = c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") `,
        q` HCHO.hv.diff = data.frame(days) `,
        q` HCHO.hv.diff["CB05"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - cb05.data["FORM + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["CBM-IV"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - cbm4.data["HCHO + hv = CO + 2.0 HO2"] `,
        q` HCHO.hv.diff["CRIv2"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - cri.data["HCHO + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["MCMv3.1"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - mcm3.1.data["HCHO + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["MOZART-4"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - mozart.data["CH2O + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["RACM"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - racm.data["HCHO + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["RACM2"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - racm2.data["HCHO + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff["RADM2"] = mcm3.2.data["HCHO + hv = CO + 2 HO2"] - radm2.data["HCHO + hv = CO + 2 HO2"] `,
        q` HCHO.hv.diff = melt(HCHO.hv.diff, id.vars = c("days"), variable.name = "Mechanism", value.name = "Rates") `,
);

#my $p = $R->run(q` print(HCHO.hv.diff) `);
#print $p, "\n";

$R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#662c91", "#12b2b2", "#b33893", "#a11d20") `); #"#ed2d2e", 
$R->run(q` plot = ggplot(HCHO.hv.diff, aes(x = days, y = Rates, group = Mechanism, colour = Mechanism)) `,
        q` plot = plot + geom_line(size = 2) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + xlab("\n") `,
        q` plot = plot + ylab(expression(bold(paste("Photolysis Rate Difference from MCM v3.2 (", s^-1, ")")))) `,
        q` plot = plot + ggtitle("HCHO + hv = CO + 2 HO2") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 12, face = "bold")) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 10)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 10)) `,
        q` plot = plot + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` plot = plot + theme(legend.title = element_text(size = 12, face = "bold")) `,
        q` plot = plot + theme(legend.text = element_text(size = 10)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(1, "cm")) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HCHO_hv_diff_from_MCM3_2.pdf", width = 10, height = 7) `,
        q` print(plot) `,
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

    my $prod_others_max = 7e7;
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        $rate = $rate(1:$NTIME-2);
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
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
        $rate = $rate(1:$NTIME-2);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        if ($rate->sum >= -$prod_others_max) {
            $consumption_reaction_rates{"Consumption Others"} += $rate;
        } else {
            $consumption_reaction_rates{$reaction_string} += $rate;
        }
    } 
    remove_common_processes(\%production_reaction_rates, \%consumption_reaction_rates);
    
    return (\%production_reaction_rates, \%consumption_reaction_rates);
}

sub get_species { #get radicals from file
    my ($file) = @_; 

    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my %hash;
    my @lines = <$in>;
    close $in;
    my @separate = map { split /\s/, $_ } @lines;
    $hash{$_} += 1 foreach (@separate);
    return keys %hash;
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
