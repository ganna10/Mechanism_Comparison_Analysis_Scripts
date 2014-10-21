#! /usr/bin/env perl
# Plot rate time series of HCHO + hv = CO + 2 HO2 reactions
# Version 0: Jane Coates 18/9/2014

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/MECCA";
my %rates;

my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
$times = $times(1:$NTIME-2);
$times -= $times->at(0);
#$times /= 86400; #for time series non-bar plot
$times /= 3600; #for bar plot
my @time_axis = map { $_ } $times->dog;
my @time_blocks;
foreach my $time (@time_axis) {
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

my @dirs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my @mech_species = qw( HCHO HCHO HCHO CH2O HCHO HCHO HCHO HCHO FORM );
#my @dirs = qw( CB05_tagging );
#my @mechanisms = qw( CB05 );
my $index = 0;

foreach my $run (@dirs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $spc_file = "$base/$run/gas.spc";
    my $tagged_mechanism_species = get_tagged_species($mech_species[$index], $spc_file);
    $rates{$mechanisms[$index]} = get_rates($tagged_mechanism_species, $mecca, $kpp);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

#$R->set('Time', [@time_axis]); #for non-bar plot
$R->set('Time', [@time_blocks]); #for bar_plot
$R->run(q` data = data.frame(Time) `);

foreach my $run (sort keys %rates) {
    $R->set('Mechanism', $run); 
    $R->set('rates', [@{$rates{$run}}]);
    $R->run(q` data[Mechanism] = rates `);
}

$R->run(q` data["CB05"] = data["MCMv3.2"] - data["CB05"] `,
        q` data["CBM-IV"] = data["MCMv3.2"] - data["CBM-IV"] `,
        q` data["CRIv2"] = data["MCMv3.2"] - data["CRIv2"] `,
        q` data["MCMv3.1"] = data["MCMv3.2"] - data["MCMv3.1"] `,
        q` data["MOZART-4"] = data["MCMv3.2"] - data["MOZART-4"] `,
        q` data["RACM"] = data["MCMv3.2"] - data["RACM"] `,
        q` data["RACM2"] = data["MCMv3.2"] - data["RACM2"] `,
        q` data["RADM2"] = data["MCMv3.2"] - data["RADM2"] `,
        q` data["MCMv3.2"] = NULL `,
        q` data = ddply(data, .(Time), colwise(sum)) `,
        q` data = data[1:7,] `,
        q` data = melt(data, id.vars = c("Time"), variable.name = "Mechanism", value.name = "rate") `, #arrange data
);
#my $p = $R->run(q` print(head(data)) `);
#print $p, "\n";

$R->set('plot.title', "HCHO + hv = CO + HO2 Rate Comparison");
$R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#662c91", "#12b2b2", "#b33893", "#a11d20") `); #"#ed2d2e", 
$R->run(q` plot = ggplot(data = data, aes(x = Time, y = rate, colour = Mechanism, group = Mechanism)) `, #plotting
        q` plot = plot + geom_line(size = 1) `,
        q` plot = plot + geom_point(size = 3) `,
        q` plot = plot + ggtitle(plot.title) `,
        #q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab(expression(bold(paste("Photolysis Rate(", s^-1, ")")))) `,
        #q` plot = plot + scale_x_continuous(limits = c(0,7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + theme_bw() `,
        #q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(0.99, 0.09)) `,
        q` plot = plot + theme(legend.justification = c(0.99, 0.09)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->set('filename', "HCHO_HO2_photolysis_rate_difference_overall.pdf");
$R->run(q` CairoPDF(file = filename, width = 10, height = 7) `, #print to file
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

sub get_rates {
    my ($all_species, $mecca, $kpp) = @_;

    my $rate;
    foreach my $species (@$all_species) {
        my $reaction = $kpp->reacting_with($species, 'hv');
        foreach (@$reaction) {
            my $reaction_string = $kpp->reaction_string($_);
            next if ($reaction_string =~ /H2\b/);
            my $r_number = $kpp->reaction_number($_);
            $rate += $mecca->rate($r_number);
        }
    }
    $rate = $rate(1:$NTIME-2);
    my @rates = map { $_ } $rate->dog;
    return \@rates;
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
