#! /usr/bin/env perl
# Compare total loss through photolysis of MCM species entered as ARGV, if it is represented in other mechanisms. Compare mixing ratios and difference between maxima.
# Version 0: Jane Coates 27/11/2014

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

my $mcm_species = $ARGV[0];
die "Need input species\n" unless (defined $mcm_species);
my ($mech_species, $runs, $mechanism) = get_species_data($mcm_species);
my $index = 0;
my (%rates, %mixing_ratio);

foreach my $run (@$runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $spc_file = "$base/$run/gas.spc";
    my $tagged_species = get_tagged_species($mech_species->[$index], $spc_file);
    foreach my $tagged (@$tagged_species) {
        my $mixing_ratio = $mecca->tracer($tagged) * 1e9;
        $mixing_ratio{$mechanism->[$index]} += $mixing_ratio(1:$ntime-2);
        my $reactions = $kpp->reacting_with($tagged, 'hv');
        foreach my $reaction (@$reactions) {
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number);
            $rates{$mechanism->[$index]} += $rate(1:$ntime-2);
        }
    }
    $index++;
}

my (%photolysis_max, %ratios_max); 
my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [map {$_} $times->dog]);
$R->run(q` rate.data = data.frame(Time) `);
foreach my $mechanism (sort keys %rates) {
    $photolysis_max{$mechanism} = $rates{$mechanism}->max;
    $R->set('mechanism', $mechanism);
    $R->set('rate', [map {$_} $rates{$mechanism}->dog]);
    $R->run(q` rate.data[mechanism] = rate `),
}
$R->run(q` rate.data = gather(rate.data, Mechanism, Rate, -Time) `);
#my $p = $R->run(q` print(rate.data) `);
#print $p, "\n";
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot.lines = function () { list( geom_line() ,
                                            theme_bw(),
                                            xlab("Time (days)"),
                                            theme(axis.title.x = element_text(face = "bold")),
                                            theme(panel.grid = element_blank()),
                                            scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)),
                                            theme(plot.title = element_text(face = "bold")) ,
                                            theme(legend.position = c(0.99, 0.99)),
                                            theme(legend.justification = c(0.99, 0.99)),
                                            theme(legend.title = element_blank()) ,
                                            theme(legend.key = element_blank()) ) } `);

$R->set('rate.filename', "${mcm_species}_photolysis_loss.pdf");
$R->set('rate.title', "$mcm_species Photolysis Loss Rate");
$R->run(q` rate.plot = ggplot(rate.data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` rate.plot = rate.plot + ggtitle(rate.title) `,
        q` rate.plot = rate.plot + ylab(expression(bold(paste("Photolysis Rate (", s^-1, ")")))) `,
        q` rate.plot = rate.plot + scale_colour_manual(values = my.colours) `,
        q` rate.plot = rate.plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = rate.filename) `,
        q` print(rate.plot) `,
        q` dev.off() `,
);

$R->run(q` mix.data = data.frame(Time) `);
foreach my $mechanism (sort keys %mixing_ratio) { 
    $ratios_max{$mechanism} = $mixing_ratio{$mechanism}->max;
    $R->set('mechanism', $mechanism);
    $R->set('mixing.ratio', [map {$_} $mixing_ratio{$mechanism}->dog]);
    $R->run(q` mix.data[mechanism] = mixing.ratio `);
}
$R->run(q` mix.data = gather(mix.data, Mechanism, Mixing.Ratio, -Time) `);
$R->set('mix.filename', "${mcm_species}_mixing_ratio.pdf");
$R->set('mix.title', "$mcm_species Mixing Ratios");
$R->run(q` mix.plot = ggplot(mix.data, aes( x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` mix.plot = mix.plot + ggtitle(mix.title) `,
        q` mix.plot = mix.plot + ylab("Mixing Ratio (ppbv)") `,
        q` mix.plot = mix.plot + scale_colour_manual(values = my.colours) `,
        q` mix.plot = mix.plot + plot.lines() `,
        q` mix.plot = mix.plot + theme(axis.title.y = element_text(face = "bold")) `,
);

$R->run(q` CairoPDF(file = mix.filename) `,
        q` print(mix.plot) `,
        q` dev.off() `,
);
$R->stop();

my $max_file = "${mcm_species}_max_differences.txt";
open my $out, '>:encoding(utf-8)', $max_file or die "Can't open $max_file : $!";
print $out "$mcm_species Maxima Differences\n";
print $out "\nPhotolysis Loss Rates\n";

my @sorted_photo_mechanisms = sort { $photolysis_max{$b} <=> $photolysis_max{$a} } keys %photolysis_max;
my @sorted_photo_max = @photolysis_max{@sorted_photo_mechanisms};
my $final = @sorted_photo_mechanisms - 1;
print $out "$sorted_photo_mechanisms[$_] => $sorted_photo_max[$_] s^-1\n" foreach (0..$final);
print $out "\nMaxima Difference = ", $sorted_photo_max[0] - $sorted_photo_max[$final], " s^-1\n";

print $out "\nMixing Ratios\n";
my @sorted_ratio_mechanisms = sort { $ratios_max{$b} <=> $ratios_max{$a} } keys %ratios_max;
my @sorted_ratio_max = @ratios_max{@sorted_ratio_mechanisms};
print $out "$sorted_ratio_mechanisms[$_] => $sorted_ratio_max[$_] ppbv\n" foreach (0..$final);
print $out "\nMaxima Difference = ", $sorted_ratio_max[0] - $sorted_ratio_max[$final], " ppbv\n";

close $out;

sub get_tagged_species {
    my ($species, $spc_file) = @_;
    open my $file, '<:encoding(utf8)', $spc_file or die "Can't open $spc_file : $!";
    my @lines = <$file>;
    close $file;
    my @tagged;
    foreach my $line (@lines) {
        next unless ($line =~ /^${species}_/);
        $line =~ s/^${species}_(.*?)\b.*\n$/${species}_$1/;
        push @tagged, $line;
    }
    return \@tagged;
}

sub get_species_data {
    my ($mcm) = @_;
    my (@species, @runs, @mechanisms);

    if ($mcm eq "HCHO") {
        @species = qw( HCHO HCHO HCHO CH2O HCHO HCHO HCHO HCHO FORM );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } elsif ($mcm eq "CH3CHO") {
        @species = qw( CH3CHO CH3CHO CH3CHO CH3CHO ALD ALD ACD ALD2 ALD2 );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } elsif ($mcm eq "CH3OOH") {
        @species = qw( CH3OOH CH3OOH CH3OOH CH3OOH OP1 OP1 OP1 MEPX );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CB05 );
    } elsif ($mcm eq "C2H5OOH") {
        @species = qw( C2H5OOH C2H5OOH C2H5OOH C2H5OOH OP2 OP2 OP2 ROOH );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CB05 );
    } elsif ($mcm eq "CH3CO3H") {
        @species = qw( CH3CO3H CH3CO3H CH3CO3H CH3COOOH PAA PAA PAA PACD );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CB05 );
    } elsif ($mcm eq "MACR") {
        @species = qw( MACR MACR MACR MACR MACR MACR );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RACM_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RACM RACM2 );
    } elsif ($mcm eq "MVK") {
        @species = qw( MVK MVK MVK MVK MVK );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RACM2 );
    } elsif ($mcm eq "CH3COCH3") {
        @species = qw( CH3COCH3 CH3COCH3 CH3COCH3 CH3COCH3 ACT );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RACM2 );
    } elsif ($mcm eq "MGLYOX") {
        @species = qw( MGLYOX MGLYOX CARB6 CH3COCHO MGLY MGLY MGLY MGLY MGLY );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } elsif ($mcm eq "GLYOX") {
        @species = qw( GLYOX GLYOX CARB3 GLYOXAL GLY GLY GLY );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } elsif ($mcm eq "MEK") {
        @species = qw( MEK MEK MEK MEK KET KET MEK );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } else {
        print "Nothing found for $mcm\n";
    }
    return (\@species, \@runs, \@mechanisms);
}
