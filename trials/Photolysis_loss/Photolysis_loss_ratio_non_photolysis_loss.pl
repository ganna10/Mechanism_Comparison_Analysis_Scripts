#!/ usr/bin/env perl
# Plot ratio of photolytic loss rates to non-photolytic loss rates
# Version 0: Jane Coates 28/11/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCM_3.2_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(2:$ntime-3);

my @species = qw( HCHO HCHO HCHO CH2O HCHO HCHO HCHO HCHO FORM );
my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my $index = 0;
my (%data);

foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $spc_file = "$base/$run/gas.spc";
    my ($photolysis, $non_photolysis);
    my $tagged_species = get_tagged_species($species[$index], $spc_file);
    foreach my $species (@$tagged_species) {
        my $consumers = $kpp->consuming($species);
        foreach my $reaction (@$consumers) {
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number);
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            if ($reaction_string =~ /hv/) {
                $photolysis += $rate(2:$ntime-3);
            } else {
                $non_photolysis += $rate(2:$ntime-3);
            }
        }
    }
    $data{$mechanisms[$index]} = $photolysis / $non_photolysis;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [map {$_} $times->dog]);
$R->run(q` data = data.frame(Time) `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->set('ratio', [map { $_ } $data{$mechanism}->dog]);
    $R->run(q` data[mechanism] = ratio `);
}
$R->run(q` data = gather(data, Mechanism, Ratio, -Time) `);
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HCHO_photolysis_non_photolysis_ratio.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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
