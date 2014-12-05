#! /usr/bin/env perl
# get RCONST of each photolysis reaction for MCM species (ARGV)
# Version 0: Jane Coates 3/12/2014

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mcm_species = $ARGV[0];
die "Need input species\n" unless (defined $mcm_species);
my ($mech_species, $runs, $mechanism) = get_species_data($mcm_species);
my $index = 0;
my %data;

my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

foreach my $run (@$runs) { 
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $spc_file = "$base/$run/gas.spc";
    my $tagged_species = get_tagged_species($mech_species->[$index], $spc_file);
    my $number_species = scalar @$tagged_species;
    foreach my $tagged (@$tagged_species) {
        my $reactions = $kpp->reacting_with($tagged, 'hv');
        foreach my $reaction (@$reactions) {
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rconst = $mecca->rconst($reaction_number);
            $data{$mechanism->[$index]} = $rconst(1:$ntime-2);
        }
    }
    #$data{$mechanism->[$index]} /= $number_species;
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame(Time) `);
$R->set('filename', "${mcm_species}_photolysis_loss_frequency.pdf");
$R->set('title', "$mcm_species Photolysis Loss Frequency");
$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#b569b3", "CRIv2" = "#ef6638", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
foreach my $mechanism (sort keys %data) {
    next if ($mechanism eq "RACM");
    $R->set('mechanism', $mechanism);
    $R->set('frequency', [map { $_ } $data{$mechanism}->dog]);
    $R->run(q` data[mechanism] = frequency `);
}
#my $p = $R->run(q` print(head(data)) `);
#print "$p\n";
$R->run(q` data = gather(data, Mechanism, Frequency, -Time) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Frequency, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + ggtitle(title) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Photolysis Frequency (s^-1)") `,
        q` plot = plot + xlab("Time (Days)") `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = filename) `,
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
        @mechanisms = qw( MCMv3.2 MCMv3.1 MOZART-4 RADM2 RACM RACM2 CB05 );
    } elsif ($mcm eq "MACR") {
        @species = qw( MACR MACR MACR MACR MACR );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates MOZART_tagging RACM_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 MOZART-4 RACM RACM2 );
    } elsif ($mcm eq "MVK") {
        @species = qw( MVK MVK MVK MVK );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates MOZART_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 MOZART-4 RACM2 );
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
        @species = qw( MEK MEK CARB11A MEK KET KET MEK );
        @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged );
        @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
    } else {
        print "Nothing found for $mcm\n";
    }
    return (\@species, \@runs, \@mechanisms);
}
