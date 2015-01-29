#! /usr/bin/env perl
# Plot mixing ratio time series of each radical in the total radical family for CBM-IV, CB05 and MCM v3.2
# Version 0: Jane Coates 29/1/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

my @mechanisms = qw( MCMv3.2 CBM-IV CB05 );
#my @mechanisms = qw( CB05 );
my %data;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $radicals_file = "$base/${mechanism}_tagged/radicals.txt";
    my $radicals = read_file($radicals_file);
    $radicals =~ s/\s+/ /g;
    my @radicals = split / /, $radicals;
    my %mixing_ratios;
    foreach my $radical (@radicals) {
        my $mixing_ratio = $mecca->tracer($radical);
        next unless (defined $mixing_ratio);
        $mixing_ratio = $mixing_ratio(1:$ntime-2) *1e9;
        (my $label = $radical) =~ s/_(.*?)\b//g;
        $mixing_ratios{$label} += $mixing_ratio;
    }
    $data{$mechanism} = sort_hash(\%mixing_ratios);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$data{$mechanism}}) {
        foreach my $radical (sort keys %$ref) {
            $R->set('radical', $radical);
            $R->set('mixing.ratio', [ map { $_ } $ref->{$radical}->dog ]);
            $R->run(q` pre[radical] = mixing.ratio `);
        }
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Radical, Mixing.Ratio, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, fill = Radical, group = Radical)) `,
        q` plot = plot + geom_area(position = "stack") `,
        q` plot = plot + facet_wrap( ~ Mechanism) `,
);

$R->run(q` CairoPDF(file = "Contributions_to_Radical_MRs.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub sort_hash {
    my ($hash) = @_;
    my $others = 0.03;
    foreach my $key (keys %$hash) {
        if ($hash->{$key}->sum < $others) {
            $hash->{"Others"} += $hash->{$key};
            delete $hash->{$key};
        }
    }
    my $sort_function = sub { $_[0]->sum };
    my @sorted = sort { &$sort_function($hash->{$b} <=> $hash->{$a}) } keys %$hash;
    my @final_sorted;
    foreach (@sorted) {
        next if ($_ =~ /Others/);
        push @final_sorted, { $_ => $hash->{$_} };
    }
    push @final_sorted, { "Others" => $hash->{"Others"} };
    return \@final_sorted;
}

sub read_file {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    close $in;
    return $data;
}
