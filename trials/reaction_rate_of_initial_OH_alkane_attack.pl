#! /usr/bin/env perl
# plot reaction rates of OH + alkane every mechanism
# Version 0: Jane Coates 6/3/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 );
my @alkanes = qw( Ethane Propane Butane 2-Methylpropane Pentane 2-Methylbutane Hexane Heptane Octane );
my %data;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn);
    foreach my $alkane (@alkanes) {
        my $species = get_species($alkane, $mechanism);
        my $reaction = $kpp->reacting_with($species, "OH");
        my $reaction_number = $kpp->reaction_number($reaction->[0]);
        $data{$mechanism}{$alkane} = $mecca->rate($reaction_number);
    }
}

###fractional_contributions of lumped species
foreach my $mechanism (sort keys %data) {
    if ($mechanism eq "RACM2") {
        $data{$mechanism}{"Propane"} *= 0.628;
        $data{$mechanism}{"Butane"} *= 0.243;
        $data{$mechanism}{"2-Methylpropane"} *= 0.129;
        $data{$mechanism}{"Pentane"} *= 0.264;
        $data{$mechanism}{"2-Methylbutane"} *= 0.615;
        $data{$mechanism}{"Hexane"} *= 0.086;
        $data{$mechanism}{"Heptane"} *= 0.035; 
    } elsif ($mechanism =~ /RA/) {
        $data{$mechanism}{"Propane"} *= 0.628;
        $data{$mechanism}{"Butane"} *= 0.243;
        $data{$mechanism}{"2-Methylpropane"} *= 0.129;
        $data{$mechanism}{"Pentane"} *= 0.264;
        $data{$mechanism}{"2-Methylbutane"} *= 0.615;
        $data{$mechanism}{"Hexane"} *= 0.086;
        $data{$mechanism}{"Heptane"} *= 0.035;
    } elsif ($mechanism =~ /MO/) { 
        $data{$mechanism}{"Butane"} *= 0.285;
        $data{$mechanism}{"2-Methylpropane"} *= 0.151;
        $data{$mechanism}{"Pentane"} *= 0.146;
        $data{$mechanism}{"2-Methylbutane"} *= 0.340;
        $data{$mechanism}{"Hexane"} *= 0.048;
        $data{$mechanism}{"Heptane"} *= 0.020;
        $data{$mechanism}{"Octane"} *= 0.010;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `);
    foreach my $alkane (sort keys %{$data{$mechanism}}) {
        $R->set('alkane', $alkane);
        $R->set('rate', [ map { $_ } $data{$mechanism}{$alkane}->dog ]);
        $R->run(q` pre[alkane] = rate `);
    }
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, Alkane, Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(pre) `);
#print $p, "\n";

$R->run(q` my.colours = c("CB05" = "#0352cb", "CBM-IV" = "#ef6638", "CRIv2" = "#b569b3", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Alkane, scales = "free")`,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "Alkanes_initial_OH_reaction_rates.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_species {
    my ($alkane, $mechanism) = @_;
    my $species;
    if ($alkane eq "Ethane") {
        if ($mechanism =~ /MC|CRI|MO/) {
            $species = "C2H6";
        } elsif ($mechanism =~ /RA/) {
            $species = "ETH";
        } elsif ($mechanism =~ /CB0/) {
            $species = "ETHA_C2H6";
        } elsif ($mechanism =~ /CBM/) {
            $species = "PAR_C2H6";
        }
    } elsif ($alkane eq "Propane") {
        if ($mechanism =~ /MC|CRI|MO/) {
            $species = "C3H8";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC3";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_C3H8";
        }
    } elsif ($alkane eq "Butane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "NC4H10";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_NC4H10";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC3";
        }
    } elsif ($alkane eq "2-Methylpropane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "IC4H10";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_IC4H10";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC3";
        }
    } elsif ($alkane eq "Pentane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "NC5H12";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_NC5H12";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC5";
        }
    } elsif ($alkane eq "2-Methylbutane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "IC5H12";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_IC5H12";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC5";
        }
    } elsif ($alkane eq "Hexane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "NC6H14";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_NC6H14";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC5";
        }
    } elsif ($alkane eq "Heptane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "NC7H16";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_NC7H16";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC5";
        }
    } elsif ($alkane eq "Octane") {
        if ($mechanism =~ /MC|CRI/) {
            $species = "NC8H18";
        } elsif ($mechanism =~ /CB/) {
            $species = "PAR_NC8H18";
        } elsif ($mechanism =~ /MO/) {
            $species = "BIGALK";
        } elsif ($mechanism =~ /RA/) {
            $species = "HC5";
        }
    } elsif ($alkane eq "Ethene") {
        if ($mechanism =~ /MC|CRI|MO/) {
            $species = "C2H4";
        } elsif ($mechanism =~ /CB/) {
            $species = "ETH_C2H4";
        } elsif ($mechanism =~ /RAD/) {
            $species = "OL2";
        } elsif ($mechanism =~ /RAC/) {
            $species = "ETE";
        }
    } elsif ($alkane eq "Toluene") {
        if ($mechanism =~ /MC|CRI|MO/) {
            $species = "TOLUENE";
        } elsif ($mechanism =~ /RA/) {
            $species = "TOL";
        } elsif ($mechanism =~ /CB/) {
            $species = "TOL_TOLUENE";
        }
    } else {
        print "No data for $alkane\n";
    }
    return $species;
}
