#! /usr/bin/env perl
# Compare % of unreacted alkanes after first day
# Version 0: Jane Coates 5/3/2015

use strict;
use diagnostics;
use PDL;
use MECCA;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my @alkanes = qw( Ethane Propane Butane 2-Methylpropane Pentane 2-Methylbutane Hexane Heptane Octane Ethene Toluene );
my %data;

my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    foreach my $alkane (@alkanes) {
        my $species = get_species($alkane, $mechanism);
        $data{$mechanism}{$alkane} = $mecca->tracer($species);
    }
}

###fractional_contributions of lumped species
foreach my $mechanism (sort keys %data) {
    if ($mechanism eq "RACM2") {
        $data{$mechanism}{"Propane"} *= 0.559;
        $data{$mechanism}{"Butane"} *= 0.228;
        $data{$mechanism}{"2-Methylpropane"} *= 0.153;
        $data{$mechanism}{"Pentane"} *= 0.256;
        $data{$mechanism}{"2-Methylbutane"} *= 0.596;
        $data{$mechanism}{"Hexane"} *= 0.100;
        $data{$mechanism}{"Heptane"} *= 0.048;
        $data{$mechanism}{"Toluene"} *= 0.852;
    } elsif ($mechanism =~ /RA/) {
        $data{$mechanism}{"Propane"} *= 0.559;
        $data{$mechanism}{"Butane"} *= 0.228;
        $data{$mechanism}{"2-Methylpropane"} *= 0.153;
        $data{$mechanism}{"Pentane"} *= 0.256;
        $data{$mechanism}{"2-Methylbutane"} *= 0.596;
        $data{$mechanism}{"Hexane"} *= 0.100;
        $data{$mechanism}{"Heptane"} *= 0.048;
        $data{$mechanism}{"Toluene"} *= 0.679;
    } elsif ($mechanism =~ /MO/) { 
        $data{$mechanism}{"Butane"} *= 0.244;
        $data{$mechanism}{"2-Methylpropane"} *= 0.129;
        $data{$mechanism}{"Pentane"} *= 0.156;
        $data{$mechanism}{"2-Methylbutane"} *= 0.364;
        $data{$mechanism}{"Hexane"} *= 0.061;
        $data{$mechanism}{"Heptane"} *= 0.029;
        $data{$mechanism}{"Octane"} *= 0.017;
        $data{$mechanism}{"Toluene"} *= 0.465;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `, 
        q` library(Cairo) `,
        q` library(scales) `,
);

$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $alkane (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Mechanism = c(mechanism)) `);
        $R->set('alkane', $alkane);
        $R->set('initial.mr', [ $data{$mechanism}{$alkane}->max ]);
        $R->set('day2.mr', [ $data{$mechanism}{$alkane}->at(73) ]);
        $R->run(q` unreacted = day2.mr / initial.mr `,
                q` pre$Alkane = alkane `,
                q` pre$Initial.MR = initial.mr `,
                q` pre$Day2.MR = day2.mr `,
                q` pre$Unreacted = unreacted `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Alkane = factor(data$Alkane, levels = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Toluene")) `,
        q` data$Mechanism = factor(data$Mechanism, levels = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05"))) `,
);

$R->run(q` plot = ggplot(data, aes(y = Unreacted, x = Mechanism)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + facet_wrap( ~ Alkane )`,
        q` plot = plot + scale_y_continuous(labels = percent) `,
);

$R->run(q` CairoPDF(file = "alkane_unreacted_after_day1.pdf") `,
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
