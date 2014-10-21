#! /usr/bin/env perl
# Show total Ox production versus Carbon number of emitted NMVOCs in all mechanisms
# Version 0: Jane Coates 25/6/2014

use strict;
use diagnostics;
use Statistics::R;

my $base_dir = "/work/users/jco/MECCA/Mechanism_Comparison/TOPP_sum_comparison";

#MCM v3.2
my $mcm_3_2_file = "$base_dir/MCMv3.2_TOPP_sums.txt";
my ($mcm_3_2_total) = get_final_TOPP_sum($mcm_3_2_file);
my ($mcm_3_2_carbons) = get_carbons('MCMv3.2');

#MCM v3.1
my $mcm_3_1_file = "$base_dir/MCMv3.1_TOPP_sums.txt";
my ($mcm_3_1_total) = get_final_TOPP_sum($mcm_3_1_file);
my ($mcm_3_1_carbons) = get_carbons('MCMv3.1');

#CRI v2
my $cri_file = "$base_dir/CRIv2_TOPP_sums.txt";
my ($cri_total) = get_final_TOPP_sum($cri_file);
my ($cri_carbons) = get_carbons('CRIv2');

#MOZART-4
my $mozart_file = "$base_dir/MOZART-4_TOPP_sums.txt";
my ($mozart_total) = get_final_TOPP_sum($mozart_file);
my ($mozart_carbons) = get_carbons('MOZART-4');

#RADM2
my $radm2_file = "$base_dir/RADM2_TOPP_sums.txt";
my ($radm2_total) = get_final_TOPP_sum($radm2_file);
my ($radm2_carbons) = get_carbons('RADM2');

#RACM
my $racm_file = "$base_dir/RACM_TOPP_sums.txt";
my ($racm_total) = get_final_TOPP_sum($racm_file);
my ($racm_carbons) = get_carbons('RACM');

plot({
        mcm3_2_data     => $mcm_3_2_total,
        mcm3_2_carbons  => $mcm_3_2_carbons,
        mcm3_1_data     => $mcm_3_1_total,
        mcm3_1_carbons  => $mcm_3_1_carbons,
        cri_data        => $cri_total,
        cri_carbons     => $cri_carbons,
        mozart_data     => $mozart_total,
        mozart_carbons  => $mozart_carbons,
        radm2_data      => $radm2_total,
        radm2_carbons   => $radm2_carbons,
        racm_data       => $racm_total,
        racm_carbons    => $racm_carbons,
});

sub plot {
    my ($args) = @_;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `,
            q` library(scales) `,
            q` library(grid) `,
            q` library(Cairo) `,
    );

    $R->run(q` Mechanism = c() `,
            q` VOC = c() `,
            q` Ox.prod = c() `,
            q` C.number = c() `,
    );

    #MCM v3.2
    $R->set('mcm3.2.rep', scalar keys %{$args->{mcm3_2_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("MCMv3.2", mcm3.2.rep)) `);
    foreach my $VOC ( sort keys %{$args->{mcm3_2_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{mcm3_2_data}}{$VOC});
        $R->set('c', ${$args->{mcm3_2_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    #MCM v3.1
    $R->set('mcm3.1.rep', scalar keys %{$args->{mcm3_1_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("MCMv3.1", mcm3.1.rep)) `);
    foreach my $VOC ( sort keys %{$args->{mcm3_1_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{mcm3_1_data}}{$VOC});
        $R->set('c', ${$args->{mcm3_1_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    #CRI v2
    $R->set('cri.rep', scalar keys %{$args->{cri_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("CRIv2", cri.rep)) `);
    foreach my $VOC ( sort keys %{$args->{cri_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{cri_data}}{$VOC});
        $R->set('c', ${$args->{cri_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    #MOZART-4
    $R->set('mozart.rep', scalar keys %{$args->{mozart_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("MOZART-4", mozart.rep)) `);
    foreach my $VOC ( sort keys %{$args->{mozart_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{mozart_data}}{$VOC});
        $R->set('c', ${$args->{mozart_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    #RADM2
    $R->set('radm2.rep', scalar keys %{$args->{radm2_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("RADM2", radm2.rep)) `);
    foreach my $VOC ( sort keys %{$args->{radm2_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{radm2_data}}{$VOC});
        $R->set('c', ${$args->{radm2_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    #RACM
    $R->set('racm.rep', scalar keys %{$args->{racm_data}});
    $R->run(q` Mechanism = c(Mechanism, rep("RACM", racm.rep)) `);
    foreach my $VOC ( sort keys %{$args->{racm_data}}) {
        $R->set('voc', $VOC);
        $R->set('Ox', ${$args->{racm_data}}{$VOC});
        $R->set('c', ${$args->{racm_carbons}}{$VOC});
        $R->run(q` VOC = cbind(VOC, voc) `,
                q` Ox.prod = cbind(Ox.prod, Ox) `,
                q` C.number = cbind(C.number, c) `,
        );
    }

    $R->run(q` VOC = c(VOC) `,
            q` Ox.prod = c(Ox.prod) `,
            q` C.number = c(C.number) `,
    );

    $R->run(q` data = data.frame(Mechanism, VOC, Ox.prod, C.number) `,
            q` data$VOC = as.character(data$VOC) `,
            q` data$VOC[data$VOC == "C2H6" | data$VOC == "C3H8" | data$VOC == "NC4H10" | data$VOC == "IC4H10" | data$VOC == "NC5H12" | data$VOC == "IC5H12" | data$VOC == "NC6H14" | data$VOC == "NC7H16" | data$VOC == "NC8H18" | data$VOC == "BIGALK" | data$VOC == "ETH" | data$VOC == "HC3" | data$VOC == "HC5" | data$VOC == "HC8"] = "Alkane" `,
            q` data$VOC[data$VOC == "C2H4" | data$VOC == "C3H6" | data$VOC == "BUT1ENE" | data$VOC == "MEPROPENE" | data$VOC == "C5H8" | data$VOC == "BIGENE" | data$VOC == "ISOP" | data$VOC == "OL2" | data$VOC == "OLI" | data$VOC == "OLT" | data$VOC == "ISO" | data$VOC == "ETE" ] = "Alkene" `,
            q` data$VOC[data$VOC == "BENZENE" | data$VOC == "TOLUENE" | data$VOC == "OXYL" | data$VOC == "PXYL" | data$VOC == "MXYL" | data$VOC == "EBENZ" | data$VOC == "TOL" | data$VOC == "XYL"] = "Aromatic" `,
            q` data$VOC = as.factor(data$VOC) `,
    );

    $R->run(q` plot.lines = function () { list(ylab("Molecules (Ox) / molecules (VOC)"),
                                               xlab("Carbon Number of VOC"),
                                               ggtitle("Total Ox Produced normalised by emissions versus Carbon number of Emitted VOC"),
                                               geom_point(size = 13), 
                                               scale_x_continuous(limits = c(1, 9), breaks = seq(1, 9, 1)), 
                                               theme_bw(), 
                                               facet_wrap( ~ VOC),
                                               theme(axis.title.x = element_text(size = 30, face = "bold")), 
                                               theme(axis.title.y = element_text(size = 30, face = "bold")), 
                                               theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), legend.title = element_text(size = 25, face = "bold"), legend.key.size = unit(3, "cm"), legend.text = element_text(size = 20), legend.key = element_blank(), plot.title = element_text(size = 40, face = "bold"), strip.text.x = element_text(size = 30, face = "bold")) ) } `);
    
    $R->run(q` plot = ggplot(data = data, aes(x = C.number, y = Ox.prod, colour = Mechanism)) `,
            q` plot = plot + plot.lines() `,
            q` CairoPNG(file = "Ox_prod_vs_C.png", width = 2000, height = 1000) `,
            q` print(plot) `,
            q` dev.off() `,
    );

    $R->stop;
}

sub get_final_TOPP_sum {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    my @lines = <$in>;
    close $in;
    
    my %total;
    foreach my $line (@lines) {
        chomp $line;
        my ($VOC, $data) = split / => /, $line; 
        my @TOPP_sums = split /\s/, $data;
        $total{$VOC} = $TOPP_sums[6];
    }
    return \%total;
}

sub get_carbons {
    my ($mech) = @_;

    my %carbons;
    if ($mech eq "MCMv3.2" or $mech eq "MCMv3.1" or $mech eq "CRIv2" or $mech eq "CBM-IV" or $mech eq "CB05") {
        $carbons{'C2H6'} = 2;
        $carbons{'C3H8'} = 3;
        $carbons{'CH4'} = 1;
        $carbons{'NC4H10'} = 4;
        $carbons{'IC4H10'} = 4;
        $carbons{'NC5H12'} = 5;
        $carbons{'IC5H12'} = 5;
        $carbons{'NC6H14'} = 6;
        $carbons{'NC7H16'} = 7;
        $carbons{'NC8H18'} = 8;
        $carbons{'C2H4'} = 2;
        $carbons{'C3H6'} = 3;
        $carbons{'BUT1ENE'} = 4;
        $carbons{'MEPROPENE'} = 4;
        $carbons{'BENZENE'} = 6;
        $carbons{'TOLUENE'} = 7;
        $carbons{'OXYL'} = 8;
        $carbons{'PXYL'} = 8;
        $carbons{'MXYL'} = 8;
        $carbons{'EBENZ'} = 8;
        $carbons{'C5H8'} = 5;
    } elsif ($mech eq "MOZART-4") {
        $carbons{'C2H6'} = 2;
        $carbons{'C3H8'} = 3;
        $carbons{'BIGALK'} = 5;
        $carbons{'C2H4'} = 2;
        $carbons{'C3H6'} = 3;
        $carbons{'BIGENE'} = 4;
        $carbons{'ISOP'} = 5;
        $carbons{'TOLUENE'} = 7;
    } elsif ($mech eq "RADM2") {
        $carbons{'ETH'} = 2;
        $carbons{'HC3'} = 2.9;
        $carbons{'HC5'} = 4.8;
        $carbons{'HC8'} = 7.9;
        $carbons{'OL2'} = 2;
        $carbons{'OLI'} = 3.8;
        $carbons{'OLT'} = 4.8;
        $carbons{'ISO'} = 5;
        $carbons{'TOL'} = 7.1;
        $carbons{'XYL'} = 8.9;
    } elsif ($mech eq "RACM") {
        $carbons{'ETH'} = 2;
        $carbons{'HC3'} = 2.9;
        $carbons{'HC5'} = 4.8;
        $carbons{'HC8'} = 7.9;
        $carbons{'ETE'} = 2;
        $carbons{'OLI'} = 5;
        $carbons{'OLT'} = 3.8;
        $carbons{'ISO'} = 5;
        $carbons{'TOL'} = 7.1;
        $carbons{'XYL'} = 8.9;
    }
    return \%carbons;
}
