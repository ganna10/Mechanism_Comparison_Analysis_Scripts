#! /usr/bin/env perl
# Compare total reactive carbon loss during pentane and toluene degradation by mixing ratios of degradation products
# Version 0: Jane Coates 27/2/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/MCMv3.2_tagged/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = ( "MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2",  "CBM-IV", "CB05" );
#my @mechanisms = qw( RADM2 RACM );
my (%n_carbon, %data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel); 
    my $spcfile = "$base/${mechanism}_tagged/gas.spc";
    my $species = get_all_species($spcfile);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $n_carbon{$mechanism} = get_carbons($mechanism, $carbon_file);
    my @VOCs = qw( Pentane Toluene );
    foreach my $NMVOC (@VOCs) {
        my $parent = get_mechanism_species($NMVOC, $mechanism);
        ($data{$mechanism}{$NMVOC}) = get_data($mecca, $species, $mechanism, $n_carbon{$mechanism}, $parent);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `);
        $R->set('voc', $VOC);
        foreach my $carbon (sort keys %{$data{$mechanism}{$VOC}}) {
            $R->set('carbon', $carbon);
            $R->set('mixing.ratio', [map { $_ } $data{$mechanism}{$VOC}{$carbon}->dog]);
            $R->run(q` pre[carbon] = mixing.ratio `);
        }
        $R->run(q` pre$VOC = rep(voc, length(Time)) `);
        $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
                q` pre = gather(pre, Carbon, Mixing.Ratio, -Time, -VOC, -Mechanism) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$VOC = factor(data$VOC, levels = c("Pentane", "Toluene")) `);
$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
$R->run(q` my.names = c("C8" = "C8 ", "C7" = "C7 ", "C6" = "C6 ", "C5" = "C5 ", "C4" = "C4 ", "C3" = "C3 ", "C2" = "C2 ", "C1" = "C1 ") `);
$R->run(q` data$Carbon = factor(data$Carbon, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) `);

$R->run(q` plot = ggplot(data, aes(y = Mixing.Ratio, x = Mechanism, fill = Carbon)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + facet_grid( Time ~ VOC ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")), expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + ylab("Day-time O3 Mixing Ratios of all Products attributed to Carbon Number of Degradation Products\n(ppbv)") `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
        q` plot = plot + theme(strip.text.y = element_text(face = "bold", angle = 0)) `,
        q` plot = plot + theme(strip.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "Mixing_ratios_for_total_reactive_carbon.pdf", width = 8.7, height = 10) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $species, $mechanism, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    my %mixing_ratios;

    foreach my $species (@$species) {
        next unless ($species =~ /$VOC/);
        my $carbon;
        (my $lookup = $species) =~ s/_(.*?)\b//;
        if (exists $carbons{$lookup}) {
            $carbon = $carbons{$lookup};
        } else {
            print "No C found for species: $species\n";
        }
        if (defined $mecca->tracer($species)) {
            my $mixing_ratio = $mecca->tracer($species) * 1e9;
            $mixing_ratios{"C$carbon"} += $mixing_ratio(1:$NTIME-2) * $carbon;
        }
    }
    delete $mixing_ratios{"C0"};

    #add carbon numbers to nearest integral carbon number
    $mixing_ratios{"C2"} += $mixing_ratios{"C2.4"} if (exists $mixing_ratios{"C2.4"});
    delete $mixing_ratios{"C2.4"};
    $mixing_ratios{"C3"} += $mixing_ratios{"C2.9"} if (exists $mixing_ratios{"C2.9"});
    delete $mixing_ratios{"C2.9"};
    $mixing_ratios{"C4"} += $mixing_ratios{"C3.5"} if (exists $mixing_ratios{"C3.5"});
    delete $mixing_ratios{"C3.5"};
    $mixing_ratios{"C4"} += $mixing_ratios{"C3.9"} if (exists $mixing_ratios{"C3.9"});
    delete $mixing_ratios{"C3.9"};
    $mixing_ratios{"C4"} += $mixing_ratios{"C4.2"} if (exists $mixing_ratios{"C4.2"});
    delete $mixing_ratios{"C4.2"};
    $mixing_ratios{"C5"} += $mixing_ratios{"C4.8"} if (exists $mixing_ratios{"C4.8"});
    delete $mixing_ratios{"C4.8"};
    $mixing_ratios{"C7"} += $mixing_ratios{"C6.6"} if (exists $mixing_ratios{"C6.6"});
    delete $mixing_ratios{"C6.6"};
    $mixing_ratios{"C7"} += $mixing_ratios{"C7.1"} if (exists $mixing_ratios{"C7.1"});
    delete $mixing_ratios{"C7.1"};
        
    foreach my $carbon (sort keys %mixing_ratios) {
        if ($VOC =~ /TOL/) {
            if ($mechanism eq "RACM2") {
                $mixing_ratios{$carbon} *= 0.852;
            } elsif ($mechanism =~ /RA/) {
                $mixing_ratios{$carbon} *= 0.679;
            } elsif ($mechanism =~ /MOZ/) {
                $mixing_ratios{$carbon} *= 0.465;
            }
        } else { #pentane
            if ($mechanism =~ /RA/) {
                $mixing_ratios{$carbon} *= 0.256;
            } elsif ($mechanism =~ /MO/) {
                $mixing_ratios{$carbon} *= 0.156;
            }
        }
        my $reshape = $mixing_ratios{$carbon}->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $mixing_ratios{$carbon} = $integrate(0:13:2);
    }

    return \%mixing_ratios;
} 

sub get_carbons {
    my ($run, $file) = @_;
    my $carbons;
    if ($run =~ /MCMv3\.1|MCMv3\.2/) {
        $carbons = mcm_n_carbon($file);
    } elsif ($run =~ /MOZART/) {
        $carbons = mozart_n_carbon($file);
    } elsif ($run =~ /CRI|RADM2|RACM|CB/) {
        $carbons = carbons_others($file);
    } else {
        print "$run doesn't match\n";
    }
    return $carbons;
}

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub mozart_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my $words = join ',', (<$in>);
    close $in;
    my ($string) = $words =~ /Solution(.*?)End\sSolution/s;
    $string =~ s/^\s+,//;
    $string =~ s/,\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/RO2(.*?)->//;
    $string =~ s/ROOH(.*?)->//;
    my @species = split ',', $string;
    my %carbons;
    foreach my $species (@species) {
        $species =~ s/^\s+|\s+$//g;
        my $C_number = 0;
        if ($species !~ /->/ and $species !~ /(C[0-9])/) {
            $C_number ++ while ($species =~ m/C/g);
            $carbons{$species} = $C_number;
        } elsif ($species !~ /->/ and $species =~ /(C[0-9])/) { 
            my ($c_nr) = $species =~ /(C[0-9]+)/s;
            $c_nr =~ s/C//; 
            $C_number = $c_nr;
            $carbons{$species} = $C_number;
        } else {
            my ($mech, $molecule) = split ' -> ', $species;
            $mech =~ s/^\s+|\s+$//g;
            if ($molecule =~ /(C[0-9]+)/) { 
                my ($c_nr) = $molecule =~ /(C[0-9]+)/s;
                $c_nr =~ s/C//; 
                $C_number = $c_nr;
                $carbons{$mech} = $C_number;
            } else {
                $C_number ++ while ($molecule =~ m/C/g);
                $carbons{$mech} = $C_number;
            }
        }
    } 
    return \%carbons;
}

sub carbons_others { #get C-number from file names that have species and C# separated by space
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Cannot open file $file: $!";
    my (@lines) = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub get_mechanism_species {
    my ($NMVOC, $run) = @_;

    my $mechanism_species;
    if ($NMVOC eq "Pentane") {
        if ($run =~ /MCM|CRI|CB/) {
            $mechanism_species = "NC5H12";
        } elsif ($run =~ /MOZART/) {
            $mechanism_species = "BIGALK";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "HC5";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } elsif ($NMVOC eq "Toluene") {
        if ($run =~ /MCM|CRI|MOZART|CB/) {
            $mechanism_species = "TOLUENE";
        } elsif ($run =~ /RADM|RACM/) {
            $mechanism_species = "TOL";
        } else {
            print "No mechanism species found for $NMVOC\n";
        }
    } else {
        print "No $NMVOC data\n";
    }
    return $mechanism_species;
}

sub get_all_species {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $lines = <$in>;
    close $in;
    my @lines  = split /\n/, $lines;
    my @species;
    foreach my $line (@lines) {
        next unless ($line =~ /IGNORE/);
        $line =~ s/\s=\sIGNORE.*$//;
        push @species, $line;
    }
    return \@species;
}
