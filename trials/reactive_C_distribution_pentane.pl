#! /usr/bin/env perl
# Plot distribution of reactive C during pentane degradation in each mechanism
# Version 0: Jane Coates 16/3/2015
#
use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$NTIME-2);

#my @mechanisms = qw( MOZART-4 );
my @mechanisms = qw( MCMv3.2 MCMv3.1 CRIv2 MOZART-4 RADM2 RACM RACM2 CBM-IV CB05 );
my (%n_carbon, %box_data, %area_data);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn);
    my $carbon_file = "$base/${mechanism}_tagged/carbons.txt";
    $n_carbon{$mechanism} = get_carbons($mechanism, $carbon_file);
    my @VOCs = qw( Pentane Toluene );
    foreach my $NMVOC (@VOCs) {
        my $parent = get_mechanism_species($NMVOC, $mechanism);
        ($box_data{$mechanism}{$NMVOC}) = get_box_data($kpp, $mecca, $mechanism, $n_carbon{$mechanism}, $parent);
        #($area_data{$mechanism}{$NMVOC}) = get_area_data($kpp, $mecca, $mechanism, $n_carbon{$mechanism}, $parent);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->set('Times', [ map { $_ } $times->dog ]);
$R->run(q` box.data = data.frame() `,
        q` area.data = data.frame() `,
        q` data = data.frame() `,
);
foreach my $mechanism (sort keys %box_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $VOC (sort keys %{$box_data{$mechanism}}) {
        $R->set('voc', $VOC);
        $R->set('loss.rate', [ map { $_ } $box_data{$mechanism}{$VOC}->dog ]);
        $R->run(q` pre[voc] = loss.rate `);
    }
    $R->set('mechanism', $mechanism);
    $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
            q` pre = gather(pre, VOC, Loss.Rate, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` my.colours = c(  "CB05" = "#0352cb", "CBM-IV" = "#ef6638", "CRIv2" = "#b569b3", "MCMv3.1" = "#000000", "MCMv3.2" = "#dc3522", "MOZART-4" = "#cc9900", "RACM" = "#6c254f", "RACM2" = "#4682b4", "RADM2" = "#035c28") `);
$R->run(q` data$Mechanism = factor(data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Loss.Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + geom_point() `,
        q` plot = plot + facet_wrap( ~ VOC) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0.2)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 3e6)) `,
        q` plot = plot + ylab("Net Carbon Loss Rate (molecules cm-3 s-1)") `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5)) `,
        q` plot = plot + theme(legend.justification = c(1, 0)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.position = c(1, 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "net_reactive_carbon_loss_pentane_toluene.pdf", width = 7.5, height = 5.3) `,
        q` print(plot) `,
        q` dev.off() `,
);

#foreach my $mechanism (sort keys %box_data) {
#    foreach my $VOC (sort keys %{$box_data{$mechanism}}) {
#        $R->run(q` pre = data.frame(Time) `);
#        foreach my $C (sort keys %{$box_data{$mechanism}{$VOC}}) {
#            $R->set('c.number', $C);
#            $R->set('reactive.c', [ map { $_ } $box_data{$mechanism}{$VOC}{$C}->dog ]);
#            $R->run(q` pre[c.number] = reactive.c `);
#        }
#        $R->set('mechanism', $mechanism);
#        $R->set('voc', $VOC);
#        $R->run(q` if("C2.4" %in% colnames(pre)) { pre$C2 = pre$C2 + pre$C2.4 ; pre$C2.4 = NULL }`,
#                q` if("C2.9" %in% colnames(pre)) { pre$C3 = pre$C3 + pre$C2.9 ; pre$C2.9 = NULL }`,
#                q` if("C3.5" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.5 ; pre$C3.5 = NULL }`,
#                q` if("C3.6" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.6 ; pre$C3.6 = NULL }`,
#                q` if("C3.9" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.9 ; pre$C3.9 = NULL }`,
#                q` if("C4.2" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C4.2 ; pre$C4.2 = NULL }`,
#                q` if("C4.5" %in% colnames(pre)) { pre$C5 = pre$C5 + pre$C4.5 ; pre$C4.5 = NULL }`,
#                q` if("C4.8" %in% colnames(pre)) { pre$C5 = pre$C4.8 ; pre$C4.8 = NULL }`,
#                q` if("C5.6" %in% colnames(pre)) { pre$C6 = pre$C5.6 ; pre$C5.6 = NULL }`,
#                q` if("C6.6" %in% colnames(pre)) { pre$C7 = pre$C6.6 ; pre$C6.6 = NULL }`,
#                q` if("C7.1" %in% colnames(pre)) { pre$C7 = pre$C7 + pre$C7.1 ; pre$C7.1 = NULL }`,
#                q` if("C7.75" %in% colnames(pre)) { pre$C8 = pre$C7.75 ; pre$C7.75 = NULL }`,
#        );
#        $R->run(q` pre$Mechanism = rep(mechanism, length(Time)) `,
#                q` pre$VOC = rep(voc, length(Time)) `,
#                q` pre = gather(pre, C.Number, Reactive.C, -Time, -Mechanism, -VOC) `,
#                q` box.data = rbind(box.data, pre) `,
#        );
#    }
#}
##my $p = $R->run(q` print(data) `);
##print $p, "\n";
#
#$R->run(q` my.colours = c("C8" = "#6db875", "C7" = "#0c3f78", "C6" = "#b569b3", "C5" = "#2b9eb3", "C4" = "#ef6638", "C3" = "#0e5628", "C2" = "#f9c500", "C1" = "#6c254f") `);
#$R->run(q` plot = ggplot(box.data, aes(x = Mechanism, y = Reactive.C, fill = C.Number)) `,
#        q` plot = plot + geom_bar(stat = "identity") `,
#        q` plot = plot + facet_grid(Time ~ VOC) `,
#        q` plot = plot + coord_flip() `,
#        q` plot = plot + scale_fill_manual(values = my.colours) `,
#        q` plot = plot + theme_bw() `,
#        q` plot = plot + scale_x_discrete(limits = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")), expand = c(0, 0)) `,
#        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
#        q` plot = plot + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) `,
#        q` plot = plot + theme(legend.margin = unit(0, "lines")) `,
#        q` plot = plot + theme(axis.title.y = element_blank()) `,
#        q` plot = plot + theme(strip.text.y = element_text(face = "bold", angle = 0)) `,
#        q` plot = plot + theme(strip.text.x = element_text(face = "bold")) `,
#        q` plot = plot + theme(strip.background = element_blank()) `,
#        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
#        #q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
#        q` plot = plot + theme(panel.grid = element_blank()) `,
#        q` plot = plot + theme(panel.margin.x = unit(5, "mm")) `,
#        q` plot = plot + theme(legend.position = "bottom") `,
#        q` plot = plot + theme(legend.key = element_blank()) `,
#        q` plot = plot + theme(legend.title = element_blank()) `,
#        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
#);
#
#$R->run(q` CairoPDF(file = "Reactive_C_pentane_box.pdf", width = 8.7, height = 10) `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

#foreach my $mechanism (sort keys %area_data) {
#    foreach my $VOC (sort keys %{$area_data{$mechanism}}) {
#        $R->run(q` pre = data.frame(Times) `);
#        foreach my $C (sort keys %{$area_data{$mechanism}{$VOC}}) {
#            $R->set('c.number', $C);
#            $R->set('reactive.c', [ map { $_ } $area_data{$mechanism}{$VOC}{$C}->dog ]);
#            $R->run(q` pre[c.number] = reactive.c `);
#        }
#        $R->set('mechanism', $mechanism);
#        $R->set('voc', $VOC);
#        $R->run(q` if("C2.4" %in% colnames(pre)) { pre$C2 = pre$C2 + pre$C2.4 ; pre$C2.4 = NULL }`,
#                q` if("C2.9" %in% colnames(pre)) { pre$C3 = pre$C3 + pre$C2.9 ; pre$C2.9 = NULL }`,
#                q` if("C3.5" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.5 ; pre$C3.5 = NULL }`,
#                q` if("C3.6" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.6 ; pre$C3.6 = NULL }`,
#                q` if("C3.9" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C3.9 ; pre$C3.9 = NULL }`,
#                q` if("C4.2" %in% colnames(pre)) { pre$C4 = pre$C4 + pre$C4.2 ; pre$C4.2 = NULL }`,
#                q` if("C4.5" %in% colnames(pre)) { pre$C5 = pre$C5 + pre$C4.5 ; pre$C4.5 = NULL }`,
#                q` if("C4.8" %in% colnames(pre)) { pre$C5 = pre$C4.8 ; pre$C4.8 = NULL }`,
#                q` if("C5.6" %in% colnames(pre)) { pre$C6 = pre$C5.6 ; pre$C5.6 = NULL }`,
#                q` if("C6.6" %in% colnames(pre)) { pre$C7 = pre$C6.6 ; pre$C6.6 = NULL }`,
#                q` if("C7.1" %in% colnames(pre)) { pre$C7 = pre$C7 + pre$C7.1 ; pre$C7.1 = NULL }`,
#                q` if("C7.75" %in% colnames(pre)) { pre$C8 = pre$C7.75 ; pre$C7.75 = NULL }`,
#        );
#        $R->run(q` pre$Mechanism = rep(mechanism, length(Times)) `,
#                q` pre$VOC = rep(voc, length(Times)) `,
#                q` pre = gather(pre, C.Number, Reactive.C, -Times, -Mechanism, -VOC) `,
#                q` area.data = rbind(area.data, pre) `,
#        );
#    }
#}
##my $p = $R->run(q` print(data) `);
##print $p, "\n";
#
#$R->run(q` area.data$Mechanism = factor(area.data$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "MOZART-4", "RADM2", "RACM", "RACM2", "CBM-IV", "CB05")) `);
#$R->run(q` plot = ggplot(area.data, aes(x = Times, y = Reactive.C, fill = C.Number)) `,
#        q` plot = plot + geom_area(position = "stack") `,
#        q` plot = plot + facet_grid(Mechanism ~ VOC) `,
#        q` plot = plot + scale_fill_manual(values = my.colours) `,
#        q` plot = plot + theme_bw() `,
#        q` plot = plot + scale_x_continuous(expand = c(0, 0)) `,
#        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
#        q` plot = plot + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) `,
#        q` plot = plot + theme(legend.margin = unit(0, "lines")) `,
#        q` plot = plot + theme(axis.title.y = element_blank()) `,
#        q` plot = plot + theme(strip.text.y = element_text(face = "bold", angle = 0)) `,
#        q` plot = plot + theme(strip.text.x = element_text(face = "bold")) `,
#        q` plot = plot + theme(strip.background = element_blank()) `,
#        q` plot = plot + theme(axis.title.x = element_text(face = "bold")) `,
#        #q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
#        q` plot = plot + theme(panel.grid = element_blank()) `,
#        q` plot = plot + theme(panel.margin.x = unit(5, "mm")) `,
#        q` plot = plot + theme(legend.position = "bottom") `,
#        q` plot = plot + theme(legend.key = element_blank()) `,
#        q` plot = plot + theme(legend.title = element_blank()) `,
#        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
#);
#
#$R->run(q` CairoPDF(file = "Reactive_C_pentane_area.pdf", width = 8.7, height = 10) `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

$R->stop();

sub get_area_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my $all_reactions = $kpp->all_reactions();
    my %reactive_C;
    
    for (0..$#$all_reactions) { #get rates for all producing reactions
        my $reaction = $all_reactions->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent =~ $VOC);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        foreach my $product (@{$kpp->products($reaction)}) {
            (my $product_nontagged = $product) =~ s/_(.*?)\b//g;
            my $C_number = get_species_carbon($product_nontagged, $carbons);
            next if ($C_number == 0);
            my @reactions = ();
            push @reactions, $reaction;
            my $yield = $kpp->yield_of($product, \@reactions);
            $reactive_C{"C$C_number"} += ${$yield}[0] * $C_number * $rate(1:$NTIME-2);
        }
    }

    foreach my $C (sort keys %reactive_C) {
        if ($VOC =~ /TOL/) {
            if ($mechanism eq "RACM2") {
                $reactive_C{$C} = $reactive_C{$C} * 0.868;
            } elsif ($mechanism =~ /RA/) {
                $reactive_C{$C} = $reactive_C{$C} * 0.667;
            } elsif ($mechanism =~ /MO/) {
                $reactive_C{$C} = $reactive_C{$C} * 0.478;
            }
        } else { #pentane
            if ($mechanism =~ /RA/){
                $reactive_C{$C} = $reactive_C{$C} * 0.264;
            } elsif ($mechanism =~ /MO/) {
                $reactive_C{$C} = $reactive_C{$C} * 0.146;
            }
        }
    } 
    return \%reactive_C;
}

sub get_box_data {
    my ($kpp, $mecca, $mechanism, $carbons, $VOC) = @_;
    my $all_reactions = $kpp->all_reactions();
    my %reactive_C;
    my $net_C;

    for (0..$#$all_reactions) { #get rates for all producing reactions
        my $reaction = $all_reactions->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent =~ $VOC);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        my $reactants = $kpp->reactants($reaction);
        my $reactant_sum = join ' + ', @$reactants;
        #$reactant_sum =~ s/_(.*?)\b//g;
        my ($org_reactant, $reactant_C);
        foreach (@$reactants) {
            next unless ($_ =~ /_/ or $_ eq $parent);
            ($org_reactant = $_) =~ s/_(.*?)\b//g;
            $reactant_C = get_species_carbon($org_reactant, $carbons);
        }
        my $product_C = 0;
        foreach my $product (@{$kpp->products($reaction)}) {
            (my $product_nontagged = $product) =~ s/_(.*?)\b//g;
            my $C_number = get_species_carbon($product_nontagged, $carbons);
            next if ($C_number == 0);
            my @reactions = ();
            push @reactions, $reaction;
            my $yield = $kpp->yield_of($product, \@reactions);
            $product_C += ${$yield}[0] * $C_number;
            $reactive_C{"C$C_number"} += ${$yield}[0] * $C_number * $rate(1:$NTIME-2);
        }
        my $net = $product_C - $reactant_C;
        next if ($net == 0);
        $net_C += $net * $rate(1:$NTIME-2);
    }
    $net_C = $net_C->reshape($N_PER_DAY, $N_DAYS);
    $net_C = $net_C->sumover;
    $net_C = $net_C(0:13:2);
    if ($VOC =~ /TOL/) {
        if ($mechanism eq "RACM2") {
            $net_C = $net_C * 0.868;
        } elsif ($mechanism =~ /RA/) {
            $net_C = $net_C * 0.667;
        } elsif ($mechanism =~ /MO/) {
            $net_C = $net_C * 0.478;
        }
    } else { #pentane
        if ($mechanism =~ /RA/){
            $net_C = $net_C * 0.264;
        } elsif ($mechanism =~ /MO/) {
            $net_C = $net_C * 0.146;
        }
    }
    return $net_C;
    
#    foreach my $C (sort keys %reactive_C) {
#        if ($VOC =~ /TOL/) {
#            if ($mechanism eq "RACM2") {
#                $reactive_C{$C} = $reactive_C{$C} * 0.868;
#            } elsif ($mechanism =~ /RA/) {
#                $reactive_C{$C} = $reactive_C{$C} * 0.667;
#            } elsif ($mechanism =~ /MO/) {
#                $reactive_C{$C} = $reactive_C{$C} * 0.478;
#            }
#        } else { #pentane
#            if ($mechanism =~ /RA/){
#                $reactive_C{$C} = $reactive_C{$C} * 0.264;
#            } elsif ($mechanism =~ /MO/) {
#                $reactive_C{$C} = $reactive_C{$C} * 0.146;
#            }
#        }
#        my $reshape = $reactive_C{$C}->reshape($N_PER_DAY, $N_DAYS);
#        my $integrate = $reshape->sumover;
#        $reactive_C{$C} = $integrate(0:13:2); 
#    } 
    #return \%reactive_C;
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

sub get_species_carbon {
    my ($species, $carbons) = @_;
    my %carbons = %$carbons;
    my $carbon;
    if (exists $carbons{$species}) {
        $carbon = $carbons{$species};
    } else {
        $carbon = 0;
        print "No C found for species: $species\n";
    }
    return $carbon;
}
