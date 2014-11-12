#! /usr/bin/env perl
# compare PPN concentration of C3-C8 alkanes and butene, ethylbenzene between MCM v3.1 and v3.2

use strict;
use diagnostics;
use MECCA;
use Statistics::R;
use PDL;
use PDL::NiceSlice;
use PDL::NetCDF;

my @VOCs = qw( C3H8 NC4H10 NC5H12 IC5H12 NC6H14 NC7H16 NC8H18 BUT1ENE EBENZ );
my %concs;

my $base_dir = "/work/users/jco/MECCA";
my $boxmodel = "$base_dir/MCM_3.2_tagged/boxmodel";
my $mecca = MECCA->new($boxmodel);
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my @time_axis = map { $_ } $times->dog;

#MCM 3.1
my $conc_file_one = "$base_dir/MCM_3.1_tagged_3.2rates/boxmodel/mecca1_tracer.nc";
my $nc_one = PDL::NetCDF->new($conc_file_one);
foreach my $VOC (@VOCs) {
    my $species = "PPN_$VOC";
    ($concs{'MCM 3.1'}{$VOC}) = get_conc($species, $nc_one);
}

#MCM 3.2
my $conc_file_two = "$base_dir/MCM_3.2_tagged/boxmodel/mecca1_tracer.nc";
my $nc_two = PDL::NetCDF->new($conc_file_two);
foreach my $VOC (@VOCs) {
    my $species = "PPN_$VOC";
    ($concs{'MCM 3.2'}{$VOC}) = get_conc($species, $nc_two);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(scales) `,
        q` library(grid) `,
        q` library(Cairo) `,
);

$R->set('time', [@time_axis]);   
$R->set('rep.number', scalar(@VOCs) * scalar(keys %concs)); 

$R->run(q` times = rep(time, rep.number) `,
        q` VOC = {} `,
        q` Mechanism = {} `,
        q` Mixing.ratio = {} `,
);

foreach my $item (sort keys %concs) {
    $R->set('mechanism', $item);
    foreach my $spec (sort keys %{$concs{$item}} ) {
        $R->set('name', $spec) ;
        $R->set('mixing.ratio', [@{$concs{$item}{$spec}}]);
        $R->run(q` VOC = cbind(VOC, rep(name, length(time))) `,
                q` Mechanism = cbind(Mechanism, rep(mechanism, length(time))) `,
                q` Mixing.ratio = cbind(Mixing.ratio, mixing.ratio) `,
        );
    }
} 

$R->run(q` VOC = c(VOC) `,
        q` Mechanism = c(Mechanism) `,
        q` Mixing.ratio = c(Mixing.ratio) `,
);

$R->run(q` plot.lines = function () {list(    ylab("Concentration (pptv)"), 
                                              xlab("Time (days)"), 
                                              ggtitle("PPN Concentration of VOC Degradation in MCM v3.1 and v3.2"), 
                                              geom_line(size = 3), 
                                              scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), 
                                              theme_bw(), 
                                              facet_wrap( ~ VOC),
                                              theme(axis.title.x = element_text(size = 30, face = "bold")), 
                                              theme(axis.title.y = element_text(size = 30, face = "bold")), 
                                              scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 50)), 
                                              theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), legend.title = element_text(size = 25, face = "bold"), legend.key.size = unit(3, "cm"), legend.text = element_text(size = 20), legend.key = element_blank(), plot.title = element_text(size = 40, face = "bold"), strip.text.x = element_text(size = 30, face = "bold")), 
                                              scale_colour_manual(values = my.colours) ) } `); 

$R->run(q` my.colours = c("black", "red") `);
$R->run(q` data = data.frame(times, VOC, Mechanism, Mixing.ratio) `,
        q` plot = ggplot(data = data, aes(x = times, y = Mixing.ratio, colour = Mechanism)) `,
        q` plot = plot + plot.lines() `,
        q` CairoPNG(file = "PPN_comparison.png", width = 2000, height = 1500) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

sub get_conc {
    my ($species, $nc) = @_;

    my $conc = $nc->get($species)->squeeze;
    $conc = $conc(1:$ntime-2) * 1e12;
    my @conc_array = map { $_ } $conc->dog;
    return \@conc_array;
}
