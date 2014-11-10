#! /usr/bin/env perl
# Plot alkane kOH versus carbon number in each mechanism
# Version 0: Jane Coates 10/11/2014

use strict;
use diagnostics;
use Statistics::R;

my @C_numbers = (2..8);
my %k_OH = ("MCMv3.2"   => [ 2.27e-13, 1.03e-12, 2.3e-12, 3.91e-12, 5.4e-12, 6.98e-12, 8.61e-12 ],
            "MCMv3.1"   => [ 2.38e-13, 1.08e-12, 2.38e-12, 3.91e-12, 5.4e-12, 6.98e-12, 8.60e-12 ],
            "CRIv2"     => [ 2.38e-13, 1.08e-12, 2.38e-12, 3.91e-12, 5.4e-12, 6.98e-12, 8.61e-12 ],
            "MOZART-4"  => [ 2.26e-13, 1.04e-12, 999, 3.5e-12, 999, 999, 999 ],
            "RADM2"     => [ 2.58e-13, 2.52e-12, 999, 4.73e-12, 999, 999, 9.95e-12 ],
            "RACM"      => [ 2.42e-13, 2.17e-12, 999, 4.73e-12, 999, 999, 1.07e-11 ],
            "RACM2"     => [ 2.27e-13, 2.17e-12, 999, 4.38e-12, 999, 999, 1.11e-11 ],
            "CBM-IV"    => [ 3.24e-13, 1.22e-12, 3.24e-12, 4.05e-12, 4.86e-12, 5.67e-12, 6.48e-12 ],
            "CB05"      => [ 2.26e-13, 1.22e-12, 3.24e-12, 4.05e-12, 4.86e-12, 5.67e-12, 6.48e-12 ],
);

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('carbons', [@C_numbers]);
$R->run(q` data = data.frame(carbons) `);

foreach my $mechanism (sort keys %k_OH) {
    $R->set('mechanism', $mechanism);
    $R->set('kOH', [@{$k_OH{$mechanism}}]);
    $R->run(q` data[mechanism] = kOH `);
}
$R->run(q` data = melt(data, id.vars = c("carbons"), variable.name = "Mechanism", value.name = "kOH") `,
        q` data = as.data.frame(lapply(data, function(x) { replace(x, x == 999, NA) })) `,
);
my $p = $R->run(q` print(data) `);
print "$p\n";

$R->run(q` plot = ggplot(data, aes(x = carbons, y = kOH, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_point() `,
        #q` plot = plot + geom_line() `,
);

$R->run(q` CairoPDF(file = "alkane_kOH_vs_C_number.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
