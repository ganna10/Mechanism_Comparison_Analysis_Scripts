#! /usr/bin/perl
# Calculate RMSE of each mechanisms's first day TOPP to MCM v3.2 first day TOPP
# Version 0: Jane Coates 18/12/2014

use strict;
use diagnostics;
use Statistics::R;

my %TOPP;
my $base = "/local/home/coates/Documents/Analysis/2014_Mechanism_comparison_paper/TOPP_plots";
opendir DIR, $base or die "Can't open $base : $!";
my @daily_TOPP_files = grep { $_ =~ /TOPP_values/ } readdir DIR;
closedir DIR;

foreach my $file (@daily_TOPP_files) {
    my $path = "$base/$file";
    my @lines = split /\n/, read_file($path);
    (my $mechanism = $file) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $VOC = get_chemical_name($VOC);
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $TOPP{$mechanism}{$VOC} = $TOPPs[0];
    }
}

my $R = Statistics::R->new();
$R->run(q` library(tidyr) `);

$R->set('VOCs', [sort keys %{$TOPP{"MCMv3.2"}}]);
$R->run(q` data = data.frame(VOC = VOCs) `); 
foreach my $mechanism (sort keys %TOPP) {
    $R->set('mechanism', $mechanism);
    $R->run(q` topp = c() `);
    foreach my $VOC (sort keys %{$TOPP{"MCMv3.2"}}) { 
        if (defined $TOPP{$mechanism}{$VOC}) {
            $R->set('value', $TOPP{$mechanism}{$VOC});
        } else {
            $R->set('value', 999);
        }
        $R->run(q` topp = c(topp, value) `);
    }
    $R->run(q` data[mechanism] = topp `);
}
$R->run(q` data = as.data.frame(lapply(data, function(x) { replace(x, x== 999, NA) })) `);
$R->run(q` rmse = function (error) { sqrt(mean(error^2)) } `);
$R->run(q` CB05.error = data$MCMv3.2 - data$CB05 `,
        q` CB05.rmse = rmse(CB05.error) `,
        q` CBM.IV.error = data$MCMv3.2 - data$CBM.IV `,
        q` CBM.IV.rmse = rmse(CBM.IV.error) `,
        q` CRIv2.error = data$MCMv3.2 - data$CRIv2 `,
        q` CRIv2.rmse = rmse(CRIv2.error) `,
        q` MCMv3.1.error = data$MCMv3.2 - data$MCMv3.1 `,
        q` MCMv3.1.rmse = rmse(MCMv3.1.error) `,
        q` MOZART.4.error = data$MCMv3.2 - data$MOZART.4 `,
        q` MOZART.4.rmse = rmse(MOZART.4.error) `,
        q` RACM2.error = data$MCMv3.2 - data$RACM2 `,
        q` RACM2.rmse = rmse(RACM2.error) `,
        q` RACM.error = data$MCMv3.2 - data$RACM `,
        q` RACM.rmse = rmse(RACM.error) `,
        q` RADM2.error = data$MCMv3.2 - data$RADM2 `,
        q` RADM2.rmse = rmse(RADM2.error) `,
);
#my $p = $R->run(q` print (CBM.IV.error) `);
#print $p, "\n";

$R->run(q` out.file = file("RMSE_values.txt") `,
        q` writeLines(c("CB05", CB05.rmse, "", "CBM-IV", CBM.IV.rmse, "", "CRIv2", CRIv2.rmse, "", "MCMv3.1", MCMv3.1.rmse, "", "MOZART-4", MOZART.4.rmse, "", "RACM2", RACM2.rmse, "", "RACM", RACM.rmse, "", "RADM2", RADM2.rmse), out.file) `,
        q` close(out.file) `,
);
$R->stop();

sub read_file {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    close $in;
    return $data;
}

sub get_chemical_name {
    my ($VOC) = @_;
    my $chemical_species;
    if ($VOC eq 'C2H6' or $VOC eq 'ETH') {
        $chemical_species = 'Ethane ';
    } elsif ($VOC eq 'C3H8' or $VOC eq 'HC3') {
        $chemical_species = 'Propane ';
    } elsif ($VOC eq 'NC4H10') {
        $chemical_species = 'Butane ';
    } elsif ($VOC eq 'IC4H10') {
        $chemical_species = '2-Methylpropane ';
    } elsif ($VOC eq 'NC5H12' or $VOC eq 'BIGALK' or $VOC eq 'HC5') {
        $chemical_species = 'Pentane ';
    } elsif ($VOC eq 'IC5H12') {
        $chemical_species = '2-Methylbutane ';
    } elsif ($VOC eq 'NC6H14') {
        $chemical_species = 'Hexane ';
    } elsif ($VOC eq 'NC7H16') {
        $chemical_species = "Heptane ";
    } elsif ($VOC eq 'NC8H18' or $VOC eq 'HC8') {
        $chemical_species = "Octane ";
    } elsif ($VOC eq 'C2H4' or $VOC eq 'OL2' or $VOC eq 'ETE') {
        $chemical_species = 'Ethene ';
    } elsif ($VOC eq 'C3H6' or $VOC eq 'OLT') {
        $chemical_species = 'Propene ';
    } elsif ($VOC eq 'BUT1ENE' or $VOC eq 'BIGENE') {
        $chemical_species = "Butene ";
    } elsif ($VOC eq 'MEPROPENE' or $VOC eq 'OLI') {
        $chemical_species = '2-Methylpropene ';
    } elsif ($VOC eq 'C5H8' or $VOC eq 'ISO' or $VOC eq 'ISOP') {
        $chemical_species = "Isoprene ";
    } elsif ($VOC eq 'BEN' or $VOC eq 'BENZENE') {
        $chemical_species = "Benzene ";
    } elsif ($VOC eq 'TOLUENE' or $VOC eq 'TOL') {
        $chemical_species = 'Toluene ';
    } elsif ($VOC eq 'MXYL' or $VOC eq 'XYM' or $VOC eq 'XYL') {
        $chemical_species = "m-Xylene ";
    } elsif ($VOC eq 'OXYL' or $VOC eq 'XYO') {
        $chemical_species = 'o-Xylene ';
    } elsif ($VOC eq 'PXYL' or $VOC eq 'XYP') {
        $chemical_species = "p-Xylene ";
    } elsif ($VOC eq 'EBENZ') {
        $chemical_species = "Ethylbenzene ";
    } else {
        print "No chemical species found for $VOC\n";
    }
}
