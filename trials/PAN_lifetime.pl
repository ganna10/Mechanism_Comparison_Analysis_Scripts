#!/usr/bin/env perl
# PAN lifetime
# Version 0: Jane Coates 10/7/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use PDL::Bad;
use Statistics::R;

my $model_run = "/local/home/coates/MECCA/MCMv3.2_tagged";
my $boxmodel = "$model_run/boxmodel";
my $mecca = MECCA->new($boxmodel);
my $eqn = "$model_run/gas.eqn";
my $kpp = KPP->new($eqn);
my $ntime = $mecca->time->nelem;

my $reactant = 'PAN_C5H8';
my $cair = $mecca->cair;
my $consumers = $kpp->consuming($reactant);
die "No reactions found for $reactant\n" if (@$consumers == 0) ;

my $lifetime;
foreach my $reaction (@$consumers) {
    my $reactants = $kpp->reactants($reaction);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if (isbad($rconst(-1)));
    my $reactivity;
    if (defined $other_reactant) {
        my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
        $reactivity = $rconst * $other_reactant_conc;
    } else {
        $reactivity = $rconst * $cair * $mecca->tracer($reactant);
    }
    $lifetime += $reactivity(1:$ntime-2)**(-1); 
}

#$lifetime /= 86400;
print $lifetime, "\n";
