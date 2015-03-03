#! /usr/bin/env perl
# ROR production budget in IC5H12 degradation, as this has highest emissions of alkanes with C >= 4
# Version 0: Jane Coates 2/3/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my @mechanisms = qw( CBM-IV CB05 );
my %data;

my $mecca = MECCA->new("$base/CB05_tagged/boxmodel");
my $NTIME = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$NTIME-2);

foreach my $mechanism (@mechanisms) {
    my $boxmodel = "$base/${mechanism}_tagged/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${mechanism}_tagged/gas.eqn";
    my $kpp = KPP->new($eqn_file); 
    my $producers = $kpp->producing("ROR_IC5H12");
    print "No producers found in $mechanism\n" if (@$producers == 0);
    
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        next unless ($reaction_string =~ /ROR_IC5H12 =/);
        my $reaction_number = $kpp->reaction_number($reaction);
        #print "$reaction : $reaction_number\n";
        my $rate = $mecca->rate($reaction_number);
        #my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        $data{$mechanism} = $rate(1:$NTIME-2);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` data = data.frame(Time) `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->set('rate', [ map { $_ } $data{$mechanism}->dog ]);
    $R->run(q` data[mechanism] = rate `);
}
$R->run(q` data = gather(data, Mechanism, Rate, -Time) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
);

$R->run(q` CairoPDF(file = "ROR_reaction_Rates.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
} 
