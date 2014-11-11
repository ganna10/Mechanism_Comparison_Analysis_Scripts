#! /usr/bin/env perl
# compare Ox production from cresol degradation in MCM v3.2, CB4 and CB5
# Version 0: Jane Coates 11/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA";
my $mecca = MECCA->new("$base/CB05_tagging/boxmodel");
my $NTIME = $mecca->time->nelem;
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT;
my $N_DAYS = int ($NTIME / $N_PER_DAY);

my @runs = qw( MCM_3.2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) CBM-IV", "(c) CB05" );
my @species = qw( CRESOL CRES CRES );
#my @runs = qw( CB05_tagging );
#my @mechanisms = ( "(c) CB05" );
#my @species = qw( CRES );
my $index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $cresol_file = "$base/$run/cresol.eqn"; 
    my $cresol_kpp = KPP->new($cresol_file);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $RO2_file = "$base/$run/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    my @cresol_no2_reservoirs = map { local $_ = $_; s/TOLUENE/$species[$index]/; $_ } grep { $_ =~ /TOLUENE/ } @no2_reservoirs;
    $families{"Ox_$run"} = [ qw( O3 NO2 O1D O HO2NO2 NO3 N2O5 ), @cresol_no2_reservoirs ]; 
    $weights{"Ox_$run"} = { NO3 => 2, N2O5 => 3 };
    ($plot_data{$mechanisms[$index]}) = get_data($mecca, $kpp, $cresol_kpp, $run, $species[$index]);
    $index++;
}

my $R = Statistics::R->new(); 
$R->run(q` library(ggplot2) `,
        q` library(grid) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(gridExtra) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plots = list() `);

$R->run(q` plotting = function (data, mechanism) {  plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) ;
                                                    plot = plot + geom_bar(stat = "identity") ;
                                                    plot = plot + theme_bw() ;
                                                    plot = plot + ggtitle(mechanism) ;
                                                    plot = plot + theme(axis.title.x = element_blank()) ;
                                                    plot = plot + theme(axis.title.y = element_blank()) ;
                                                    plot = plot + theme(axis.text.x = element_text(size = 55, angle = 45, vjust = 0.5)) ;
                                                    plot = plot + theme(axis.text.y = element_text(size = 50)) ;
                                                    plot = plot + theme(plot.title = element_text(size = 70, face = "bold")) ;
                                                    plot = plot + theme(panel.grid.major = element_blank()) ;
                                                    plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                    plot = plot + theme(axis.ticks.length = unit(0.5, "cm")) ;
                                                    plot = plot + theme(axis.ticks.margin = unit(0.3, "cm")) ;
                                                    plot = plot + theme(legend.title = element_blank()) ;
                                                    plot = plot + theme(legend.key = element_blank()) ;
                                                    plot = plot + theme(legend.text = element_text(size = 30)) ;
                                                    plot = plot + theme(legend.key.size = unit(2.5, "cm")) ;
                                                    plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                    plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                    return(plot) } `);

foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) { 
            next if ($reaction =~ /^C2O3 |MEO2|CH3O2|CH3CO3/);
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = rate `);
        }
    }
    $R->set('mechanism', $run);
    $R->run(q` data$Mechanism = rep(mechanism, length(Time)) `,
            q` data = melt(data, id.vars = c("Time", "Mechanism"), variable.name = "Reaction", value.name = "Rate") `,
            q` plot = plotting(data, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
#my $p = $R->run(q` print(plot) `);
#print "$p\n";
}

$R->run(q` CairoPDF(file = "Cresol_Ox_production.pdf", width = 50, height = 35) `,
        q` multiplot = grid.arrange(arrangeGrob(plots[[1]],
                                                plots[[2]],
                                                plots[[3]],
                                                nrow = 1),
                                    nrow = 1, ncol = 1,
                                    left = textGrob(expression(bold(paste("Production Rates (molecules ", cm^-3, s^-1, ")"))), rot = 90, gp = gpar(fontsize = 85), vjust = 0.5) )`,
        q` print(plots[[1]]) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $cresol_kpp, $run, $cresol) = @_;
    $families{"HO2x"} = [ qw(HO2 HO2NO2) ];
    my @loop = ("Ox_$run");

    my (%production_rates, %consumption_rates, %plot_data);
    foreach my $species (@loop) { #get cresol producers and yields from cresol.eqn file
        my ($producers, $producer_yields, $consumers, $consumer_yields);
        if (exists $families{$species}) {
            $cresol_kpp->family({
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
            });
            $producers = $cresol_kpp->producing($species);
            $producer_yields = $cresol_kpp->effect_on($species, $producers);
            $consumers = $cresol_kpp->consuming($species);
            $consumer_yields = $cresol_kpp->effect_on($species, $consumers);
        } else {
            print "No family found for $species\n";
        }
        print "No producers found for $species\n" if (@$producers == 0);
        print "No consumers found for $species\n" if (@$consumers == 0);

        for (0..$#$producers) { #change CRES tag to TOLUENE tag so that the reaction data can be extracted from gas.eqn file, need to do this otherwise reaction_number doesn't match for mecca
            my $reaction = $producers->[$_];
            next unless ($reaction =~ /_/);
            $reaction =~ s/$cresol/TOLUENE/;
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $products) = split / = /, $reaction_string;
            if ($reactants =~ /XO2/) {
                my $operator = "XO2_CRES";
                my $op_producers = $cresol_kpp->producing($operator);
                my $op_producer_yields = $cresol_kpp->effect_on($operator, $op_producers);
                print "No producers for $operator\n" if (@$op_producers == 0);

                for (0..$#$op_producers) {
                    my $reaction = $op_producers->[$_];
                    $reaction =~ s/$cresol/TOLUENE/;
                    my $reaction_number = $kpp->reaction_number($reaction);
                    my $rate = $mecca->rate($reaction_number) * $op_producer_yields->[$_];
                    next if ($rate->sum == 0);
                    my $new_string = $kpp->reaction_string($reaction);
                    $new_string =~ s/_(.*?)\b//g;
                    my ($op_reactants, $op_products) = split / = /, $new_string;
                    $production_rates{$op_reactants} += $rate(1:$NTIME-2);
                }
            } else {
                $production_rates{$reactants} += $rate(1:$NTIME-2);
            }
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            next unless ($reaction =~ /_/);
            $reaction =~ s/$cresol/TOLUENE/;
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            $reaction_string =~ s/_(.*?)\b//g;
            my ($reactants, $consumts) = split / = /, $reaction_string;
            $consumption_rates{$reactants} += $rate(1:$NTIME-2);
        }
        remove_common_processes(\%production_rates, \%consumption_rates);
    }

    my $prod_max = 3e6;
    foreach my $reaction (keys %production_rates) {
        if ($production_rates{$reaction}->sum < $prod_max) {
            $production_rates{"Production Others"} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
    }
    
    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrated = $reshape->sumover;
        $integrated = $integrated(0:13:2);
        $production_rates{$reaction} = $integrated;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_data = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    
    my @final_sorted_data;
    foreach (@sorted_data) { 
        next if ($_ eq 'Production Others') ;
        push @final_sorted_data, { $_ => $production_rates{$_} };
    } 
    push @final_sorted_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); 

    return \@final_sorted_data;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 

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
