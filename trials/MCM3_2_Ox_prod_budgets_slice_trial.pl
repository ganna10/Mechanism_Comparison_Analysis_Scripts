#!/usr/bin/perl -w
# calculate Ox production from tagged and non-tagged model runs and plot as stacked bar plot
# Version 0: Jane Coates 31/05/2014
# Version 1: Jane Coates 11/08/2014 including consumption calculation and plotting net Ox production

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my %families = (    'tagged_HO2x' => [ qw( HO2 HO2NO2 ) ],
                    'HO2x' => [ qw( HO2 HO2NO2 ) ] );
my (%weights, %production_reaction_rates) ;

####tagged run
my $tagged_run = "/work/users/jco/MECCA/MCM_3.2_tagged/boxmodel";
my $tagged_mecca = MECCA->new($tagged_run); 
my $tagged_eqnfile = "/work/users/jco/MECCA/MCM_3.2_tagged/gas.eqn";
my $tagged_kpp = KPP->new($tagged_eqnfile); 
my $tagged_ro2file = "/work/users/jco/MECCA/MCM_3.2_tagged/RO2_species.txt";
my @tagged_no2_reservoirs = get_no2_reservoirs($tagged_kpp, $tagged_ro2file);
$families{'tagged_Ox'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @tagged_no2_reservoirs ];
$weights{'tagged_Ox'} = { NO3 => 2, N2O5 => 3 };
my $ntime = $tagged_mecca->time->nelem;

foreach my $species (qw( tagged_Ox tagged_HO2x )) { #get all production and consumption rates
    my ($producers, $producer_yields);
    if (exists $families{$species}) { #get family reaction numbers and yields
        $tagged_kpp->family({ 
                        name    => $species,
                        members => $families{$species},
                        weights => $weights{$species},
        });
        $producers = $tagged_kpp->producing($species);
        $producer_yields = $tagged_kpp->effect_on($species, $producers);  
    } else {
        print "No family found for $species\n";
    }
    die "No producers found for $species\n" if (@$producers == 0);#check that species reactions are found
    
    for (0..$#$producers) { #get rates for all producing reactions
        my $reaction = $producers->[$_];
		my $reaction_number = $tagged_kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $tagged_mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        my $string;
        if (defined $parent) { # for tagged reactions
            $string = $tagged_kpp->reaction_string($reaction);
            $string =~ s/_$parent//g; #removing tag from reaction strings
            $string = $parent; # in order to merge all production rates from all parent species reactions into 1 pdl
        } else { # for non-tagged reactions
            $string = $tagged_kpp->reaction_string($reaction);
        }
        $production_reaction_rates{$species}{$string} += $rate(1:$ntime-2); #attribute rates to each parent tagged species and the non-tagged reactions
    }
} 

my $tagged_ho2x_total_production = zeroes(PDL::float, $ntime-2);
$tagged_ho2x_total_production += $production_reaction_rates{'tagged_HO2x'}{$_} for (keys %{ $production_reaction_rates{'tagged_HO2x'} });

foreach my $reaction ( keys %{ $production_reaction_rates{'tagged_HO2x'} }) {
    $production_reaction_rates{'tagged_Ox'}{$reaction} += $production_reaction_rates{'tagged_Ox'}{'HO2 + NO = NO2 + OH'} * $production_reaction_rates{'tagged_HO2x'}{$reaction} / $tagged_ho2x_total_production;
}
delete $production_reaction_rates{'tagged_Ox'}{'HO2 + NO = NO2 + OH'};

my $tagged_prod_others_max = 1e8;
my ($tagged_plot_data) = get_plot_data($tagged_prod_others_max, \%{$production_reaction_rates{'tagged_Ox'}});

######non-tagged run
my $run = "/work/users/jco/MECCA/MCM_3.2_no_tagging/boxmodel";
my $mecca = MECCA->new($run); 
my $eqnfile = "/work/users/jco/MECCA/MCM_3.2_no_tagging/gas.eqn";
my $kpp = KPP->new($eqnfile); 
my $ro2file = "/work/users/jco/MECCA/MCM_3.2_no_tagging/RO2_species.txt";
my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);
$families{'Ox'} = [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
$weights{'Ox'} = { NO3 => 2, N2O5 => 3 };

foreach my $species (qw( Ox HO2x )) { #get all production and consumption rates
    my ($producers, $producer_yields);
    if (exists $families{$species}) { #get family reaction numbers and yields
        $kpp->family({ 
                        name    => $species,
                        members => $families{$species},
                        weights => $weights{$species},
        });
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers);  
    } else {
        print "No family found for $species\n";
    }
    die "No producers found for $species\n" if (@$producers == 0);#check that species reactions are found
    
    for (0..$#$producers) { #get rates for each producing reactions
        my $reaction = $producers->[$_];
		my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my $string = $kpp->reaction_string($reaction);
        $production_reaction_rates{$species}{$string} += $rate(1:$ntime-2);
    }
} 

my $ho2x_total_production = zeroes(PDL::float, $ntime-2);
$ho2x_total_production += $production_reaction_rates{'HO2x'}{$_} for (keys %{ $production_reaction_rates{'HO2x'} });

foreach my $reaction ( keys %{ $production_reaction_rates{'HO2x'} }) {
    $production_reaction_rates{'Ox'}{$reaction} += $production_reaction_rates{'Ox'}{'HO2 + NO = NO2 + OH'} * $production_reaction_rates{'HO2x'}{$reaction} / $ho2x_total_production;
}
delete $production_reaction_rates{'Ox'}{'HO2 + NO = NO2 + OH'};

my $prod_others_max = 1e8;
my ($plot_data) = get_plot_data($prod_others_max, \%{$production_reaction_rates{'Ox'}});

#Create x-axis for plot in hours
my $times = $tagged_mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog;
#map to day and night
my @time_blocks;
foreach my $time (@time_axis) {
    if ($time <= 12) {
        push @time_blocks, "Day 1";
    } elsif ($time > 12 and $time <= 24) {
        push @time_blocks, "Night 1";
    } elsif ($time > 24 and $time <= 36) {
        push @time_blocks, "Day 2";
    } elsif ($time > 36 and $time <= 48) {
        push @time_blocks, "Night 2";
    } elsif ($time > 48 and $time <= 60) {
        push @time_blocks, "Day 3",
    } elsif ($time > 60 and $time <= 72) {
        push @time_blocks, "Night 3";
    } elsif ($time > 72 and $time <= 84) {
        push @time_blocks, "Day 4";
    } elsif ($time > 84 and $time <= 96) {
        push @time_blocks, "Night 4";
    } elsif ($time > 96 and $time <= 108) {
        push @time_blocks, "Day 5";
    } elsif ($time > 108 and $time <= 120) {
        push @time_blocks, "Night 5";
    } elsif ($time > 120 and $time <= 132) {
        push @time_blocks, "Day 6";
    } elsif ($time > 132 and $time <= 144) {
        push @time_blocks, "Night 6";
    } elsif ($time > 144 and $time <= 156) {
        push @time_blocks, "Day 7";
    } else {
        push @time_blocks, "Night 7";
    }
}

my ($reaction_plot_data) = plot(\@time_blocks, $tagged_plot_data, $plot_data);

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open FILE, $file or die $!; 
    my @ro2;
    for (<FILE>) {
        push @ro2, split /\s+/, $_; 
    }
    close FILE;
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

sub get_plot_data {
    my ($prod_others_max, $production_rates) = @_;
    my %production_rates = %$production_rates;

    my (@production_others, %production_plot_rates, @sorted_plot_data, @plot_data);

    foreach my $item (keys %production_rates) {#sort production
        if ($production_rates{$item}->sum < $prod_others_max) { #get production others
            push @production_others, $production_rates{$item};
            my $prod_other_rates = cat(@production_others);
            $prod_other_rates = $prod_other_rates->xchg(0,1)->sumover;
            $production_rates{'Production Others'} = $prod_other_rates;
            delete $production_rates{$item};
        }
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;

    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    } 
    push @sorted_plot_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); #add Production Others to the beginning 

    foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
        foreach my $item (keys %{$ref}) {
            my $data = $ref->{$item};
            my $day1 = $data(0:35);
            my $night1 = $data(36:71);
            my $day2 = $data(72:107);
            my $night2 = $data(108:143);
            my $day3 = $data(144:179);
            my $night3 = $data(180:215);
            my $day4 = $data(216:251);
            my $night4 = $data(252:287);
            my $day5 = $data(288:323);
            my $night5 = $data(324:359);
            my $day6 = $data(360:395);
            my $night6 = $data(396:431);
            my $day7 = $data(432:467);
            my $night7 = $data(468:503);
            my @vals = ($day1->sum, $night1->sum, $day2->sum, $night2->sum, $day3->sum, $night3->sum, $day4->sum, $night4->sum, $day5->sum, $night5->sum, $day6->sum, $night6->sum, $day7->sum, $night7->sum);
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@vals};
        }
    } 

    return \@plot_data;
}

sub plot { #create dataframe and then create plot
    my ($time, $tagged_plot_data, $plot_data) = @_;
    
    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(Cairo) `);
    
    $R->run(q` time = c("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4 ", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7") `);
    $R->run(q` plot.function = function(data, colours, names, plot.title, legend.title, y.max, y.breaks) {  
                                                plot.data = melt( data = data, id = names(data)[1], measured = names(data)[-1] ) ;
                                                colnames(plot.data) = c("time", "NMVOC", "rate") ; 
                                                NMVOC.levels = (levels(factor(plot.data$NMVOC))) ;
                                                plot.data$NMVOC = ordered(plot.data$NMVOC, levels = NMVOC.levels) ;
                                                scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } ;
                                                plot = ggplot(plot.data, aes(x = time, y = rate, fill = NMVOC)) ;
                                                plot = plot + geom_bar(stat = "identity", width = 0.6) ;
                                                plot = plot + theme_bw() ;
                                                plot = plot + theme(axis.text.x = element_text(size = 75)) ;
                                                plot = plot + theme(axis.text.y = element_text(size = 90)) ;
                                                plot = plot + theme(axis.title.y = element_blank()) ;
                                                plot = plot + theme(axis.title.x = element_blank()) ;
                                                plot = plot + theme(legend.key.size = unit(6.5, "cm")) ;
                                                plot = plot + theme(legend.text = element_text(size = 100)) ;
                                                plot = plot + theme(legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95)) ; 
                                                plot = plot + theme(legend.title = element_text(size = 120, face = "bold")) ;
                                                plot = plot + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) ;
                                                plot = plot + theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) ;
                                                plot = plot + theme(axis.ticks.length = unit(1.5, "cm"));
                                                plot = plot + theme(axis.ticks.margin = unit(.9, "cm"));
                                                plot = plot + theme(axis.ticks = element_line(size = 5)) ;
                                                plot = plot + scale_y_continuous(limits=c(0, y.max), breaks=seq(0, y.max, y.breaks), label = scientific_10);
                                                plot = plot + theme(legend.key = element_blank()) ;
                                                plot = plot + scale_fill_manual( name = legend.title,
                                                                                 limits = rev(NMVOC.levels),
                                                                                 labels = names,
                                                                                 values = colours) ;
                                                return(plot) } `,
    );

    ###tagged plot
    $R->run(q` tagged.data = data.frame(time)`);
    foreach my $ref (@$tagged_plot_data) {
        for my $key (keys %$ref) {
            #next if ($key eq "O3 + OH = HO2");
            my @value = @{ $ref->{$key} };
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@value]);
            $R->run(q` tagged.data[name] = rate`); 
        }
    }
    $R->run(q` tagged.colours = c( "Production Others" = "#696537", "NC7H16" = "#f9c600", "EBENZ" = "#76afca", "BENZENE" = "#dc3522", "OXYL" = "#8c6238", "PXYL" = "#9bb08f", "NC6H14" = "#8b1537", "IC4H10" = "#e7e85e", "C3H6" = "#0352cb", "C2H6" = "#86b650", "MXYL" = "#f3aa7f", "C5H8" = "#ee6738", "C2H4" = "#c65d6c", "NC5H12" = "#8ed6d5", "C3H8" = "#6c254f", "TOLUENE" = "#58691b", "NC4H10" = "#888a87", "IC5H12" = "#0e5c28", "CH4" = "#b569b3", "CO + OH = HO2" = "#2c9def", "O3 + OH = HO2" = "#603912" ) `,
            q` tagged.names = c( "NC7H16" = "Heptane", "EBENZ" = "Ethylbenzene", "BENZENE" = "Benzene", "OXYL" = "o-Xylene", "PXYL" = "p-Xylene", "NC6H14" = "Hexane", "IC4H10" = "2-Methylpropane", "C3H6" = "Propene", "C2H6" = "Ethane", "MXYL" = "m-Xylene", "C5H8" = "Isoprene", "C2H4" = "Ethene", "NC5H12" = "Pentane", "C3H8" = "Propane", "TOLUENE" = "Toluene", "NC4H10" = "Butane", "IC5H12" = "2-Methylbutane", "CH4" = "Methane" ) `,
    );
    $R->set('tagged.title', "(b) MCM v3.2 Ox Budget: Tagged Model Run");
    $R->set('tagged.legend', "MCM v3.2 Tagged");
    $R->set('tagged.y.max', 1.6e9);
    $R->set('tagged.y.breaks', 2e8);
    $R->run(q` tagged.plot = plot.function(tagged.data, tagged.colours, tagged.names, tagged.title, tagged.legend, tagged.y.max, tagged.y.breaks) `);
    #my $p = $R->run(q` print(tagged.data) `);
    #print "$p\n";

    ###non-tagged plot
    $R->run(q` non.tagged.data = data.frame(time)`);
    foreach my $ref (@$plot_data) {
        for my $key (keys %$ref) {
            #next if ($key eq "O3 + OH = HO2");
            my @value = @{ $ref->{$key} };
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@value]);
            $R->run(q` non.tagged.data[name] = rate`); 
        }
    }
    $R->run(q` non.tagged.colours = c( "Production Others" = "#696537", "C2H5O = CH3CHO + HO2" = "#f9c600", "C2H5O2 + NO = C2H5O + NO2" = "#76afca", "HCHO + hv = CO + HO2 + HO2" = "#dc3522", "CH3CO3 + NO = CH3O2 + NO2" = "#8c6238", "HCHO + OH = CO + HO2" = "#9bb08f", "CH3O = HCHO + HO2" = "#8b1537", "CH3O2 + NO = CH3O + NO2" = "#e7e85e", "CO + OH = HO2" = "#2c9def", "O3 + OH = HO2" = "#603912"  ) `,
            q` non.tagged.names = c( "Production Others" = "Production Others", "C2H5O = CH3CHO + HO2" = "C2H5O = CH3CHO + HO2", "C2H5O2 + NO = C2H5O + NO2" = "C2H5O2 + NO = C2H5O + NO2", "HCHO + hv = CO + HO2 + HO2" = "HCHO + hv = CO + HO2 + HO2", "CH3CO3 + NO = CH3O2 + NO2" = "CH3CO3 + NO = CH3O2 + NO2", "HCHO + OH = CO + HO2" = "HCHO + OH = CO + HO2", "CH3O = HCHO + HO2" = "CH3O = HCHO + HO2", "CH3O2 + NO = CH3O + NO2" = "CH3O2 + NO = CH3O + NO2", "CO + OH = HO2" = "CO + OH = HO2" ) `);
    $R->set('non.tagged.title', "(a) MCM v3.2 Ox Budget: Non-Tagged Model Run");
    $R->set('non.tagged.legend', "MCM v3.2 Non-Tagged");
    $R->set('non.tagged.y.max', 1.6e9);
    $R->set('non.tagged.y.breaks', 2e8);
    $R->run(q` non.tagged.plot = plot.function(non.tagged.data, non.tagged.colours, non.tagged.names, non.tagged.title, non.tagged.legend, non.tagged.y.max, non.tagged.y.breaks) `);
    #my $p1 = $R->run(q` print(non.tagged.data) `);
    #print "$p1\n";

     $R->run(q` CairoPDF(file = "Slice_MCM_Ox_tagged_non_tagged_production.pdf", width = 120, height = 60) `, 
             q` y.label = textGrob(expression(bold(paste("\nRate (molecules ", cm^-3, s^-1, ")"))), rot = 90, gp = gpar(fontsize = 110))`,
             q` main.plot = grid.arrange(y.label, 
                                         arrangeGrob(non.tagged.plot, 
                                                     tagged.plot + theme(axis.text.y = element_blank()) + scale_y_continuous(breaks = NULL), 
                                                     nrow = 1), 
                                        nrow = 1, ncol = 2,
                                        sub = textGrob("\n", gp = gpar(fontsize = 20)), 
                                        widths=unit.c(unit(10, "lines"), unit(1, "npc") - unit(10, "lines") )) `,
             q` print(main.plot) `,
             q` dev.off() `,
     );

    $R->stop(); 
}
