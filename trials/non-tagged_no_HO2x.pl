#!/usr/bin/perl -w
# calculate Ox budget from non-tagged model run and plot as stacked bar plot - no HO2x
# Version 0: Jane Coates 05/05/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA_TIM;
use KPP;
use Statistics::R;

my $run = "/work/users/jco/MECCA/MCM_3.2_7_days_no_tagging/boxmodel";
my $mecca = MECCA_TIM->new($run); 
my $eqnfile = "/work/users/jco/MECCA/MCM_3.2_7_days_no_tagging/gas.eqn";
my $kpp = KPP->new($eqnfile); 
my $ro2file = "/work/users/jco/MECCA/MCM_3.2_7_days_no_tagging/RO2_species.txt";
my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2file);

my %families = (
    'Ox' => [ qw(O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ],
);

my %weights = (
    'Ox' => { NO3 => 2, N2O5 => 3},
);
my %production_reaction_rates;
my $ntime = $mecca->time->nelem;

foreach my $species (qw( Ox )) { #get all production and consumption rates
    my ($producers, $producer_yields);
    if (exists $families{$species}) { #get family reaction numbers and yields
        $kpp->family({ 
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers);  
    } else { #get reaction numbers and yields
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers);  
    }

    #check that species reactions are found
    die "No producers found for $species\n" if (@$producers == 0);
    
    for (0..$#$producers) { #get rates for each producing reactions
        my $reaction = $producers->[$_];
		my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur
        my $string = $kpp->reaction_string($reaction);
        $production_reaction_rates{$species}{$string} += $rate(1:$ntime-2);
    }
} 

my (@production_others, %production_plot_rates);
my $prod_others_max = 1e8;

#create hash with production of the reactions
my (%prod_hash, @sorted_plot_data);
$prod_hash{$_} += $production_reaction_rates{'Ox'}{$_} for (keys %{ $production_reaction_rates{'Ox'} });

#sort production
my $sort_function = sub { $_[0]->sum };
foreach my $item (keys %prod_hash) { 
    if ($prod_hash{$item}->sum < $prod_others_max) { #get production others
        $prod_hash{'Production Others'} += $prod_hash{$item};
        delete $prod_hash{$item};
    }
}

$prod_hash{$_}->where($prod_hash{$_} < 0) .= 0 foreach (keys %prod_hash); 
my @sorted_prod = sort { &$sort_function($prod_hash{$b}) <=> &$sort_function($prod_hash{$a}) } keys %prod_hash;


foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
    next if ($_ eq 'Production Others');
    push @sorted_plot_data, { $_ => $prod_hash{$_} };
}

push @sorted_plot_data, { 'Production Others' => $prod_hash{'Production Others'} } if (defined $prod_hash{'Production Others'}); #add Production Others to the beginning 

my @plot_data;
foreach my $ref (@sorted_plot_data) {#extract reaction and rates for each plot
    foreach my $item (keys %{$ref}) {
        my $rate_list = join ":", $ref->{$item}->dog;
        my @rate_array = split /:/, $rate_list;
        push @plot_data, { $item => \@rate_array };
    }
} 

#Create x-axis for plot in hours
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 3600;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list;
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

my ($reaction_plot_data) = plot(\@time_blocks, \@plot_data); 

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

sub plot { #create dataframe and then create plot
    my ($time, $plot_data) = @_;
    
    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(plyr) `);
    $R->run(q` library(reshape2) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(RColorBrewer) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(Cairo) `);
    
    $R->set('time', [@$time]);
    my $data_frame = $R->run(q` data = data.frame(time)`);
    foreach my $ref (@$plot_data) {
        for my $key (keys %$ref) {
            next if ($key eq "O3 + OH = HO2");
            my @value = @{ $ref->{$key} };
            my $R_name = $R->set('name', $key);
            my $R_data = $R->set('rate', [@value]);
            $R->run(q` data[name] = rate`); 
        }
    } 

    $R->run(q` data = ddply(data, .(time), colwise(sum))`); 
    $R->run(#reshape data frames
            q` plot.data = melt( data = data, 
                                 id = names(data)[1],
                                 measured = names(data)[-1] )`,
            q` colnames(plot.data) = c("time", "reaction", "rate") `,

            #order factor levels for plot
            q` reaction.levels = (levels(factor(plot.data$reaction))) `,
            q` plot.data$reaction = ordered(plot.data$reaction, levels = reaction.levels) `,
            q` plot.data = ddply( plot.data, .(reaction)) `,

            #set colours
            q` my.colours = c( "Production Others" = "#696537", "C2H5O = CH3CHO + HO2" = "#f9c600", "C2H5O2 + NO = C2H5O + NO2" = "#76afca", "HCHO + hv = CO + HO2 + HO2" = "#dc3522", "CH3CO3 + NO = CH3O2 + NO2" = "#8c6238", "HCHO + OH = CO + HO2" = "#9bb08f", "CH3O = HCHO + HO2" = "#8b1537", "CH3O2 + NO = CH3O + NO2" = "#e7e85e", "CO + OH = HO2" = "#2c9def" ) `,

            q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis

            #plot 
            q` CairoPDF(file = "non-tagged_no_HO2x.pdf", width = 70, height = 40) `,
            q` plot = ggplot(plot.data, aes(x = time, y = rate, fill = reaction)) `,
            q` plot = plot + geom_bar(stat = "identity", width = 0.5) `,
            q` plot = plot + theme(legend.text = theme_text(size = 60)) `,
            q` plot = plot + theme_bw() `,
            q` plot = plot + theme(axis.text.x = element_text(size = 70)) `,
            q` plot = plot + theme(axis.text.y = element_text(size = 70)) `,
            q` plot = plot + theme(axis.title.x = element_blank()) `,
            q` plot = plot + theme(legend.key.size = unit(3, "cm")) `,
            q` plot = plot + theme(axis.title.y = element_text(size = 80)) `,
            q` plot = plot + ylab(expression(bold(paste("Rate (molecules ", cm^-3, s^-1, " x ", 10^6, ")")))) `,
            q` plot = plot + theme(legend.text = element_text(size = 60)) `,
            q` plot = plot + theme(legend.title = element_text(size = 80, face = "bold")) `,
            q` plot = plot + scale_x_discrete(limits = c("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")) `,
            q` plot = plot + scale_y_continuous(limits=c(0, 1.4e9), breaks=seq(0, 1.4e9, 2e8), label = scientific_10)`,
            q` plot = plot + scale_fill_manual( name = "MCM v3.2 non-Tagged", 
                                                limits = rev(reaction.levels),
                                                values = my.colours) `,
            q` print(plot) `,
            q` dev.off() `,
    );

    $R->stop(); 
}

