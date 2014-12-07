#! /usr/bin/env perl
# Compare First day and final cumulative TOPP for each VOC in each mechanism in the constant emissions and constant mixing ratio runs
# Version 0: Jane Coates 7/12/2014

use strict;
use diagnostics;
use Statistics::R;
use Cwd qw (cwd);

my $base_dir = cwd();
#get constant mixing ratios first day  TOPP
opendir DIR, $base_dir or die "Can't open $base_dir : $!";
my @daily_constant_MR = grep { $_ =~ /_TOPPs_MR/ } readdir DIR;
closedir DIR;

my %MR_TOPP_first_day;
foreach my $daily_MR (@daily_constant_MR) {
    my @lines = split /\n/, read_file($daily_MR);
    (my $mechanism = $daily_MR) =~ s/^(.*?)_TOPPs_MR\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $TOPPs =~ s/\[|\]//g;
        $TOPPs =~ s/^\s+|\s+$//g;
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $MR_TOPP_first_day{$mechanism}{$VOC} = $TOPPs[0];
    }
}

#get constant mixing ratios cumulative last day  TOPP
opendir DIR, $base_dir or die "Can't open $base_dir : $!";
my @cumulative_constant_MR = grep { $_ =~ /_cumulative/ } readdir DIR;
closedir DIR;

my %MR_TOPP_cumulative;
foreach my $cumulative_MR (@cumulative_constant_MR) {
    my @lines = split /\n/, read_file($cumulative_MR);
    (my $mechanism = $cumulative_MR) =~ s/^(.*?)_cumulative\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $TOPPs =~ s/\[|\]//g;
        $TOPPs =~ s/^\s+|\s+$//g;
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $MR_TOPP_cumulative{$mechanism}{$VOC} = $TOPPs[6];
    }
}

#get constant emissions first day and calculate cumulative sum TOPP
(my $dir = $base_dir) =~ s/(.*?)Emissions_vs_Mixing_ratio/$1/;
opendir DIR, $dir or die "Can't open $dir : $!";
my @daily_constant_emissions = grep { $_ =~ /_TOPP_values/ } readdir DIR;
closedir DIR;

my (%emissions_TOPP_first_day, %emissions_TOPP_cumulative);
foreach my $daily_emissions (@daily_constant_emissions) {
    my @lines = split /\n/, read_file("$dir/$daily_emissions");
    (my $mechanism = $daily_emissions) =~ s/^(.*?)_TOPP_values\.txt/$1/;
    foreach my $line (@lines) {
        next if ($line =~ /^Working|^CH4/);
        my ($VOC, $TOPPs) = split / => /, $line;
        $TOPPs =~ s/\[|\]//g;
        $TOPPs =~ s/^\s+|\s+$//g;
        $TOPPs =~ s/\s+/:/g;
        my @TOPPs = split /:/, $TOPPs;
        $emissions_TOPP_first_day{$mechanism}{$VOC} = $TOPPs[0];
        my $sum = 0;
        $sum += $_ foreach (@TOPPs);
        $emissions_TOPP_cumulative{$mechanism}{$VOC} = $sum;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` plotting = function (data, filename, title) {plot = ggplot(data, aes(x = Emissions, y = Mixing.Ratio, colour = VOC, group = VOC));
                                                        plot = plot + geom_point();
                                                        plot = plot + xlab("Constant Emissions TOPP");
                                                        plot = plot + ylab("Constant Mixing Ratio TOPP");
                                                        plot = plot + ggtitle(title);
                                                        plot = plot + geom_abline(intercept = 0, slope = 1);
                                                        plot = plot + theme_bw();
                                                        plot = plot + theme(legend.key = element_blank());
                                                        plot = plot + theme(axis.title = element_text(face = "bold"));
                                                        plot = plot + theme(plot.title = element_text(face = "bold"));
                                                        plot = plot + theme(panel.border = element_rect(colour = "black"));
                                                        CairoPDF(file = filename);
                                                        print(plot);
                                                        dev.off() } `);

foreach my $mechanism (sort keys %emissions_TOPP_first_day) {
    $R->set('daily.filename', "${mechanism}_first_day.pdf");
    $R->set('daily.title', "$mechanism : First Day TOPP Values");
    $R->set('cumulative.filename', "${mechanism}_cumulative.pdf");
    $R->set('cumulative.title', "$mechanism : Cumulative TOPP Values");
    my @VOCs = sort keys %{$emissions_TOPP_first_day{$mechanism}};
    $R->set('voc', [@VOCs]);
    $R->run(q` daily.data = data.frame(Mixing.Ratio = as.numeric(0), Emissions = as.numeric(0)) `,
            q` cumulative.data = data.frame(Mixing.Ratio = as.numeric(0), Emissions = as.numeric(0)) `,
    );
    foreach my $VOC (sort keys %{$emissions_TOPP_first_day{$mechanism}}) {
        $R->set('daily.MR', $MR_TOPP_first_day{$mechanism}{$VOC});
        $R->set('daily.emissions', $emissions_TOPP_first_day{$mechanism}{$VOC});
        $R->set('cumulative.MR', $MR_TOPP_cumulative{$mechanism}{$VOC});
        $R->set('cumulative.emissions', $emissions_TOPP_cumulative{$mechanism}{$VOC}); 
        $R->run(q` daily.data = rbind(daily.data, c(daily.MR, daily.emissions)) `);
        $R->run(q` cumulative.data = rbind(cumulative.data, c(cumulative.MR, cumulative.emissions)) `);
    }
    $R->run(q` daily.data = daily.data[-1,] `,
            q` daily.data$VOC = voc `,
            q` plotting(daily.data, daily.filename, daily.title) `,
            q` cumulative.data = cumulative.data[-1,] `,
            q` cumulative.data$VOC = voc `,
            q` plotting(cumulative.data, cumulative.filename, cumulative.title) `,
    );
    #my $p = $R->run(q` print(cumulative.data) `);
    #print $p, "\n";
}

$R->stop();

sub read_file {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    close $in;
    return $data;
}
