#! /usr/bin/env perl
# compare radical production as prescribed by NO source calculation
# Version 0: Jane Coates 3/10/2014

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
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 3600;
my @time_axis = map { $_ } $times->dog; 
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

my @runs = qw( MCM_3.2_tagged MCM_3.1_tagged_3.2rates CRI_tagging MOZART_tagging RADM2_tagged RACM_tagging RACM2_tagged CBM4_tagging CB05_tagging );
my @mechanisms = ( "(a) MCMv3.2", "(b) MCMv3.1", "(c) CRIv2", "(g) MOZART-4", "(d) RADM2", "(e) RACM", "(f) RACM2", "(h) CBM-IV", "(i) CB05" );
#my @runs = qw( CBM4_tagging CB05_tagging );
#my @mechanisms = qw( CBM-IV CB05 );
my $index = 0;

my (%families, %weights, %plot_data, %legend);
foreach my $run (@runs) {
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $radical_file = "$base/$run/radicals.txt";
    my @radicals = get_species($radical_file);
    $families{$mechanisms[$index]} = [ @radicals ];
    ($plot_data{$mechanisms[$index]}, $legend{$mechanisms[$index]}) = get_data($mecca, $kpp, $mechanisms[$index]);
    $index++;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(plyr) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(gridExtra) `,
        q` library(grid) `,
);

$R->run(q` my.colours = c(  "Production Others" = "#696537",
                            "MGLYOX + hv" = "#f9c500", "CH3COCHO + hv" = "#f9c500", "CARB6 + hv" = "#f9c500", "MGLY + hv" = "#f9c500", 
                            "NO + TCO3" = "#1c3e3d",
                            "HC5 + OH" = "#c9a415",
                            "EPX + O3" = "#6db875",
                            "DCB1 + O3" = "#ae4901",
                            "BIGALD + hv" = "#8fd5d3",
                            "OH + PAR" = "#bb8a01", 
                            "C2O3 + NO" = "#e7e85e",
                            "ROR" = "#0352cb", 
                            "C2H4 + OH" = "#9bb08f", 
                            "CXO3 + NO" = "#ba8b01",
                            "C2H6 + OH" = "#8ed6d2",
                            "O1D" = "#6c254f",
                            "HCHO + hv" = "#0e5c28", "CH2O + hv" = "#0e5c28", "FORM + hv" = "#0e5c28",
                            "DCB + hv" = "#76afca",
                            "KET + hv" = "#a67c52",
                            "CH4 + OH" = "#0d3e76",
                            "MEK + hv" = "#86b650" ) `,
);

$R->run(q` plotting = function (data, legend, mechanism) {  plot = ggplot(data, aes(x = time, y = Rate, fill = Reaction)) ;
                                                            plot = plot + geom_bar(stat = "identity") ;
                                                            plot = plot + scale_y_continuous(limits = c(0, 6.5e8), breaks = seq(0, 6.5e8, 1e8)) ;
                                                            plot = plot + theme_bw() ;
                                                            plot = plot + ggtitle(mechanism) ;
                                                            plot = plot + xlab("\n");
                                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                                            plot = plot + theme(axis.title.x = element_text(size = 100)) ;
                                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                                            plot = plot + theme(legend.title = element_blank()) ;
                                                            plot = plot + theme(legend.key = element_blank()) ;
                                                            plot = plot + theme(plot.title = element_text(size = 200, face = "bold")) ;
                                                            plot = plot + theme(axis.text.x = element_text(size = 140)) ;
                                                            plot = plot + theme(axis.text.y = element_text(size = 140)) ;
                                                            plot = plot + theme(legend.text = element_text(size = 120)) ;
                                                            plot = plot + theme(legend.key.size = unit(7, "cm")) ;
                                                            plot = plot + theme(legend.justification = c(0.99, 0.99)) ;
                                                            plot = plot + theme(legend.position = c(0.99, 0.99)) ;
                                                            plot = plot + scale_fill_manual(values = my.colours, limits = rev(legend)) ;
                                                            return(plot) } `,
);

$R->set('time', [@time_blocks]);
$R->run(q` plots = c() `);
foreach my $run (sort keys %plot_data) {
    $R->run(q` data = data.frame(time) `);
    $R->set('legend', [@{$legend{$run}}] );
    $R->set('mechanism', $run );
    foreach my $ref (@{$plot_data{$run}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [map { $_ } $ref->{$reaction}->dog]);
            $R->run(q` data[reaction] = rate `);
        }
    }

    $R->run(q` data = ddply(data, .(time), colwise(sum)) `,
            q` data = data[1:7,] `,
            q` data = melt(data, id.vars = c("time"), variable.name = "Reaction", value.name = "Rate") `,
            q` plot = plotting(data, legend, mechanism) `,
            q` plots = c(plots, list(plot)) `,
    );
}

$R->run(q` CairoPDF(file = "radical_production_analysis.pdf", width = 141, height = 200) `,
        q` multiplot = grid.arrange(arrangeGrob(plots[[1]],
                                                plots[[2]],
                                                plots[[3]],
                                                plots[[4]],
                                                plots[[5]],
                                                plots[[6]],
                                                plots[[7]],
                                                plots[[8]],
                                                plots[[9]],
                                                nrow = 3),
                                    nrow = 1, ncol = 1,
                                    left = textGrob(expression(bold(paste("Reaction Rate (molecules ", cm^-3, s^-1, ")"))), rot = 90, gp = gpar(fontsize = 140), vjust = 0.5)) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $family) = @_;
    
    $families{'NOx'} = [ qw( NO NO2 NO3 N2O5 ) ];
    $weights{'NOx'} = { N2O5 => 2 };
    $kpp->family({ #NOx family
                    name    => 'NOx',
                    members =>$families{'NOx'},
                    weights => $weights{'NOx'},
    });

    my ($radical_producers, $radical_producer_yields, $NOx_yields, %production_rates);
    if (exists $families{$family}) { #radicals family
        $kpp->family({                                                                                                                           
                        name    => $family,
                        members => $families{$family},
                        weights => $weights{$family},
        }); 
        $radical_producers = $kpp->producing($family);
        $radical_producer_yields = $kpp->effect_on($family, $radical_producers);
        $NOx_yields = $kpp->effect_on('NOx', $radical_producers);
    } else {
        print "No radical family found for $family\n";
    }
    
    die "No producers found for $family\n" if (@$radical_producers == 0);

    my $string_max = 30;
    for (0..$#$radical_producers) { #get rates for all radical producing reactions
        next if ($NOx_yields->[$_] < 0);
        my $reaction = $radical_producers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $net_radical_yield = $radical_producer_yields->[$_] - $NOx_yields->[$_];
        next if ($net_radical_yield == 0);
        my $rate = $net_radical_yield * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $reaction_string =~ s/_(.*?)\b//g;
        $reaction_string =~ s/(.{1,$string_max})/$1\n/gs;
        $reaction_string =~ s/\n$//;
        $reaction_string = string_mapping($reaction_string);
        my ($reactants, $products) = split / = /, $reaction_string;
        $production_rates{$reactants} += $rate(1:$NTIME-2);
    }

    sub string_mapping {
        my ($string) = @_;
        if ($string =~ /MGLYOX \+ hv/) {
            $string = "MGLYOX + hv = CH3CO3 + CO\n+ HO2";
        } elsif ($string =~ /NO \+ TCO3/) {
            $string = "NO + TCO3 = 0.05 ACO3 + 0,95 CO\n+ 0.89 GLY + 0.92 HO2 + 0.11 MGLY\n+ NO2 + 2 XO2";
        } elsif ($string =~ /HC5 \+ OH/) {
            $string = "HC5 + OH = H2O + HC5P\n+ 0.25 XO2";
        } elsif ($string =~ /EPX \+ O3/) {
            $string = "EPX + O3 = 0.85 BALD + 1.5 CO\n+ 0.5 CO2 + GLY + 1.5 HO2 + 0.05 OH";
        } elsif ($string =~ /DCB1 \+ O3/) {
            $string = "DCB1 + O3 = 1.5 CO + 0.5 CO2 + 0.05 GLY\n+ 0.05 HCHO + HO2 + 0.08 MGLY + 0.05 OH\n+ 0.65 OP2 + 0.6 RCO3 + 0.6 XO2";
        } elsif ($string =~ /CH3COCHO \+ hv/) {
            $string = "CH3COCHO + hv = CH3CO3 + CO + HO2";
        } elsif ($string =~ /BIGALD \+ hv/) {
            $string = "BIGALD + hv = .13 CH3CO3 + .18 CH3COCHO\n+ .45 CO + .13 GLYOXAL + .56 HO2";
        } elsif ($string =~ /OH \+ PAR = 0\.11 ALD2/) {
            $string = "OH + PAR = 0.11 ALD2 + 0.11 HO2\n+ 0.11 PAROP + 0.76 ROR + 0.87 XO2 + 0.13 XO2N";
        } elsif ($string =~ /C2O3 \+ NO/) {
            $string = "C2O3 + NO = HCHO + HO2\n+ NO2 + XO2";
        } elsif ($string =~ /ROR = 1\.1 ALD2/) {
            $string = "ROR = 1.1 ALD2 + 0.94 HO2 + 2.1 PAROP\n+ 0.02 ROR + 0.96 XO2 + 0.04 XO2N";
        } elsif ($string =~ /ETH \+ OH = 0\.22 ALD2/) {
            $string = "C2H4 + OH = 0.22 ALD2 + 1.56 HCHO\n+ HO2 + XO2";
        } elsif ($string =~ /OH \+ PAR = 0\.060 ALD2/) {
            $string = "OH + PAR = 0.06 ALD2 + 0.05 ALDX + 0.11 HO2\n+ 0.11 PAROP + 0.76 ROR + 0.87 XO2 + 0.13 XO2N";
        } elsif ($string =~ /ROR = 0\.600 ALD2/) {
            $string = "ROR = 0.6 ALD2 + 0.5 ALDX + 0.94 HO2 + 2.1 PAROP\n+ 0.02 ROR + 0.96 XO2 + 0.04 XO2N";
        } elsif ($string =~ /CXO3 \+ NO/) {
            $string = "CXO3 + NO = ALD2 + HO2\n+ NO2 + XO2";
        } elsif ($string =~ /ETH \+ OH = 0\.220 ALDX/) {
            $string = "C2H4 + OH = 0.22 ALDX + 1.56 FORM\n+ HO2 + XO2";
        } elsif ($string =~ /ETHA \+ OH/) {
            $string = "C2H6 + OH = 0.991 ALD2 + HO2\n+ 0.991 XO2 + 0.009 XO2N";
        } elsif ($string =~ /HCHO \+ hv = CO \+ 2\.0 HO2/) {
            $string = "HCHO + hv = CO + 2 HO2";
        }
        return $string;
    }

    my $others = 3e7;
    foreach my $reaction (keys %production_rates) {
        if ($production_rates{$reaction}->sum < $others) {
            $production_rates{'Production Others'} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
    }

    my $sort_function = sub { $_[0]->sum };
    my (@sorted_plot_data, @legend);
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @legend, $_;
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); #add Production Others to the beginning 
    push @legend, 'Production Others';
    return (\@sorted_plot_data, \@legend);
}

sub get_species {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die $!;
    local $/ = undef;
    my $lines = <$in>;
    close $in;
    my @species = split /\s+/, $lines;
}
