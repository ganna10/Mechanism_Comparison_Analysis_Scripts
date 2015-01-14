#R script to plot first day TOPP per carbon number vs alkane C number
# Version 0: Jane Coates 5/1/2015

#maximum TOPP for all alkanes
cb05.data = data.frame(TOPP = c(0.0422586239874363, 0.121386567751567, 0.238201901316643, 0.238716972345505, 0.233166781832192, 0.233782935142518, 0.229173560937245, 0.225463160341231, 0.222113206982612), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.data = data.frame(TOPP = c(0.0434719733893872, 0.107176880041758, 0.210127040743828, 0.210603311657905, 0.20547103881836, 0.206040835380554, 0.201778173446655, 0.19834666076005, 0.195248052477836), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
cri.data = data.frame(TOPP = c(0.057076256722212, 0.118185053269068, 0.243190437555313, 0.241874828934669, 0.385513019561768, 0.31553225517273, 0.379698912302653, 0.339673059492897, 0.32830798625946), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.data = data.frame(TOPP = c(0.06194157153368, 0.127277582883835, 0.259183406829835, 0.25218290090561, 0.34785361289978, 0.321570491790772, 0.370329062143962, 0.348299688624053, 0.339634329080581), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.data = data.frame(TOPP = c(0.057087603956461, 0.115475565195084, 0.242281585931778, 0.244841560721398, 0.331647205352784, 0.31221923828125, 0.358553489049275, 0.336237172053613, 0.324726074934006), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mozart.data = data.frame(TOPP = c(0.069627322256565, 0.142314900954564, 0.274010121822358, 0.274010121822358, 0.274010133743286, 0.274010133743286, 0.274010133743286, 0.274010133743286, 0.274010133743286), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
radm2.data = data.frame(TOPP = c(0.064120053903007, 0.312896624378774, 0.312896617527665, 0.312896617527665, 0.31000804901123, 0.31000804901123, 0.31000804901123, 0.31000804901123, 0.330835284918191), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
racm.data = data.frame(TOPP = c(0.0658913031220435, 0.252559267241379, 0.25255926724138, 0.25255926724138, 0.292924404144288, 0.292924404144288, 0.292924404144288, 0.292924404144288, 0.270566041764641), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
racm2.data = data.frame(TOPP = c(0.055301286280155, 0.230745209587945, 0.230745209587945, 0.230745209587945, 0.23941831929343, 0.23941831929343, 0.23941831929343, 0.23941831929343, 0.31167953027883), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))

full = rbind(cb05.data, cbm4.data, cri.data, mcm1.data, mcm2.data, mozart.data, radm2.data, racm.data, racm2.data)
full$Mechanism = factor(full$Mechanism, levels = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")))
full$VOC = factor(full$VOC, levels = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))

#plotting 
library(ggplot2)
library(plyr)
library(Cairo)

lm_eqn = function(df){
    m = lm(TOPP ~ C, df);
    eq <- substitute(italic(TOPP) == a + b %.% italic(C)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 3), 
              b = format(coef(m)[2], digits = 3), 
             r2 = format(summary(m)$r.squared, digits = 4)))
    as.character(as.expression(eq));
}
eqn = ddply(full, .(Mechanism), lm_eqn)

my.colours = c( "Ethane" = "#696537", "Propane" = "#f9c600", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb" )

plot = ggplot(full, aes(x = Mechanism, y = TOPP, colour = VOC))
plot = plot + geom_point()
plot = plot + coord_flip()
plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC)) per Carbon Number")
plot = plot + scale_colour_manual(values = my.colours)
#plot = plot + scale_x_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1))
plot = plot + theme_bw()
plot = plot + theme(axis.title = element_text(face = "bold"))
plot = plot + theme(panel.border = element_rect(colour = "black"))
plot = plot + theme(panel.grid = element_blank())
plot = plot + theme(legend.title = element_blank())
plot = plot + theme(legend.key = element_blank())
plot = plot + theme(axis.title.y = element_blank())

CairoPDF(file = "Alkanes_vs_C.pdf", width = 6, height = 4.2)
print(plot)
dev.off()
