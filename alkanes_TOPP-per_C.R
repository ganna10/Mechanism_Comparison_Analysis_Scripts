#R script to plot first day TOPP per carbon number vs alkane C number
# Version 0: Jane Coates 5/1/2015

#first day  TOPP for all alkanes
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
library(dplyr)
library(Cairo)

my.colours = c( "Ethane" = "#696537", "Propane" = "#f9c600", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb" )

plot = ggplot(full, aes(x = Mechanism, y = TOPP, colour = VOC))
plot = plot + geom_point()
plot = plot + coord_flip()
plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC)) per Carbon Number")
plot = plot + scale_colour_manual(values = my.colours)
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

#cumulative  TOPP for all alkanes
cb05.sums = data.frame(TOPP = c(0.60335508361459, 0.923657203714053, 1.84673887491226, 1.84680387749839, 1.84610499031539, 1.84618226289749, 1.8456027507782, 1.84513571958547, 1.84471476078034), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.sums = data.frame(TOPP = c(0.257339833304286, 0.642852172255517, 1.28430606424809, 1.28446382284164, 1.28276613950729, 1.28295435905457, 1.28154492378235, 1.28041032181179, 1.27938625961543), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
cri.sums = data.frame(TOPP = c(0.90744724497199, 1.19261922438939, 2.02547438442707, 1.30156747996807, 2.15941095352172, 1.48752684593201, 2.19520086050033, 2.16198409962616, 2.15395130589604), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.sums = data.frame(TOPP = c(0.95420465618372, 1.24286951621373, 2.04508712887764, 1.33807555586099, 2.09400859475136, 1.63178697228432, 2.10139422118663, 2.09759336875754, 2.03783712722361), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.sums = data.frame(TOPP = c(0.865593794733285, 1.11035437385241, 1.9816093146801, 1.25637395679951, 2.0586550951004, 1.56194756031036, 2.07305474579333, 2.04698729770051, 1.97013106569648), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
mozart.sums = data.frame(TOPP = c(0.88885911554098, 1.05732961495717, 1.66577024757862, 1.66577024757862, 1.66577024757862, 1.66577024757862, 1.66577024757862, 1.66577024757862, 1.66577024757862), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
radm2.sums = data.frame(TOPP = c(1.01430642907959, 1.7917226785901, 1.791722640909, 1.791722640909, 1.47968599696954, 1.47968599696954, 1.47968599696954, 1.47968599696954, 1.24018539668885), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
racm.sums = data.frame(TOPP = c(1.00219408422708, 1.82024583049204, 1.82024583049204, 1.82024583049204, 1.62936080495516, 1.62936080495516, 1.62936080495516, 1.62936080495516, 1.00811153726678), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))
racm2.sums = data.frame(TOPP = c(0.909124590456485, 1.42925610696828, 1.42925610696828, 1.42925610696828, 1.1370189424072, 1.1370189424072, 1.1370189424072, 1.1370189424072, 1.05259203890672), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane"))

sum = rbind(cb05.sums, cbm4.sums, cri.sums, mcm1.sums, mcm2.sums, mozart.sums, radm2.sums, racm.sums, racm2.sums)
sum = filter(sum, VOC == c("Hexane", "Heptane", "Octane"))
sum$Mechanism = factor(sum$Mechanism, levels = rev(c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05")))
sum$VOC = factor(sum$VOC, levels = c("Hexane", "Heptane", "Octane"))

plot = ggplot(sum, aes(x = Mechanism, y = TOPP, colour = VOC))
plot = plot + geom_point()
plot = plot + coord_flip()
plot = plot + ylab("TOPP (molecules(Ox)/molecules(VOC)) per Carbon Number")
plot = plot + scale_colour_manual(values = my.colours)
plot = plot + theme_bw()
plot = plot + theme(legend.direction = "horizontal")
plot = plot + theme(legend.position = c(0.27, 0.93))
plot = plot + theme(axis.title = element_text(face = "bold"))
plot = plot + theme(panel.border = element_rect(colour = "black"))
plot = plot + theme(panel.grid = element_blank())
plot = plot + theme(legend.title = element_blank())
plot = plot + theme(legend.key = element_blank())
plot = plot + theme(axis.title.y = element_blank())

CairoPDF(file = "Alkanes_vs_C_sums.pdf", width = 6, height = 4.2)
print(plot)
dev.off()
