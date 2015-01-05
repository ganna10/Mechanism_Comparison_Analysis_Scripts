#R script to plot first day TOPP vs carbon number of VOC
# Version 0: Jane Coates 5/1/2015

#first day TOPPs for all alkanes
cb05.data = data.frame(TOPP = c(0.0845172479748726, 0.3641597032547, 0.952807605266571, 0.95486788938202, 1.16583390916096, 1.16891467571259, 1.37504136562347, 1.57824212238862, 1.7769056558609), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.data = data.frame(TOPP = c(0.0869439467787743, 0.321530640125275, 0.840508162975311, 0.842413246631622, 1.0273551940918, 1.03020417690277, 1.21066904067993, 1.38842662532035, 1.56198441982269), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cri.data = data.frame(TOPP = c(0.114152513444424, 0.354555159807205, 0.972761750221252, 0.967499315738678, 1.92756509780884, 1.57766127586365, 2.27819347381592, 2.37771141645028, 2.62646389007568), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.data = data.frame(TOPP = c(0.12388314306736, 0.381832748651505, 1.03673362731934, 1.00873160362244, 1.7392680644989, 1.60785245895386, 2.22197437286377, 2.43809782036837, 2.71707463264465), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.data = data.frame(TOPP = c(0.114175207912922, 0.346426695585251, 0.969126343727112, 0.97936624288559, 1.65823602676392, 1.56109631061554, 2.15132093429565, 2.35366020437529, 2.59780859947205), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mozart.data = data.frame(TOPP = c(0.13925464451313, 0.426944702863693, 1.71256327629089, 1.71256327629089, 1.37005066871643, 1.37005066871643, 1.1417088508606, 0.978607594966888, 0.856281638145447), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
radm2.data = data.frame(TOPP = c(0.128240107806014, 0.877153517802556, 0.657865138351917, 0.657865138351917, 1.42851711273193, 1.42851711273193, 1.19043092727661, 1.0203693662371, 2.58092876646804), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm.data = data.frame(TOPP = c(0.131782606244087, 0.7080078125, 0.531005859375, 0.531005859375, 1.34979560852051, 1.34979560852051, 1.12482967376709, 0.964139720371791, 2.11075333331641), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm2.data = data.frame(TOPP = c(0.11060257256031, 0.996819305419922, 0.747614479064941, 0.747614479064941, 1.5016316986084, 1.5016316986084, 1.25135974884033, 1.07259407043457, 2.43148993558772), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))

full = rbind(cb05.data, cbm4.data, cri.data, mcm1.data, mcm2.data, mozart.data, radm2.data, racm.data, racm2.data)
full$Mechanism = factor(full$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05"))

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

plot = ggplot(full, aes(x = C, y = TOPP))
plot = plot + geom_point()
plot = plot + geom_smooth(method = "lm", se = FALSE, colour = "black")
plot = plot + geom_text(data = eqn, aes(x = 4.7, y = 3, label = V1), size = 3, parse = TRUE)
plot = plot + facet_wrap( ~ Mechanism)
plot = plot + ylab("TOPP (molecules(Ox) / molecules(VOC))")
plot = plot + xlab("VOC Carbon Number")
plot = plot + scale_x_continuous(limits = c(2, 8), breaks = seq(2, 8, 1))
plot = plot + theme_bw()
plot = plot + theme(axis.title = element_text(face = "bold"))
plot = plot + theme(strip.background = element_blank())
plot = plot + theme(strip.text = element_text(face = "bold"))
plot = plot + theme(panel.border = element_rect(colour = "black"))
plot = plot + theme(panel.grid = element_blank())

CairoPDF(file = "Alkanes_vs_C.pdf", width = 9, height = 12.7)
print(plot)
dev.off()
