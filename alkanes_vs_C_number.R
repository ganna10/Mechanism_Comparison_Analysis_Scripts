#R script to plot first day TOPP vs carbon number of VOC
# Version 0: Jane Coates 5/1/2015

#maximum TOPP for all alkanes
cb05.data = data.frame(TOPP = c(0.0845172479748726, 0.3641597032547, 0.952807605266571, 0.952807605266571, 1.16583390916096, 1.16583390916096, 1.37504136562347, 1.57824206352234, 1.7769056558609), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.data = data.frame(TOPP = c(0.0869439467787743, 0.321530640125275, 0.840508162975311, 0.842413246631622, 1.0273551940918, 1.03020417690277, 1.21066904067993, 1.38842666149139, 1.56198441982269), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cri.data = data.frame(TOPP = c(0.114152513444424, 0.354555159807205, 0.972761750221252, 0.967499315738678, 1.92756509780884, 1.57766127586365, 2.27819347381592, 2.37771153450012, 2.62646389007568), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.data = data.frame(TOPP = c(0.12388314306736, 0.381832748651505, 1.03673362731934, 1.00873160362244, 1.7392680644989, 1.60785245895386, 2.22197437286377, 2.43809771537781, 2.71707463264465), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.data = data.frame(TOPP = c(0.114175207912922, 0.346426695585251, 0.969126343727112, 0.97936624288559, 1.65823602676392, 1.56109619140625, 2.15132093429565, 2.35366010665894, 2.59780859947205), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mozart.data = data.frame(TOPP = c(0.13925464451313, 0.426944702863693, 1.09604048728943, 1.09604048728943, 1.37005066871643, 1.37005066871643, 1.64406085014343, 1.91807103157043, 2.19208097457886), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
radm2.data = data.frame(TOPP = c(0.128240107806014, 0.938689873136323, 1.25158647011066, 1.25158647011066, 1.55004024505615, 1.55004024505615, 1.86004837354024, 2.17005630334218, 2.64668227934553), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm.data = data.frame(TOPP = c(0.131782606244087, 0.757677801724138, 1.01023706896552, 1.01023706896552, 1.46462202072144, 1.46462202072144, 1.75754626592, 2.05047070980072, 2.16452833411713), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm2.data = data.frame(TOPP = c(0.11060257256031, 0.692235628763835, 0.922980838351779, 0.922980838351779, 1.19709159646715, 1.19709159646715, 1.43650991576059, 1.67592832020351, 2.49343624223064), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))

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
plot = plot + scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1))
plot = plot + theme_bw()
plot = plot + theme(axis.title = element_text(face = "bold"))
plot = plot + theme(strip.background = element_blank())
plot = plot + theme(strip.text = element_text(face = "bold"))
plot = plot + theme(panel.border = element_rect(colour = "black"))
plot = plot + theme(panel.grid = element_blank())

CairoPDF(file = "Alkanes_vs_C.pdf", width = 9, height = 12.7)
print(plot)
dev.off()

#cumulative data
cb05.sum = data.frame(TOPP = c(1.20671016722918, 2.77097161114216, 7.38695549964905, 7.38695549964905, 9.23091131448746, 9.23091131448746, 11.0736165046692, 12.9159500370983, 14.7577180862427), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9))
cbm4.sum = data.frame(TOPP = c(0.514679666608572, 1.92855651676655, 5.13722425699234, 5.13722425699234, 6.41477179527283, 6.41477179527283, 7.68926954269408, 8.96287225268254, 10.2350900769234), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9))
cri.sum = data.frame(TOPP = c(1.81489448994398, 3.57785767316818, 8.10189753770828, 8.10189753770828, 7.43763422966003, 7.43763422966003, 13.171205163002, 15.1338886973831, 17.2316104471683), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9))
mcm1.sum = data.frame(TOPP = c(1.90840931236744, 3.7286085486412, 8.18034851551057, 8.18034851551057, 10.4700429737568, 10.4700429737568, 12.6083653271198, 14.6831535813028, 16.3026970177889), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9))
mcm2.sum = data.frame(TOPP = c(1.73118758946657, 3.33106312155723, 7.9264372587204, 7.9264372587204, 10.293275475502, 10.293275475502, 12.43832847476, 14.3289110839036, 15.7610485255718), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9))
mozart.sum = data.frame(TOPP = c(1.77771823108196, 3.17198884487152, 6.66308099031449, 6.66308099031449, 8.32885110378265, 8.32885110378265, 9.99462133646011, 11.6603919863701, 13.326161980629), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MOZART-4", 9))
radm2.sum = data.frame(TOPP = c(2.02861285815918, 5.37516803577029, 7.16689056363599, 7.16689056363599, 7.3984299848477, 7.3984299848477, 8.87811621030172, 10.3578020632267, 9.92148317351077), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RADM2", 9))
racm.sum = data.frame(TOPP = c(2.00438816845417, 5.46073749147613, 7.28098334937261, 7.28098334937261, 8.14680402477582, 8.14680402477582, 9.77616469065348, 11.4055256545544, 8.06489229813425), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9))
racm2.sum = data.frame(TOPP = c(1.81824918091297, 4.28776832090484, 5.71702430645625, 5.71702430645625, 5.685094712036, 5.685094712036, 6.82211367147309, 7.95913251382964, 8.42073631125379), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9))

sum = rbind(cb05.sum, cbm4.sum, cri.sum, mcm1.sum, mcm2.sum, mozart.sum, radm2.sum, racm.sum, racm2.sum)
sum$Mechanism = factor(sum$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05"))

sum.eqn = ddply(sum, .(Mechanism), lm_eqn)

plot = ggplot(sum, aes(x = C, y = TOPP))
plot = plot + geom_point()
plot = plot + geom_smooth(method = "lm", se = FALSE, colour = "black")
plot = plot + geom_text(data = sum.eqn, aes(x = 4.7, y = 16.5, label = V1), size = 3, parse = TRUE)
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

CairoPDF(file = "Alkanes_vs_C_sum.pdf", width = 9, height = 12.7)
print(plot)
dev.off()
