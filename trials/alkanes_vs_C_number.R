#R script to plot first day TOPP vs carbon number of VOC
# Version 0: Jane Coates 5/1/2015

#first day TOPPs for all alkanes
cb05.data = data.frame(TOPP = c(0.0845172479748726, 0.3641597032547, 0.952807605266571, 0.95486788938202, 1.16583390916096, 1.16891467571259, 1.37504136562347, 1.57824212238862, 1.7769056558609), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.data = data.frame(TOPP = c(0.0869439467787743, 0.321530640125275, 0.840508162975311, 0.842413246631622, 1.0273551940918, 1.03020417690277, 1.21066904067993, 1.38842662532035, 1.56198441982269), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cri.data = data.frame(TOPP = c(0.114152513444424, 0.354555159807205, 0.972761750221252, 0.967499315738678, 1.92756509780884, 1.57766127586365, 2.27819347381592, 2.37771141645028, 2.62646389007568), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.data = data.frame(TOPP = c(0.12388314306736, 0.381832748651505, 1.03673362731934, 1.00873160362244, 1.7392680644989, 1.60785245895386, 2.22197437286377, 2.43809782036837, 2.71707463264465), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.data = data.frame(TOPP = c(0.114175207912922, 0.346426695585251, 0.969126343727112, 0.97936624288559, 1.65823602676392, 1.56109631061554, 2.15132093429565, 2.35366020437529, 2.59780859947205), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mozart.data = data.frame(TOPP = c(0.13925464451313, 0.426944702863693, 1.37005066871643, 1.37005066871643, 1.37005066871643, 1.37005066871643, 1.37005066871643, 1.37005066871643, 1.37005066871643), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
radm2.data = data.frame(TOPP = c(0.128240107806014, 0.907400190830231, 0.907400190830231, 0.907400190830231, 1.48803865909576, 1.48803865909576, 1.48803865909576, 1.48803865909576, 2.61359875085371), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm.data = data.frame(TOPP = c(0.131782606244087, 0.732421875, 0.732421875, 0.732421875, 1.40603709220886, 1.40603709220886, 1.40603709220886, 1.40603709220886, 2.13747172994067), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm2.data = data.frame(TOPP = c(0.11060257256031, 0.830682754516602, 0.830682754516602, 0.830682754516602, 1.34074258804321, 1.34074258804321, 1.34074258804321, 1.34074258804321, 2.46226828920276), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))

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

#cumulative data
cb05.sum = data.frame(TOPP = c(1.20671016722918, 2.77097161114216, 7.38695549964905, 7.38721550999356, 9.23052495157695, 9.23052495157695, 11.0736165046692, 12.9159500370983, 14.7577180862427), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9))
cbm4.sum = data.frame(TOPP = c(0.514679666608572, 1.92855651676655, 5.13722425699234, 5.13722425699234, 6.41383069753647, 6.41383069753647, 7.68926954269408, 8.96287225268254, 10.2350900769234), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9))
cri.sum = data.frame(TOPP = c(1.81489448994398, 3.57785767316818, 8.10189753770828, 5.20626991987229, 10.7970547676086, 7.43763422966003, 13.171205163002, 15.1338886973831, 17.2316104471683), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9))
mcm1.sum = data.frame(TOPP = c(1.90840931236744, 3.7286085486412, 8.18034851551057, 5.35230222344398, 10.4700429737568, 8.15893486142159, 12.6083653271198, 14.6831535813028, 16.3026970177889), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9))
mcm2.sum = data.frame(TOPP = c(1.73118758946657, 3.33106312155723, 7.9264372587204, 5.02549582719803, 10.293275475502, 7.80973780155182, 12.43832847476, 14.3289110839036, 15.7610485255718), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9))
mozart.sum = data.frame(TOPP = c(1.77771823108196, 3.17198884487152, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MOZART-4", 9))
radm2.sum = data.frame(TOPP = c(2.02861285815918, 5.19599565863609, 5.19599565863609, 5.19599565863609, 7.10249280929565, 7.10249280929565, 7.10249280929565, 7.10249280929565, 9.79746463384187), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RADM2", 9))
racm.sum = data.frame(TOPP = c(2.00438816845417, 5.27871292829513, 5.27871292829513, 5.27871292829513, 7.8209319114685, 7.8209319114685, 7.8209319114685, 7.8209319114685, 7.96408114440758), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9))
racm2.sum = data.frame(TOPP = c(1.81824918091297, 5.14532187581062, 5.14532187581062, 5.14532187581062, 6.36730608344079, 6.36730608344079, 6.36730608344079, 6.36730608344079, 8.31547710736311), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9))

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
