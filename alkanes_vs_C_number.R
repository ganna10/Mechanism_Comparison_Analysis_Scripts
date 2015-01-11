#R script to plot first day TOPP vs carbon number of VOC
# Version 0: Jane Coates 5/1/2015

#maximum TOPP for all alkanes
cb05.data = data.frame(TOPP = c(0.215276703238487, 0.523665487766266, 1.39628636837006, 1.39628636837006, 1.74514518959709, 1.74517107009888, 2.09397292137146, 2.44274997711182, 2.79148983955383), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cbm4.data = data.frame(TOPP = c(0.106132298707962, 0.398037999868393, 1.06159269809723, 1.06159269809723, 1.32720911502838, 1.32720911502838, 1.59285974502563, 1.85856354236603, 2.12430930137634), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
cri.data = data.frame(TOPP = c(0.301193416118622, 0.631050407886505, 1.69316434860229, 1.19700396060944, 2.94067287445068, 2.07821154594421, 3.78355765342712, 4.77706861495972, 5.56714296340942), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm1.data = data.frame(TOPP = c(0.316490352153778, 0.653761506080627, 1.73543918132782, 1.22345685958862, 2.97586512565613, 2.18296837806702, 3.96669507026672, 4.72498798370361, 5.42283296585083), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mcm2.data = data.frame(TOPP = c(0.285142958164215, 0.569781184196472, 1.62394320964813, 1.17484652996063, 2.82610464096069, 2.10391139984131, 3.80315542221069, 4.42267608642578, 4.93928337097168), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
mozart.data = data.frame(TOPP = c(0.287124693393, 0.608259558677673, 1.99833223951006, 2.01251309187801, 1.85423988092579, 1.87466563591183, 1.76157933476884, 1.68114180189661, 1.58341844549298), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("MOZART-4", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
radm2.data = data.frame(TOPP = c(0.342958648488, 1.67143841991616, 1.5708976838747, 1.58514964596511, 1.99016983776016, 2.01276544767989, 1.87340886945383, 1.74628670267072, 2.4296772480011), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8),Mechanism = rep("RADM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm.data = data.frame(TOPP = c(0.329033464193, 1.49119046317878, 1.40149204236976, 1.4142070725483, 2.05412249571477, 2.07744419909902, 1.93360949874895, 1.80240235374113, 1.90360593795776), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))
racm2.data = data.frame(TOPP = c(0.306748777627945, 1.27652395072158, 1.19973819777321, 1.21062281711449, 1.79557052939662, 1.81595673488313, 1.69022647800119, 1.57553434872759, 2.67722725868225), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9), VOC = c("Ethane", "Propane", "Butane", "2-methylpropane", "Pentane", "2-methylbutane", "Hexane", "Heptane", "Octane"))

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
plot = plot + geom_text(data = eqn, aes(x = 4.7, y = 5, label = V1), size = 3, parse = TRUE)
plot = plot + facet_wrap( ~ Mechanism)
plot = plot + ylab("TOPP (molecules(Ox) / molecules(VOC))")
plot = plot + xlab("VOC Carbon Number")
plot = plot + scale_x_continuous(limits = c(2, 8), breaks = seq(2, 8, 1))
plot = plot + scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1))
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
#cb05.sum = data.frame(TOPP = c(1.20671016722918, 2.77097161114216, 7.38695549964905, 7.38721550999356, 9.23052495157695, 9.23052495157695, 11.0736165046692, 12.9159500370983, 14.7577180862427), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CB05", 9))
#cbm4.sum = data.frame(TOPP = c(0.514679666608572, 1.92855651676655, 5.13722425699234, 5.13722425699234, 6.41383069753647, 6.41383069753647, 7.68926954269408, 8.96287225268254, 10.2350900769234), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CBM-IV", 9))
#cri.sum = data.frame(TOPP = c(1.81489448994398, 3.57785767316818, 8.10189753770828, 5.20626991987229, 10.7970547676086, 7.43763422966003, 13.171205163002, 15.1338886973831, 17.2316104471683), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("CRIv2", 9))
#mcm1.sum = data.frame(TOPP = c(1.90840931236744, 3.7286085486412, 8.18034851551057, 5.35230222344398, 10.4700429737568, 8.15893486142159, 12.6083653271198, 14.6831535813028, 16.3026970177889), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.1", 9))
#mcm2.sum = data.frame(TOPP = c(1.73118758946657, 3.33106312155723, 7.9264372587204, 5.02549582719803, 10.293275475502, 7.80973780155182, 12.43832847476, 14.3289110839036, 15.7610485255718), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MCMv3.2", 9))
#mozart.sum = data.frame(TOPP = c(1.77771823108196, 3.17198884487152, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426, 8.32885125279426), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("MOZART-4", 9))
#radm2.sum = data.frame(TOPP = c(2.02861285815918, 5.19599565863609, 5.19599565863609, 5.19599565863609, 7.10249280929565, 7.10249280929565, 7.10249280929565, 7.10249280929565, 9.79746463384187), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RADM2", 9))
#racm.sum = data.frame(TOPP = c(2.00438816845417, 5.27871292829513, 5.27871292829513, 5.27871292829513, 7.8209319114685, 7.8209319114685, 7.8209319114685, 7.8209319114685, 7.96408114440758), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM", 9))
#racm2.sum = data.frame(TOPP = c(1.81824918091297, 5.14532187581062, 5.14532187581062, 5.14532187581062, 6.36730608344079, 6.36730608344079, 6.36730608344079, 6.36730608344079, 8.31547710736311), C = c(2, 3, 4, 4, 5, 5, 6, 7, 8), Mechanism = rep("RACM2", 9))
#
#sum = rbind(cb05.sum, cbm4.sum, cri.sum, mcm1.sum, mcm2.sum, mozart.sum, radm2.sum, racm.sum, racm2.sum)
#sum$Mechanism = factor(sum$Mechanism, levels = c("MCMv3.2", "MCMv3.1", "CRIv2", "RADM2", "RACM", "RACM2", "MOZART-4", "CBM-IV", "CB05"))
#
#sum.eqn = ddply(sum, .(Mechanism), lm_eqn)
#
#plot = ggplot(sum, aes(x = C, y = TOPP))
#plot = plot + geom_point()
#plot = plot + geom_smooth(method = "lm", se = FALSE, colour = "black")
#plot = plot + geom_text(data = sum.eqn, aes(x = 4.7, y = 16.5, label = V1), size = 3, parse = TRUE)
#plot = plot + facet_wrap( ~ Mechanism)
#plot = plot + ylab("TOPP (molecules(Ox) / molecules(VOC))")
#plot = plot + xlab("VOC Carbon Number")
#plot = plot + scale_x_continuous(limits = c(2, 8), breaks = seq(2, 8, 1))
#plot = plot + theme_bw()
#plot = plot + theme(axis.title = element_text(face = "bold"))
#plot = plot + theme(strip.background = element_blank())
#plot = plot + theme(strip.text = element_text(face = "bold"))
#plot = plot + theme(panel.border = element_rect(colour = "black"))
#plot = plot + theme(panel.grid = element_blank())
#
#CairoPDF(file = "Alkanes_vs_C_sum.pdf", width = 9, height = 12.7)
#print(plot)
#dev.off()
