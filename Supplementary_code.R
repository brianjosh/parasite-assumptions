
#Code supporting: "Host and parasite identity interact in scale-dependent
#fashion to determine parasite community structure" (J.I. Brian, D.C. Aldridge) 

#Load packages

library(ggplot2)
library(ggpubr)
library(cowplot)
library(svglite)
library(vegan)
library(MuMIn)
library(dplyr)
library(MASS)

#First, check for independence of Echinoparyphium and Echinostoma in Viviparus 

vivip <- read.csv("viviponly.csv", header=T)
vivippa <- vivip
vivippa$Echinoparyphium <- decostand(vivippa$Echinoparyphium, method = "pa")
vivippa$Echinostoma <- decostand(vivippa$Echinostoma, method = "pa")
count(vivippa, Echinoparyphium, Echinostoma)

#Read in data

paryph <- read.csv("Echinoparyphium.csv", header=T)
stoma <- read.csv("Echinostoma.csv", header=T)
paryphdens <- read.csv("paryphdens.csv", header=T)
stomadens <- read.csv("stomadens.csv", header=T)
paryphdata <- read.csv("paryphmodelling.csv", header=T)
  paryphdata$month <- as.factor(paryphdata$month)
  paryphdata$trem <- as.factor(paryphdata$trem)
stomadata <- read.csv("stomamodelling.csv", header=T)
  stomadata$month <- as.factor(stomadata$month)
  stomadata$trem <- as.factor(stomadata$trem)

#Make versions with just presence-absence data

paryphdatapa <- paryphdata
paryphdatapa$Echinoparyphium <- decostand(paryphdatapa$Echinoparyphium, method = "pa")
stomadatapa <- stomadata
stomadatapa$Echinostoma <- decostand(stomadatapa$Echinostoma, method = "pa")

########################################################################################

#Analysis of prevalence data

#Logistic regression: presence-absence of Echinoparyphium

paryphmodpa <- glm(Echinoparyphium ~ host + month + trem + length + host:month +
                     host:trem + month:trem + host:length + trem:length +
                   month:length, family="binomial", data=paryphdatapa)

options(na.action = "na.fail")
dredge(paryphmodpa)
paryphmodpafinal <- glm(Echinoparyphium ~ host + month + trem + length +
                     host:trem, family="binomial", data=paryphdatapa)
summary(paryphmodpafinal)

#Compare the final full model (paryphmodpafinal) with other models leaving
#out or adding in terms (i.e. goodness of fit)

#No month 
paryphmodpanomonth <- glm(Echinoparyphium ~ host + trem + length +
                          host:trem, family="binomial", data=paryphdatapa)
anova(paryphmodpafinal, paryphmodpanomonth)

#Adding interaction between host and month
paryphmodpaplushostmonth <- glm(Echinoparyphium ~ host + month + trem + length +
                          host:trem + host:month, family="binomial", data=paryphdatapa)
anova(paryphmodpafinal, paryphmodpaplushostmonth)

#No host terms
paryphmodpanohost <- glm(Echinoparyphium ~ month + trem + length,
                          family="binomial", data=paryphdatapa)
anova(paryphmodpafinal, paryphmodpanohost)

#No length term
paryphmodpanolength <- glm(Echinoparyphium ~ host + month + trem +
                             host:trem, family="binomial", data=paryphdatapa)
anova(paryphmodpafinal, paryphmodpanolength)

#No host:castrator interaction
paryphmodpanohostcastrator <- glm(Echinoparyphium ~ host + month + trem + length,
                             family="binomial", data=paryphdatapa)
anova(paryphmodpafinal, paryphmodpanohostcastrator)

species <- split(paryphdatapa, paryphdatapa$host)
anatina <- species$`A. anatina`
vivip <- species$`V. viviparus`
count(anatina, Echinoparyphium, trem)
count(vivip, Echinoparyphium, trem)
min(anatina$length)
min(vivip$length)

#Logistic regression: presence-absence of Echinostoma

stomamodpa <- glm(Echinostoma ~ month + length + trem + length:month + length:trem
                  + month:trem, family="binomial", data=stomadatapa)
options(na.action = "na.fail")
dredge(stomamodpa)
stomamodpafinal <- glm(Echinostoma ~ month + length + trem + length:trem
                  + month:trem, family="binomial", data=stomadatapa)

summary(stomamodpafinal)

#Compare the final full model (stomamodpafinal) with other models leaving
#out or adding in terms (i.e. goodness of fit)

#No month 
stomamodpanomonth <- glm(Echinostoma ~ length + trem + length:trem,
                       family="binomial", data=stomadatapa)
anova(stomamodpafinal, stomamodpanomonth)

#No length:castrator interaction
stomamodpanolengthcastrator <- glm(Echinostoma ~ month + length + trem
                         + month:trem, family="binomial", data=stomadatapa)
anova(stomamodpafinal, stomamodpanolengthcastrator)

#No month:castrator interaction
stomamodpanomonthcastrator <- glm(Echinostoma ~ month + length + trem + length:trem,
                                   family="binomial", data=stomadatapa)
anova(stomamodpafinal, stomamodpanomonthcastrator)

#####################################################################################

#Analysis of intensity data

mean(paryphcount$Echinoparyphium)
var(paryphcount$Echinoparyphium)
#i.e. highly overdispersed - variance 12.5 times larger than mean
mean(stomacount$Echinostoma)
var(stomacount$Echinostoma)
#i.e. highly overdispersed - variance 7.1 times larger than mean

#Echinoparyphium negative binomial 
paryphmodnegbin2 <- glm.nb(Echinoparyphium ~ host + month + trem + length + host:month +
                            host:trem + month:trem + host:length, data=paryphcount)
#Looking at deviance, is a good fit to data 

#Echinoparyphium quasi-poisson
paryphmodqpois2 <- glm(Echinoparyphium ~ host + month + trem + length + host:month +
                        host:trem + month:trem + host:length, family="quasipoisson", data=paryphcount)
#Looking at deviance, is bad fit to data

#Echinostoma negative binomial
stomamodnegbin2 <- glm.nb(Echinostoma ~ month + length + trem
                         + month:trem, data=stomacount)
#Looking at deviance, is good fit to data

#Echinostoma quasipoisson
stomamodqpois2 <- glm(Echinostoma ~ month + length + trem
                     + month:trem, family="quasipoisson", data=stomacount)
#Looking at deviance, is a bad fit to data 

#In both cases, use negative binomial regression

options(na.action = "na.fail")

#Echinoparyphium

paryphmodnegbin2 <- glm.nb(Echinoparyphium ~ host + month + trem + length + host:month +
                     host:trem + month:trem + host:length + trem:length +
                     month:length, data=paryphcount)
dredge(paryphmodnegbin2)
paryphmodfinal2 <- glm.nb(Echinoparyphium ~ host + length + trem +
                           host:trem, data=paryphcount)
summary(paryphmodfinal2)

#Compare the final full model (paryphmodpafinal2) with other models leaving
#out or adding in terms (goodness of fit)

#Addition of month 
paryphmod2month <- glm.nb(Echinoparyphium ~ host + length + trem + month +
                            host:trem, data=paryphcount)
summary(paryphmod2month)
anova(paryphmod2month, paryphmodfinal2)

#No host factors
paryphmod2nohost <- glm.nb(Echinoparyphium ~ length + trem,
                            data=paryphcount)
summary(paryphmod2nohost)
anova(paryphmod2nohost, paryphmodfinal2)

#No length
paryphmod2nolength <- glm.nb(Echinoparyphium ~ host + trem +
                               host:trem, data=paryphcount)
anova(paryphmod2nolength, paryphmodfinal2)

#Include host:length interaction
paryphmod2hostlength <- glm.nb(Echinoparyphium ~ host + trem + length +
                               host:trem + host:length, data=paryphcount)
anova(paryphmod2hostlength, paryphmodfinal2)

#Lack of host:castrator interaction
paryphmod2nohostcastrator <- glm.nb(Echinoparyphium ~ host + trem + length,
                                 data=paryphcount)
anova(paryphmod2nohostcastrator, paryphmodfinal2)

#Echinostoma

stomamodnegbin2 <- glm(Echinostoma ~ month + length + trem + length:month + length:trem
                  + month:trem, data=stomacount)
dredge(stomamodnegbin2)
stomamodfinal2 <- glm.nb(Echinostoma ~ length + month + length:month, data=stomacount)
summary(stomamodfinal2)

#Analyse effect with no interaction
stomamod2nointeraction <- glm.nb(Echinostoma ~ length + month, data=stomacount)
summary(stomamod2nointeraction)
anova(stomamod2nointeraction, stomamodfinal2)

#####################################################################################

#Just exploring some of the interactions further

#See co-occurence of Echinostoma and castrators

month <- split(stomadatapa, stomadatapa$month)
one <- month$`190306`
two <- month$`190403`
three <- month$`190507`
four <- month$`190604`
five <- month$`190625`
six <- month$`190812`
seven <- month$`190905`
eight <- month$`191002`
nine <- month$`191107`
ten <- month$`191202`
eleven <- month$`200120`
twelve <- month$`200224`

count(one, Echinostoma, trem)
count(two, Echinostoma, trem)
count(three, Echinostoma, trem)
count(four, Echinostoma, trem)
count(five, Echinostoma, trem)
count(six, Echinostoma, trem)
count(seven, Echinostoma, trem)
count(eight, Echinostoma, trem)
count(nine, Echinostoma, trem)
count(ten, Echinostoma, trem)
count(eleven, Echinostoma, trem)
count(twelve, Echinostoma, trem)

#12x Mann-Whitney: compare Echinostoma mean vs. pres of trematodes for 
#each of the 12 months separately 

month2 <- split(stomadata, stomadata$month)
one2 <- month2$`190306`
two2 <- month2$`190403`
three2 <- month2$`190507`
four2 <- month2$`190604`
five2 <- month2$`190625`
six2 <- month2$`190812`
seven2 <- month2$`190905`
eight2 <- month2$`191002`
nine2 <- month2$`191107`
ten2 <- month2$`191202`
eleven2 <- month2$`200120`
twelve2 <- month2$`200224`

wilcox.test(one2$Echinostoma ~ one2$trem)
wilcox.test(two2$Echinostoma ~ two2$trem)
wilcox.test(three2$Echinostoma ~ three2$trem)
wilcox.test(four2$Echinostoma ~ four2$trem)
wilcox.test(five2$Echinostoma ~ five2$trem)
wilcox.test(six2$Echinostoma ~ six2$trem)
wilcox.test(seven2$Echinostoma ~ seven2$trem)
wilcox.test(eight2$Echinostoma ~ eight2$trem)
wilcox.test(nine2$Echinostoma ~ nine2$trem)
wilcox.test(ten2$Echinostoma ~ ten2$trem)
wilcox.test(eleven2$Echinostoma ~ eleven2$trem)
wilcox.test(twelve2$Echinostoma ~ twelve2$trem)

#Explore length:castrator interaction - are there more trematodes at increasing lengths?
trempres <- glm(trem ~ length, family="binomial", data=stomadatapa)
summary(trempres)

#See if effect of castrators relates to prevalence of Echinostoma
hypoth <- read.csv("hypothesis.csv", header=T)
hypothmod <- lm(hypoth$diff ~ hypoth$stomaprev)
summary(hypothmod)

#Explore castrator-metacercariae intensity interaction in the two host species
paryphcount <- filter(paryphdata, Echinoparyphium > 0)
stomacount <- filter(stomadata, Echinostoma > 0)
paryphcount2 <- split(paryphcount, paryphcount$host)
anatina3 <- paryphcount2$`A. anatina`
vivip3 <- paryphcount2$`V. viviparus`

anatina3 %>%
  group_by(trem) %>%
  summarize(Mean = mean(Echinoparyphium))
wilcox.test(anatina3$Echinoparyphium ~ anatina3$trem)

vivip3 %>%
  group_by(trem) %>%
  summarize(Mean = mean(Echinoparyphium))
wilcox.test(vivip3$Echinoparyphium ~ vivip3$trem)

# Inspect length effect for Echinostoma on a month-by-month basis

month3 <- split(stomacount, stomacount$month)
one3 <- month3$`190306`
two3 <- month3$`190403`
three3 <- month3$`190507`
four3 <- month3$`190604`
five3 <- month3$`190625`
six3 <- month3$`190812`
seven3 <- month3$`190905`
eight3 <- month3$`191002`
nine3 <- month3$`191107`
ten3 <- month3$`191202`
eleven3 <- month3$`200120`
twelve3 <- month3$`200224`

#####################################################################################

## FIGURES 

#Figure 1 components:

paryphprev <- ggplot(paryph, aes(x=Time, y=Prevalence, color=Host, fill=Host)) +
  geom_point(size=3) + geom_smooth(alpha=0.2) + 
  theme_classic() + xlab("Days since start of 2019") + 
  ylab("E. recurvatum prevalence (%)") + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(legend.position = 'none') 
paryphprev

stomaprev <- ggplot(stoma, aes(x=Time, y=Prevalence)) +
  geom_point(size=3, colour="#E69F00") + geom_smooth(alpha=0.2, colour="#E69F00", fill="#E69F00") + 
  theme_classic() + xlab("Days since start of 2019") + ylab("Echinostoma sp. prevalence (%)")
stomaprev

densityplot <- ggplot(paryphdens, aes(x=Echinoparyphium_sp., color=Host, fill=Host)) +
  geom_density(alpha=0.4) + theme_classic() +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00")) + 
  xlab("E. recurvatum per infected individual") + ylab("Density") +
  theme(legend.position = c(0.8, 0.8))
densityplot

densitystomaplot <- ggplot(stomadens, aes(x=Echinostoma_sp.)) +
  geom_density(alpha=0.4, color="#E69F00", fill="#E69F00") + theme_classic() +
  xlab("Echinostoma sp. per infected individual") + ylab("Density")
densitystomaplot

Figure_1 <- ggarrange(paryphprev, densityplot, stomaprev, 
                      densitystomaplot, nrow=2, ncol=2, labels=c("A", "B", "C", "D"))
Figure_1
ggsave("Figure1.svg", width=160, height=160, units="mm")

# Figure 2 components:

paryphpalengthana <- ggplot(anatina, aes(x=length, y=Echinoparyphium)) +
  geom_count() + geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                             color="gray") + xlab("A. anatina length (mm)") +
  ylab("Probability of E. recurvatum") +
  theme_classic() + theme(legend.position = 'none')
paryphpalengthana

paryphpalengthvivip <- ggplot(vivip, aes(x=length, y=Echinoparyphium)) +
  geom_count() + geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                             color="gray") + xlab("V. viviparus length (mm)") +
  ylab("Probability of E. recurvatum") +
  theme_classic() + theme(legend.position = 'none')
paryphpalengthvivip

negbin <- negative.binomial(theta=1.665, link="log")
paryphlengthana <- ggplot(anatina3, aes(x=length, y=Echinoparyphium)) +
  geom_point() + geom_smooth(method = "glm", method.args = list(family = negbin),  
                             color="gray") + xlab("A. anatina length (mm)") +
  ylab("E. recurvatum intensity") +
  theme_classic() + scale_y_continuous(trans='log10')
paryphlengthana

paryphlengthvivip <- ggplot(vivip3, aes(x=length, y=Echinoparyphium)) +
  geom_point() + geom_smooth(method = "glm", method.args = list(family = negbin),  
                             color="gray") + xlab("V. viviparus length (mm)") +
  ylab("E. recurvatum intensity") +
  theme_classic() + scale_y_continuous(trans='log10')
paryphlengthvivip

plotparyphbar <- ggplot(data=paryphbar, aes(x=Species, y=proportion,
                                            fill=obs.exp)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  ylab("Proportion of E. recurvatum and \ncastrator co-occurences") + theme_classic() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  ylim(0, 0.125) + theme(legend.position = c(0.25, 0.8)) 
plotparyphbar

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-(sd(x)/sqrt(length(x)))
  ymax <- m+(sd(x)/sqrt(length(x)))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
paryphviolin <- ggplot(paryphcount, aes(x=host, y=Echinoparyphium, fill=trem)) +
  geom_violin() + theme_classic() +
  stat_summary(fun.data=data_summary, geom="pointrange", 
               color="black", position = position_dodge(0.9)) +
  scale_fill_manual(name = "Castrators", labels = c("Absent", "Present"), 
                    values=c("#999999", "#E69F00")) +
  ylab("E. recurvatum intensity") + theme(axis.title.x = element_blank()) + 
  scale_y_continuous(trans='log10') + theme(legend.position = c(0.4, 0.8))
paryphviolin

Figure_2 <- ggarrange(paryphpalengthana, paryphpalengthvivip, paryphlengthana,  
                       paryphlengthvivip, plotparyphbar, paryphviolin,
                      nrow=3, ncol=2, labels=c("A", "B", "C", "D", "E", "F"))
Figure_2
ggsave("Figure2.svg", width=160, height=240, units="mm")


#Figure 3 components: 

stomatrem <- split(stomadatapa, stomadatapa$trem)
stomatremabsent <- stomatrem$`0`
stomatrempresent <- stomatrem$`1`

stomapalengthnotrem <- ggplot(stomatremabsent, aes(x=length, y=Echinostoma)) +
  geom_count() + geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                             color="gray") +
  xlab("V. viviparus length (mm)") + ylab("Probability of Echinostoma sp.") +
  theme_classic() + theme(legend.position = 'none')
stomapalengthnotrem

stomapalengthtrem <- ggplot(stomatrempresent, aes(x=length, y=Echinostoma)) +
  geom_count() + geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                             color="gray") +
  xlab("V. viviparus length (mm)") + ylab("Probability of Echinostoma sp.") +
  theme_classic() + theme(legend.position = 'none')
stomapalengthtrem

str(stomadatatrem)
tremplot <- ggplot(stomadatatrem, aes(x=length, y=trem)) +
  geom_count() + scale_size_area() + xlab("V. viviparus length (mm)") +
  ylab("Probability of castrators") + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color="gray") +
  theme_classic() + theme(legend.position = 'none')
tremplot

hypothplot <- ggplot(hypoth, aes(x=stomaprev, y=diff)) +
  geom_point() + geom_smooth(method="lm", color="gray", linetype="dashed") + 
  theme_classic() + xlab("Prevalance of Echinostoma sp.") + 
  ylab("Difference between observed \nand expected co-occurence")
hypothplot

Figure_3 <- ggarrange(stomapalengthnotrem, stomapalengthtrem, hypothplot, 
                      tremplot, nrow=2, ncol=2, labels=c("A", "B", "C", "D"))
Figure_3
ggsave("Figure3.svg", width=160, height=160, units="mm")

#Supplementary figure

month1plot <- ggplot(one3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() + 
  theme(axis.title.x = element_blank()) + ylab("Echinostoma sp. intensity") +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month1plot

month2plot <- ggplot(two3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month2plot

month3plot <- ggplot(three3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month3plot

month4plot <- ggplot(four3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title.x = element_blank()) + ylab("Echinostoma sp. intensity") +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month4plot

month5plot <- ggplot(five3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month5plot

month6plot <- ggplot(six3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month6plot

month7plot <- ggplot(seven3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title.x = element_blank()) + ylab("Echinostoma sp. intensity") +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month7plot

month8plot <- ggplot(eight3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month8plot

month9plot <- ggplot(nine3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  theme(axis.title = element_blank()) + scale_y_continuous(trans='log10') + xlim(5.5, 15)
month9plot

month10plot <- ggplot(ten3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  xlab("V. viviparus length (mm)") + ylab("Echinostoma sp. intensity") +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month10plot

month11plot <- ggplot(eleven3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  xlab("V. viviparus length (mm)") + theme(axis.title.y = element_blank()) +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month11plot

month12plot <- ggplot(twelve3, aes(x=length, y=Echinostoma)) +
  geom_point() + geom_smooth(method="lm", color="gray") + theme_classic() +
  xlab("V. viviparus length (mm)") + theme(axis.title.y = element_blank()) +
  scale_y_continuous(trans='log10') + xlim(5.5, 15)
month12plot

Supp_Fig <- ggarrange(month1plot, month2plot, month3plot, month4plot, 
                      month5plot, month6plot, month7plot, month8plot, 
                      month9plot, month10plot, month11plot, month12plot, 
                      nrow=4, ncol=3, labels=c("A", "B", "C", "D", "E", "F",
                                               "G", "H", "I", "J", "K", "L"))
Supp_Fig
ggsave("Supplementary_Fig.svg", width=233, height=300, units="mm")






