## eidolon helvum: paramyxovirus analyses
## maya juman
## updated 02/09/24

## clean up space
rm(list=ls()) 
graphics.off()

## set up

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(readxl)
library(GGally)
library(dplyr)
library(factoextra)
library(lme4)
library(lmerTest)
library(mgcv)
library(Hmisc)
library(corrplot)

## set up data

setwd("~/Documents/PhD/eidolon/paramyxo submission")
MFI.data <- read_excel("Supplementary Data.xlsx",sheet="Data", guess_max = 5000) %>% 
  mutate(Site = factor(Site)) %>%
  dplyr::rename(Weight = `Weight (g)`,
                Forearm = `Forearm (mm)`,
                Age = `Age Category`,
                Repro = `Reproductive Status`)

#defining body condition as weight/forearm length (Plowright et al. 2008)
MFI.data$BodyCondition <- MFI.data$Weight/MFI.data$Forearm

bead.col <- which(str_starts(names(MFI.data),"Analyte"))
bead.names <- names(MFI.data)[bead.col]
virus.names <- str_split_fixed(bead.names," ",3)[,3]

MFI.data.long <- MFI.data %>% 
  pivot_longer(starts_with("Analyte"), names_to="Analyte", values_to="MFI") %>%
  mutate(Analyte = str_split_fixed(Analyte," ",3)[,3],
         log.MFI = log(MFI))

bead.key <- read_excel("Supplementary Data.xlsx",sheet="Key")
controls <- read_excel("Supplementary Data.xlsx",sheet="Controls",na = "NaN")

mAb.controls <- controls %>% 
  filter(str_starts(Sample,"mAb")) %>%
  separate(Sample, into=c("mAb","Dilution"),sep=" 1:", extra="drop") %>%
  separate(Dilution, into=c("Dilution"),sep=" ", extra="drop") %>%
  mutate(Dilution = as.integer(Dilution))

mAb.controls.long <- mAb.controls %>%
  pivot_longer(starts_with("Analyte"), names_to="Analyte", values_to="MFI") %>%
  mutate(Analyte = str_split_fixed(Analyte," ",3)[,3],
         log.MFI = log(MFI),
         Control = sapply(Analyte, function(x){with(bead.key, Control_mAB[which(Abbreviation==x)])}),
         Match = (mAb==Control)
  )

#generate corrected log MFI values by subtracting background (mock) MFI

MFI.mock.corrected <- MFI.data.long %>% 
  group_by(Sample) %>%
  mutate(log.MFI.mc = log.MFI - log.MFI[Analyte=="MOCK"]) %>%
  ungroup() %>%
  filter(Analyte != "MOCK")

MFI.mock.corrected$Analyte <- factor(MFI.mock.corrected$Analyte, 
                                     levels = c("GHV","NIV","HEV","CEDV","MOJV","MENV"))

MFI.mock.corrected$Site <- factor(MFI.mock.corrected$Site, 
                                     levels = c("Accra urban roost","Akosombo rural roost",
                                                "Kumasi urban roost","Captive colony"))

MFI.mock.corrected.wide <- MFI.mock.corrected %>% dplyr::select(-MFI, -log.MFI) %>% 
  pivot_wider(names_from = Analyte, values_from = log.MFI.mc)

#monoclonal antibodies as positive controls

bead.key %>% filter(!is.na(Control_mAB)) %>% dplyr::select(c(4,8))

#mock correction for mAb
mAb.controls.mc <- mAb.controls.long %>% 
  group_by(`Luminex Date`, mAb, Dilution, Analyte) %>%
  dplyr::summarise(log.MFI = mean(log.MFI), .groups = "keep") %>%
  group_by(`Luminex Date`, mAb, Dilution) %>%
  mutate(log.MFI.mc = log.MFI - log.MFI[Analyte=="MOCK"],
         Control = sapply(Analyte, function(x){with(bead.key, Control_mAB[which(Abbreviation==x)])}),
         Match = (mAb==Control)
  )

#monoclonal antibody cutoffs for all possible viruses
#Inflection point: x = e/b, y = (d-c)/2

bead.key %>% filter(!is.na(Control_mAB)) %>% dplyr::select(c(4,8))

mAb.controls.mc$`Luminex Date` <- as.factor(mAb.controls.mc$`Luminex Date`)

#NIV
sig.NIV.control <- nls(log.MFI.mc ~ c+(d-c)/(1+(exp(b*log10(Dilution)-e))), 
                       data = mAb.controls.mc %>% 
                         filter(Analyte=="NIV" & Match & `Luminex Date` != "2021-03-10"),  
                       start = list(b=1, c=0, d=2, e=5))

NIV.ctr.par <- sig.NIV.control$m$getPars()
NIVcut <- (NIV.ctr.par['d']-NIV.ctr.par['c'])/2

NIVsero <- MFI.mock.corrected %>% 
  filter(Analyte=="NIV") %>%
  mutate(Positive = (log.MFI.mc>NIVcut)) %>%
  group_by(Site, `Sampling month`, Analyte) %>% 
  dplyr::summarise(Prevalence = mean(as.integer(Positive)))

#HEV

sig.hev.control <- nls(log.MFI.mc ~ c+(d-c)/(1+(exp(b*log10(Dilution)-e))), 
                       data = mAb.controls.mc %>% 
                         filter(Analyte=="HEV" & Match & `Luminex Date` != "2021-03-10"),  
                       start = list(b=1, c=0, d=2, e=5))

hev.ctr.par <- sig.hev.control$m$getPars()
HEVcut <- (hev.ctr.par['d']-hev.ctr.par['c'])/2 

HEVsero <- MFI.mock.corrected %>% 
  filter(Analyte=="HEV") %>%
  mutate(Positive = (log.MFI.mc>HEVcut)) %>%
  group_by(Site, `Sampling month`, Analyte) %>% 
  dplyr::summarise(Prevalence = mean(as.integer(Positive)))

#using an average of NIV and HEV for GHV
GHVsero <- MFI.mock.corrected %>% 
  filter(Analyte=="GHV") %>%
  mutate(Positive = (log.MFI.mc>((NIVcut+HEVcut)/2))) %>%
  group_by(Site, `Sampling month`, Analyte) %>% 
  dplyr::summarise(Prevalence = mean(as.integer(Positive)))

#CEDV

sig.cedv.control <- nls(log.MFI.mc ~ c+(d-c)/(1+(exp(b*log10(Dilution)-e))), 
                        data = mAb.controls.mc %>% filter(Analyte=="CEDV" & Match),  
                        start = list(b=1, c=0, d=2, e=5))

cedv.ctr.par <- sig.cedv.control$m$getPars()
CEDVcut <- (cedv.ctr.par['d']-cedv.ctr.par['c'])/2 #2.137574 

CEDVsero <- MFI.mock.corrected %>% 
  filter(Analyte=="CEDV") %>%
  mutate(Positive = (log.MFI.mc>CEDVcut)) %>%
  group_by(Site, `Sampling month`, Analyte) %>% 
  dplyr::summarise(Prevalence = mean(as.integer(Positive)))

#pulling seroprevalences together
seroprevalences <- bind_rows(NIVsero, HEVsero, GHVsero)

#alternative cutoffs for sensitivity analyses

paramyxo_cutoffs2 <- data.frame(matrix(ncol=4, nrow=5, 
                                      dimnames = list(c(1:5),c("NIV","HEV","MOJV","GHV"))))
paramyxo_cutoffs2[3,c(1:2)] <- c(NIVcut, HEVcut)
paramyxo_cutoffs2[3,c(3,4)] <- mean(c(NIVcut, HEVcut))
paramyxo_cutoffs2[c(1,2,4,5),] <- c(0.2,1,2.5,3.5)

###plots of mAb curves
nivmab <- mAb.controls.mc %>% 
  filter(Analyte=="NIV" & Match & `Luminex Date` != "2021-03-10") %>% 
  ggplot(aes(x=log10(Dilution), y=log.MFI.mc)) +
  geom_line(aes(col=factor(`Luminex Date`))) +
  geom_hline(yintercept=NIVcut, linetype="solid", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[1,1], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[2,1], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[4,1], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[5,1], linetype="dashed", col="red") + 
  theme(legend.position = "none") + ylab("log corrected MFI value")+
  ylim(-0.3,4.4)

hevmab <- mAb.controls.mc %>% 
  filter(Analyte=="HEV" & Match & `Luminex Date` != "2021-03-10") %>% 
  ggplot(aes(x=log10(Dilution), y=log.MFI.mc)) +
  geom_line(aes(col=factor(`Luminex Date`))) +
  geom_hline(yintercept=HEVcut, linetype="solid", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[1,2], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[2,2], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[4,2], linetype="dashed", col="red") + 
  geom_hline(yintercept=paramyxo_cutoffs2[5,2], linetype="dashed", col="red") + 
  ylab("log corrected MFI value") + theme(legend.position = "none") +
  ylim(-0.3,4.4)

ggarrange(nivmab, hevmab, labels = c("A", "B"), align="v",ncol=2)
#ggsave("suppFig1.png", width = 250, height = 100, units= "mm", bg="white")

rownames(paramyxo_cutoffs2) <- c("Cutoff 1", "Cutoff 2", "Cutoff 3", "Cutoff 4", "Cutoff 5")
paramyxo_sero <- data.frame()
paramyxo_sero_sex <- data.frame()

for (i in 1:5) {
  for (j in 1:4) {
    x <- MFI.mock.corrected %>% 
      filter(Analyte==colnames(paramyxo_cutoffs2[j])) %>%
      mutate(Positive = (log.MFI.mc>paramyxo_cutoffs2[i,j])) %>%
      group_by(Site, `Sampling month`, Analyte) %>% 
      dplyr::summarise(Prevalence = mean(as.integer(Positive)),
                Positives = sum(as.integer(Positive)),
                n=n())
    x$cutoff <- rownames(paramyxo_cutoffs2[i,])
    paramyxo_sero <- rbind(paramyxo_sero, x)
    
    y <- MFI.mock.corrected %>% drop_na(Sex) %>%
      filter(Analyte==colnames(paramyxo_cutoffs2[j])) %>%
      mutate(Positive = (log.MFI.mc>paramyxo_cutoffs2[i,j])) %>%
      group_by(Site, `Sampling month`, Analyte, Sex) %>% 
      dplyr::summarise(Prevalence = mean(as.integer(Positive)),
                Positives = sum(as.integer(Positive)),
                n=n())
    y$cutoff <- rownames(paramyxo_cutoffs2[i,])
    paramyxo_sero_sex <- rbind(paramyxo_sero_sex, y)
  }
}

paramyxo_sero$`Sampling month` <- as.Date(paramyxo_sero$`Sampling month`)
paramyxo_sero_sex$`Sampling month` <- as.Date(paramyxo_sero_sex$`Sampling month`)

conf <- binconf(paramyxo_sero$Positives, paramyxo_sero$n, alpha=0.05)[,2:3]
paramyxo_sero <- cbind(paramyxo_sero, conf)
confsex <- binconf(paramyxo_sero_sex$Positives, paramyxo_sero_sex$n, alpha=0.05)[,2:3]
paramyxo_sero_sex <- cbind(paramyxo_sero_sex, confsex)

paramyxo_sero %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=cutoff)) +
  geom_line() + facet_wrap(vars(Analyte, Site))

####summary table (supplementary table S2)

MFI.mock.corrected.wide %>% group_by(Site, `Sampling month`) %>%
  summarise(n=n(),
            females=sum(Sex == "F", na.rm=TRUE),
            propf=sum(Sex == "F", na.rm=TRUE)/n()*100,
            males=sum(Sex == "M", na.rm=TRUE), 
            propm=sum(Sex == "M", na.rm=TRUE)/n()*100,
            adults=sum(Age == "A"),
            propad=sum(Age == "A")/n()*100,
            subadults=sum(Age == "SA"),
            propsub=sum(Age == "SA")/n()*100,
            juveniles=sum(Age == "J"),
            propj=sum(Age == "J")/n()*100)

##supplementary table S4
round((sum(MFI.mock.corrected.wide$NIV > paramyxo_cutoffs2[5,1])/2090)*100,2)

### Multivariate analyses

###all MFI values

levels(MFI.mock.corrected$Analyte) <- c("GhV","NiV","HeV","CedV","MojV","MenV")

all <- ggplot(MFI.mock.corrected) + 
  geom_density(aes(x=log.MFI.mc, fill=Site, col=Site), alpha=0.4) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(vars(Analyte),scales = "free") + 
  xlab("Log mock-corrected MFI values") +
  ylab("Frequency") +
  theme(plot.margin=unit(c(0.2,-1,0.2,0.2), "cm")) +
  theme_minimal()

#correlation plot

MFI.mock.corrected.wide <- MFI.mock.corrected.wide %>% 
  dplyr::select(Site, `Sampling month`, Sample,
                Age, `Age in Years`, Sex, Repro, Weight, Forearm, `Bat ID`, 
                BodyCondition, GHV, HEV, NIV, CEDV, MOJV, MENV)

MFI.mock.corrected.wide.caps <- MFI.mock.corrected.wide %>% 
  rename("GhV" = "GHV", "HeV" = "HEV", "NiV" = "NIV",
         "CedV" = "CEDV", "MojV" = "MOJV", "MenV" = "MENV")

library(ggcorrplot)

cor <- corrplot(cor(MFI.mock.corrected.wide.caps[,12:17]))
ggcor <- ggcorrplot(cor(MFI.mock.corrected.wide.caps[,17:12]))

ggarrange(all, ggcor, labels = c("A", "B"), widths = c(1, 0.7), ncol=2)
#ggsave("newFig1.png", width = 320, height = 100, units= "mm", bg="white")

#pca

paramyxo.names <- bead.key %>% filter(Family=="Paramyxoviridae") %>% pull(Abbreviation)
paramyxo.columns <- which(names(MFI.mock.corrected.wide) %in% paramyxo.names)

#paramyxos only

paramyxo.PCA <- MFI.mock.corrected.wide.caps %>% 
  dplyr::select(all_of(paramyxo.columns)) %>% prcomp(scale. = TRUE)
fviz_pca_biplot(paramyxo.PCA,habillage=MFI.mock.corrected.wide.caps$Site, 
                addEllipses=TRUE, ellipse.level=0.95, label="var")
#ggsave("PCAparamyxo.jpg", width=10,height=6)

## WILD ROOSTS

######seroprevalence patterns

#GHV

GHVpatterns <- paramyxo_sero_sex %>% filter(Analyte == "GHV", cutoff == 3) %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_point(aes(size=n)) + scale_size_continuous(range = c(0.1, 3)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Sex), 
              alpha=0.1,
              color=NA) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("GhV seroprevalence")

#sensitivity analyses

GHVsens <- paramyxo_sero_sex %>% filter(Analyte == "GHV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  geom_line() + 
  facet_wrap(vars(Site, cutoff), labeller = label_wrap_gen(multi_line=FALSE)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("GhV seroprevalence")

#NIV

NIVpatterns <- paramyxo_sero_sex %>% filter(Analyte == "NIV", cutoff == 3) %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_point(aes(size=n)) + scale_size_continuous(range = c(0.1, 3)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Sex), 
              alpha=0.1,
              color=NA) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("NiV seroprevalence")

#sensitivity analyses

NIVsens <- paramyxo_sero_sex %>% filter(Analyte == "NIV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  geom_line() + facet_wrap(vars(Site, cutoff)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("NiV seroprevalence")

#HEV

HEVpatterns <- paramyxo_sero_sex %>% filter(Analyte == "HEV", cutoff == 3) %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_point(aes(size=n)) + scale_size_continuous(range = c(0.1, 3)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Sex), 
              alpha=0.1,
              color=NA) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("HeV seroprevalence")


#sensitivity analyses

HEVsens <- paramyxo_sero_sex %>% filter(Analyte == "HEV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  geom_line() + facet_wrap(vars(Site, cutoff)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("HeV seroprevalence")


###MOJV

MOJVpatterns <- paramyxo_sero_sex %>% filter(Analyte == "MOJV", cutoff == 3) %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_point(aes(size=n)) + scale_size_continuous(range = c(0.1, 3)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Sex), 
              alpha=0.1,
              color=NA) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_sero_sex, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("MojV seroprevalence")

#sensitivity analyses

MOJVsens <- paramyxo_sero_sex %>% filter(Analyte == "MOJV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex)) +
  geom_line() + facet_wrap(vars(Site, cutoff)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("MojV seroprevalence")


fig2 <- ggarrange(GHVpatterns, NIVpatterns, labels = c("A", "B"), 
                  widths = c(0.7, 1), align="v",ncol=1)
#ggsave("Fig2.png", width = 150, height = 234, units= "mm", bg="white")

ggarrange(GHVsens, NIVsens, HEVsens, MOJVsens, labels = c("A", "B", "C", "D"), 
                      align="v",ncol=2,nrow=2)
#ggsave("suppFig2.png", width = 400, height = 300, units= "mm", bg="white")

suppfig4 <- ggarrange(HEVpatterns, MOJVpatterns, labels = c("A", "B"),
                      widths = c(0.7, 1), align="v",ncol=1)
#ggsave("FigS4.png", width = 150, height = 234, units= "mm", bg="white")

######linear models

#creating a continuous time variable + time squared
MFI.mock.corrected.wide$date <- as.POSIXct(MFI.mock.corrected.wide$`Sampling month`)
MFI.mock.corrected.wide$time <- as.numeric(MFI.mock.corrected.wide$date)
MFI.mock.corrected.wide$time <- (MFI.mock.corrected.wide$time-min(MFI.mock.corrected.wide$time))/(max(MFI.mock.corrected.wide$time)-min(MFI.mock.corrected.wide$time))
MFI.mock.corrected.wide$timesq <- (MFI.mock.corrected.wide$time*MFI.mock.corrected.wide$time)

###combining reproduction + sex variables
MFI.mock.corrected.wide$ReproSex <- MFI.mock.corrected.wide$Repro
MFI.mock.corrected.wide$ReproSex[which(MFI.mock.corrected.wide$Sex == "M")] <- "Male"
MFI.mock.corrected.wide$ReproSex[which(is.na(MFI.mock.corrected.wide$ReproSex))] <- "Non-preg female"
MFI.mock.corrected.wide$ReproSex <- factor(MFI.mock.corrected.wide$ReproSex, 
                                           levels = c("Non-preg female","Pregnant","Lactating","Male"))

#releveling age categories if necessary
#MFI.mock.corrected.wide$Age <- factor(MFI.mock.corrected.wide$Age, levels = c("SA","J","A"))
 
####wild roosts only: comparing models with and without time^2

#GHV

timsqghv <- (glm(GHV ~ Age + Sex + BodyCondition + Site + time + timesq, 
            data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))

timghv <- (glm(GHV ~ Age + Sex + BodyCondition + Site + time, 
                 data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))
anova(timghv, timsqghv, test = "Chisq")

summary(timsqghv)

#NIV
timsqniv <- (glm(NIV ~ Age + Sex + BodyCondition + Site + time + timesq, 
                 data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))

timniv <- (glm(NIV ~ Age + Sex + BodyCondition + Site + time, 
               data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))
anova(timniv, timsqniv, test = "Chisq")

summary(timsqniv)

#HEV

timsqhev <- (glm(HEV ~ Age + Sex + BodyCondition + Site + time + timesq, 
                 data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))

timhev <- (glm(HEV ~ Age + Sex + BodyCondition + Site + time, 
               data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))
anova(timhev, timsqhev, test = "Chisq")

summary(timsqhev)

#MOJV

timsqmojv <- (glm(MOJV ~ Age + Sex + BodyCondition + Site + time + timesq, 
                 data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))

timmojv <- (glm(MOJV ~ Age + Sex + BodyCondition + Site + time, 
               data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony")))
anova(timmojv, timsqmojv, test = "Chisq")

summary(timsqmojv)


#LMEMs of captive col
summary(lmer(GHV ~ Age + BodyCondition + ReproSex + (1|`Bat ID`), 
             data=MFI.mock.corrected.wide %>% filter(Site == "Captive colony")))

summary(lmer(NIV ~ Age + BodyCondition + ReproSex + (1|`Bat ID`), 
             data=MFI.mock.corrected.wide %>% filter(Site == "Captive colony")))

summary(lmer(HEV ~ Age + BodyCondition + ReproSex + (1|`Bat ID`), 
             data=MFI.mock.corrected.wide %>% filter(Site == "Captive colony")))

summary(lmer(MOJV ~ Age + BodyCondition + ReproSex + (1|`Bat ID`), 
             data=MFI.mock.corrected.wide %>% filter(Site == "Captive colony")))

######CAPTIVE COLONY: short term

#length(unique(AZ$`Bat ID`)) #179 unique bats
AZall <- MFI.mock.corrected.wide %>% filter(Site=="Captive colony")

###GHV

##trajectories for five sampling points
#by sex
GHVfem <- ggplot(AZall %>% filter(Age == "A", Sex == "F"), 
                 aes(x=`Sampling month`, y=GHV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="red") + geom_point(colour="red") + 
  ylim(-0.5,4.1) +
  theme(legend.position = "none") + ylab("GhV log MFI values") + theme_minimal()

GHVmal <- ggplot(AZall %>% filter(Age == "A", Sex == "M"), 
                 aes(x=`Sampling month`, y=GHV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="blue") + geom_point(colour="blue") + 
  ylim(-0.5,4.1) +
  theme(legend.position = "none") + ylab("GhV log MFI values") + theme_minimal()

##GAMMs

AZall <- MFI.mock.corrected.wide %>% filter(Site=="Captive colony")
AZall$ID <- AZall$`Bat ID`
AZall$date <- as.numeric(as.Date(AZall$`Sampling month`))

ghv2<-((gamm(GHV~s(date, bs="cr", k=5, by=Sex)+Sex,
             data=AZall %>% filter(Age == "A"),random=list(ID=~1))))
summary(ghv2$gam)

#gam.check(ghv2$gam)
#k.check(gamm1$gam) #edf < k', and non-significant resids, so this is OK

newdata <- data.frame(date=rep(17928:18536, 2))
newdata$Sex[1:609] <- "M"
newdata$Sex[609:1218] <- "F"

predGHV <- bind_cols(newdata, as.data.frame(predict(ghv2$gam, newdata=newdata, type="link", se.fit=TRUE)))
predGHV$date <- as.POSIXct(as.Date(predGHV$date, origin = "1970-01-01"))
predGHV$upr <- predGHV$fit + (1.96 * predGHV$se.fit)
predGHV$lwr <- predGHV$fit - (1.96 * predGHV$se.fit)

GHVgamplot <- AZall %>% filter(!is.na(Sex), Age == "A") %>% ggplot(aes(x=`Sampling month`)) +
  geom_line(data=predGHV, aes(x=date, y=fit, colour=Sex)) +
  geom_ribbon(data=predGHV, aes(x=date, ymin=lwr, ymax=upr, fill=Sex), alpha=0.2) +
  ylim(-0.5,4.1) + theme_bw() + ylab("GhV log MFI values") +
  scale_color_manual(values=c("red","blue")) + theme_minimal() +
  geom_segment(aes(x = as.POSIXct("2020-07-01"), y = 2.3, 
                   xend = as.POSIXct("2020-07-01"), yend = 1.5),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = as.POSIXct("2019-07-01"), y = 2.3, 
                   xend = as.POSIXct("2019-07-01"), yend = 1.5),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm")))

###NIV

##trajectories for five sampling points

#by sex
NIVfem <- ggplot(AZall %>% filter(Age == "A", Sex == "F"), 
                 aes(x=`Sampling month`, y=NIV, group=`Bat ID`)) +
    geom_line(linewidth=0.1, colour="red") + geom_point(colour="red") + 
    ylim(-0.2,4.3) +
    theme(legend.position = "none") + ylab("NiV log MFI values") + theme_minimal()
  
NIVmal <- ggplot(AZall %>% filter(Age == "A", Sex == "M"), 
                 aes(x=`Sampling month`, y=NIV, group=`Bat ID`)) +
    geom_line(linewidth=0.1, colour="blue") + geom_point(colour="blue") + 
    ylim(-0.2,4.3) +
    theme(legend.position = "none") + ylab("NiV log MFI values") + theme_minimal()

##GAMMs

niv2<-((gamm(NIV~s(date, bs="cr", k=5, by=Sex)+Sex,
              data=AZall %>% filter(Age == "A"),random=list(ID=~1))))

summary(niv2$gam)
#gam.check(niv2$gam)
#k.check(niv2$gam) #edf < k', and non-significant resids, so this is OK

predNIV <- bind_cols(newdata, as.data.frame(predict(niv2$gam, newdata=newdata, type="link", se.fit=TRUE)))
predNIV$date <- as.POSIXct(as.Date(predNIV$date, origin = "1970-01-01"))
predNIV$upr <- predNIV$fit + (1.96 * predNIV$se.fit)
predNIV$lwr <- predNIV$fit - (1.96 * predNIV$se.fit)

NIVgamplot <- AZall %>% filter(!is.na(Sex), Age == "A") %>% ggplot(aes(x=`Sampling month`)) +
  geom_line(data=predNIV, aes(x=date, y=fit, colour=Sex)) +
  geom_ribbon(data=predNIV, aes(x=date, ymin=lwr, ymax=upr, fill=Sex), alpha=0.2) +
  ylim(-0.2,4.3) + theme_bw() + ylab("NiV log MFI values") +
  scale_color_manual(values=c("red","blue")) + theme_minimal() +
  geom_segment(aes(x = as.POSIXct("2020-07-01"), y = 2.8, 
                   xend = as.POSIXct("2020-07-01"), yend = 2),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = as.POSIXct("2019-07-01"), y = 2.8, 
                   xend = as.POSIXct("2019-07-01"), yend = 2),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm")))

ggarrange(GHVfem, GHVmal, GHVgamplot, NIVfem, NIVmal, NIVgamplot, 
          ncol=3, nrow=2, labels = c("A","B","C","D","E","F"))
#ggsave("Fig2.png", width = 300, height = 180, units= "mm", bg="white")

###HEV

##trajectories for five sampling points
#by sex
HEVfem <- ggplot(AZall %>% filter(Age == "A", Sex == "F"), 
                 aes(x=`Sampling month`, y=HEV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="red") + geom_point(colour="red") + 
  ylim(-0.9,3.2) +
  theme(legend.position = "none") + ylab("HeV log MFI values") + theme_minimal()

HEVmal <- ggplot(AZall %>% filter(Age == "A", Sex == "M"), 
                 aes(x=`Sampling month`, y=HEV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="blue") + geom_point(colour="blue") + 
  ylim(-0.9,3.2) +
  theme(legend.position = "none") + ylab("HeV log MFI values") + theme_minimal()

##GAMMs

hev2<-((gamm(HEV~s(date, bs="cr", k=5, by=Sex)+Sex,
             data=AZall %>% filter(Age == "A"),random=list(ID=~1))))

summary(hev2$gam)
#gam.check(niv2$gam)
#k.check(niv2$gam) #edf < k', and non-significant resids, so this is OK

predHEV <- bind_cols(newdata, as.data.frame(predict(hev2$gam, newdata=newdata, type="link", se.fit=TRUE)))
predHEV$date <- as.POSIXct(as.Date(predHEV$date, origin = "1970-01-01"))
predHEV$upr <- predHEV$fit + (1.96 * predHEV$se.fit)
predHEV$lwr <- predHEV$fit - (1.96 * predHEV$se.fit)

HEVgamplot <- AZall %>% filter(!is.na(Sex), Age == "A") %>% ggplot(aes(x=`Sampling month`)) +
  geom_line(data=predHEV, aes(x=date, y=fit, colour=Sex)) +
  geom_ribbon(data=predHEV, aes(x=date, ymin=lwr, ymax=upr, fill=Sex), alpha=0.2) +
  ylim(-0.9,3.2) + theme_bw() + ylab("HeV log MFI values") +
  scale_color_manual(values=c("red","blue")) + theme_minimal() +
  geom_segment(aes(x = as.POSIXct("2020-07-01"), y = 1.9, 
                   xend = as.POSIXct("2020-07-01"), yend = 1.1),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = as.POSIXct("2019-07-01"), y = 1.9, 
                   xend = as.POSIXct("2019-07-01"), yend = 1.1),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm")))


#ggarrange(allHEV, HEVgamplot, labels = c("A","B"))
#ggsave("HEV short term.png", width = 200, height = 70, units= "mm", bg="white")

####MOJV

#by sex
MOJVfem <- ggplot(AZall %>% filter(Age == "A", Sex == "F"), aes(x=`Sampling month`, y=MOJV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="red") + geom_point(colour="red") + 
  ylim(-0.1,3) +
  theme(legend.position = "none") + ylab("MojV log MFI values") + theme_minimal()

MOJVmal <- ggplot(AZall %>% filter(Age == "A", Sex == "M"), aes(x=`Sampling month`, y=MOJV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="blue") + geom_point(colour="blue") + 
  ylim(-0.1,3) +
  theme(legend.position = "none") + ylab("MojV log MFI values") + theme_minimal()

##GAMMs

mojv2<-((gamm(MOJV~s(date, bs="cr", k=5, by=Sex)+Sex,
             data=AZall %>% filter(Age == "A"),random=list(ID=~1))))

summary(mojv2$gam)
#gam.check(niv2$gam)
#k.check(niv2$gam) #edf < k', and non-significant resids, so this is OK

predMOJV <- bind_cols(newdata, as.data.frame(predict(mojv2$gam, newdata=newdata, type="link", se.fit=TRUE)))
predMOJV$date <- as.POSIXct(as.Date(predMOJV$date, origin = "1970-01-01"))
predMOJV$upr <- predMOJV$fit + (1.96 * predMOJV$se.fit)
predMOJV$lwr <- predMOJV$fit - (1.96 * predMOJV$se.fit)

MOJVgamplot <- AZall %>% filter(!is.na(Sex), Age == "A") %>% ggplot(aes(x=`Sampling month`)) +
  geom_line(data=predMOJV, aes(x=date, y=fit, colour=Sex)) +
  geom_ribbon(data=predMOJV, aes(x=date, ymin=lwr, ymax=upr, fill=Sex), alpha=0.2) +
  ylim(-0.1,3) + theme_bw() + ylab("MojV log MFI values") +
  scale_color_manual(values=c("red","blue")) + theme_minimal() +
  geom_segment(aes(x = as.POSIXct("2020-07-01"), y = 2.8, 
                   xend = as.POSIXct("2020-07-01"), yend = 2),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = as.POSIXct("2019-07-01"), y = 2.8, 
                   xend = as.POSIXct("2019-07-01"), yend = 2),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm")))

ggarrange(HEVfem, HEVmal, HEVgamplot, HEVfem, HEVmal, HEVgamplot, 
          ncol=3, nrow=2, labels = c("A","B","C","D","E","F"))
#ggsave("FigS5.png", width = 300, height = 180, units= "mm", bg="white")


###tracking individual seroconversions

#look at zoo only, remove bats which were only captured once
AZ <- MFI.mock.corrected.wide %>% filter(Site=="Captive colony") %>%
  group_by(`Bat ID`) %>%
  dplyr::filter(!n() == 1) %>%
  ungroup()

#bats sampled at all 5 sampling points
fivetimes <- names(sort(table(AZ$`Bat ID`), decreasing=TRUE)[1:89])
AZ5 <- subset(AZ, `Bat ID` %in% fivetimes)
AZallage <- AZall %>% drop_na(`Age in Years`) %>% 
  dplyr::rename(age = `Age in Years`)
AZallage$age <- as.factor(AZallage$age)

AZ5$ID <- AZ5$`Bat ID`
AZ5$date <- as.numeric(as.Date(AZ5$`Sampling month`))

#GHV

AZ5_GHV <- AZ5 %>% dplyr::select(`Sampling month`, Sex, `Bat ID`, GHV) %>% 
  pivot_wider(names_from = `Sampling month`, values_from = GHV) %>%
  dplyr::select(Sex, `Bat ID`, `2019-02-01`, `2019-07-01`, `2020-02-01`,
                `2020-07-01`, `2020-10-01`)

AZ5_GHV$Jul19 <- NA
AZ5_GHV$Feb20 <- NA
AZ5_GHV$Jul20 <- NA
AZ5_GHV$Oct20 <- NA

for (i in 1:nrow(AZ5_GHV)) {
  AZ5_GHV$Jul19[i] <- as.numeric(AZ5_GHV[i,4]-AZ5_GHV[i,3])
  AZ5_GHV$Feb20[i] <- as.numeric(AZ5_GHV[i,5]-AZ5_GHV[i,4])
  AZ5_GHV$Jul20[i] <- as.numeric(AZ5_GHV[i,6]-AZ5_GHV[i,5])
  AZ5_GHV$Oct20[i] <- as.numeric(AZ5_GHV[i,7]-AZ5_GHV[i,6])
}

#seroconversion cutoff for nipah from Jax et al.: 0.97824008
##sensitivity analysis done at 0.5 and 1.5

for (i in 1:nrow(AZ5_GHV)) {
  for (j in 8:11) {
    if (AZ5_GHV[i,j] > 1.5) {
      AZ5_GHV[i,j] <- 1
    } else {
      AZ5_GHV[i,j] <- 0
    }
  }
}

table(AZ5_GHV$Jul19)
table(AZ5_GHV$Feb20)
table(AZ5_GHV$Jul20)
table(AZ5_GHV$Oct20)

print(AZ5_GHV[which(AZ5_GHV$Oct20 == 1),],n=23) #find individuals

#NIV 

AZ5_NIV <- AZ5 %>% dplyr::select(`Sampling month`, Sex, `Bat ID`, NIV) %>% 
  pivot_wider(names_from = `Sampling month`, values_from = NIV) %>%
  dplyr::select(Sex, `Bat ID`, `2019-02-01`, `2019-07-01`, `2020-02-01`,
                `2020-07-01`, `2020-10-01`)

AZ5_NIV$Jul19 <- NA
AZ5_NIV$Feb20 <- NA
AZ5_NIV$Jul20 <- NA
AZ5_NIV$Oct20 <- NA

for (i in 1:nrow(AZ5_NIV)) {
  AZ5_NIV$Jul19[i] <- as.numeric(AZ5_NIV[i,4]-AZ5_NIV[i,3])
  AZ5_NIV$Feb20[i] <- as.numeric(AZ5_NIV[i,5]-AZ5_NIV[i,4])
  AZ5_NIV$Jul20[i] <- as.numeric(AZ5_NIV[i,6]-AZ5_NIV[i,5])
  AZ5_NIV$Oct20[i] <- as.numeric(AZ5_NIV[i,7]-AZ5_NIV[i,6])
}


for (i in 1:nrow(AZ5_NIV)) {
  for (j in 8:11) {
    if (AZ5_NIV[i,j] > 1.5) {
      AZ5_NIV[i,j] <- 1
    } else {
      AZ5_NIV[i,j] <- 0
    }
  }
}

table(AZ5_NIV$Jul19)
table(AZ5_NIV$Feb20)
table(AZ5_NIV$Jul20)
table(AZ5_NIV$Oct20)

AZ5_NIV[which(AZ5_NIV$Oct20 == 1),] #find individuals


#HEV

AZ5_HEV <- AZ5 %>% dplyr::select(`Sampling month`, Sex, `Bat ID`, HEV) %>% 
  pivot_wider(names_from = `Sampling month`, values_from = HEV) %>%
  dplyr::select(Sex, `Bat ID`, `2019-02-01`, `2019-07-01`, `2020-02-01`,
                `2020-07-01`, `2020-10-01`)

AZ5_HEV$Jul19 <- NA
AZ5_HEV$Feb20 <- NA
AZ5_HEV$Jul20 <- NA
AZ5_HEV$Oct20 <- NA

for (i in 1:nrow(AZ5_HEV)) {
  AZ5_HEV$Jul19[i] <- as.numeric(AZ5_HEV[i,4]-AZ5_HEV[i,3])
  AZ5_HEV$Feb20[i] <- as.numeric(AZ5_HEV[i,5]-AZ5_HEV[i,4])
  AZ5_HEV$Jul20[i] <- as.numeric(AZ5_HEV[i,6]-AZ5_HEV[i,5])
  AZ5_HEV$Oct20[i] <- as.numeric(AZ5_HEV[i,7]-AZ5_HEV[i,6])
}

for (i in 1:nrow(AZ5_HEV)) {
  for (j in 8:11) {
    if (AZ5_HEV[i,j] > 1.5) {
      AZ5_HEV[i,j] <- 1
    } else {
      AZ5_HEV[i,j] <- 0
    }
  }
}

table(AZ5_HEV$Jul19)
table(AZ5_HEV$Feb20)
table(AZ5_HEV$Jul20)
table(AZ5_HEV$Oct20)

AZ5_HEV[which(AZ5_HEV$Oct20 == 1),] #find individuals


#MOJV seroconversions?

AZ5_MOJV <- AZ5 %>% dplyr::select(`Sampling month`, Sex, `Bat ID`, MOJV) %>% 
  pivot_wider(names_from = `Sampling month`, values_from = MOJV) %>%
  dplyr::select(Sex, `Bat ID`, `2019-02-01`, `2019-07-01`, `2020-02-01`,
                `2020-07-01`, `2020-10-01`)

AZ5_MOJV$Jul19 <- NA
AZ5_MOJV$Feb20 <- NA
AZ5_MOJV$Jul20 <- NA
AZ5_MOJV$Oct20 <- NA

for (i in 1:nrow(AZ5_MOJV)) {
  AZ5_MOJV$Jul19[i] <- as.numeric(AZ5_MOJV[i,4]-AZ5_MOJV[i,3])
  AZ5_MOJV$Feb20[i] <- as.numeric(AZ5_MOJV[i,5]-AZ5_MOJV[i,4])
  AZ5_MOJV$Jul20[i] <- as.numeric(AZ5_MOJV[i,6]-AZ5_MOJV[i,5])
  AZ5_MOJV$Oct20[i] <- as.numeric(AZ5_MOJV[i,7]-AZ5_MOJV[i,6])
}

for (i in 1:nrow(AZ5_MOJV)) {
  for (j in 8:11) {
    if (AZ5_MOJV[i,j] > 1.5) {
      AZ5_MOJV[i,j] <- 1
    } else {
      AZ5_MOJV[i,j] <- 0
    }
  }
}

table(AZ5_MOJV$Jul19)
table(AZ5_MOJV$Feb20)
table(AZ5_MOJV$Jul20)
table(AZ5_MOJV$Oct20)

AZ5_MOJV[which(AZ5_MOJV$Oct20 == 1),] #find individuals; mostly female

#####CAPTIVE COLONY: LONG TERM

#age bins and age grouped data frames
AZallage <- AZallage %>% mutate(Agebin = case_when(age == 0 ~ "0-1",
                                                   age == 1 ~ "0-1",
                                                   age == 2 ~ "2-3",
                                                   age == 3 ~ "2-3",
                                                   age == 4 ~ "4-5",
                                                   age == 5 ~ "4-5",
                                                   age == 6 ~ "6-7",
                                                   age == 7 ~ "6-7",
                                                   age == 8 ~ "8-9",
                                                   age == 9 ~ "8-9",
                                                   age == 10 ~ "10-12",
                                                   age == 11 ~ "10-12",
                                                   age == 12 ~ "10-12"))

AZallage$Agebin <- factor(AZallage$Agebin, levels = 
                            c("0-1","2-3","4-5","6-7","8-9","10-12"))

#age boxplots

hevage <- ggplot(AZallage, aes(x=Agebin, y=HEV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex)) +
  ylab("HeV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

nivage <- ggplot(AZallage, aes(x=Agebin, y=NIV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex)) +
  ylab("NiV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ghvage <- ggplot(AZallage, aes(x=Agebin, y=GHV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex)) +
  ylab("GhV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mojvage <- ggplot(AZallage, aes(x=Agebin, y=MOJV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex)) +
  ylab("MojV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fig4 <- ggarrange(ghvage, nivage, labels = c("A", "B"), align="v",ncol=1)
#ggsave("Fig4.png", width = 150, height = 180, units= "mm", bg="white")

figs6 <- ggarrange(hevage, mojvage, labels = c("A", "B"), align="v",ncol=1)
#ggsave("FigS6.png", width = 150, height = 180, units= "mm", bg="white")

####LMEMs to replace catalytic model

AZagemod <- AZallage %>% mutate(Agebin = case_when(age == 0 ~ "0-1",
                                                   age == 1 ~ "0-1",
                                                   age == 2 ~ "2-3",
                                                   age == 3 ~ "2-3",
                                                   age == 4 ~ "4-5",
                                                   age == 5 ~ "4-5",
                                                   age == 6 ~ "6-7",
                                                   age == 7 ~ "6-7",
                                                   age == 8 ~ "8-12",
                                                   age == 9 ~ "8-12",
                                                   age == 10 ~ "8-12",
                                                   age == 11 ~ "8-12",
                                                   age == 12 ~ "8-12"))

AZagemod$Agebin <- factor(AZagemod$Agebin, levels = 
                            c("0-1","2-3","4-5","6-7","8-12"))

summary(lmer(GHV ~ Agebin*Sex + (1|`Bat ID`), data=AZagemod))
summary(lmer(NIV ~ Agebin*Sex + (1|`Bat ID`), data=AZagemod))
summary(lmer(HEV ~ Agebin*Sex + (1|`Bat ID`), data=AZagemod))
summary(lmer(MOJV ~ Agebin*Sex + (1|`Bat ID`), data=AZagemod))
