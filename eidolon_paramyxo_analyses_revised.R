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
library(stringr)
library(ggcorrplot)
library(ggdendro)
library(dendextend)

## set up data

setwd("~/Documents/PhD/eidolon/paramyxo submission/round 2/round 3")
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

#wilcoxon paired tests

wilcox.test(MFI.data$`Analyte 35 GHV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")
wilcox.test(MFI.data$`Analyte 46 NIV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")
wilcox.test(MFI.data$`Analyte 43 HEV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")
wilcox.test(MFI.data$`Analyte 29 MOJV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")
wilcox.test(MFI.data$`Analyte 53 CEDV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")
wilcox.test(MFI.data$`Analyte 28 MENV`, MFI.data$`Analyte 26 MOCK`, paired=TRUE, alternative="greater")

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

### Multivariate analyses

###all MFI values

MFI.mock.corrected$Analyte <- factor(MFI.mock.corrected$Analyte, levels = c("NIV","HEV","CEDV","GHV","MOJV","MENV"))
levels(MFI.mock.corrected$Analyte) <- c("NiV","HeV","CedV","GhV","MojV","MenV")

all <- ggplot(MFI.mock.corrected) + 
  geom_histogram(aes(x=log.MFI.mc), alpha=0.4) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(vars(Analyte)) + 
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

MFI.mock.corrected.wide.cor <- MFI.mock.corrected.wide.caps %>%
  select(MenV, MojV, GhV, CedV, HeV, NiV)
ggcor <- ggcorrplot(cor(MFI.mock.corrected.wide.cor))

hc <- hclust(dist(c(1,2,4,8,16,29)), "ave")
labels(hc) <- c("MenV","MojV","GhV","CedV","HeV","NiV")
ddata <- dendro_data(hc, type = "rectangle")
dendro<-ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  geom_text(data = ddata$labels, 
            aes(x = x-0.1, y = y-3, label = label), size = 4, vjust = 0) +
  scale_y_reverse(expand = c(0.2, 0)) +
  ylab(NULL) + xlab(NULL) + theme_dendro()

ggarrange(all, ggcor, dendro, labels = c("A", "B", "C"), widths = c(1, 0.7, 0.4), ncol=3)
#ggsave("Fig1.png", width = 320, height = 100, units= "mm", bg="white")

#pca

paramyxo.names <- bead.key %>% filter(Family=="Paramyxoviridae") %>% pull(Abbreviation)
paramyxo.columns <- which(names(MFI.mock.corrected.wide) %in% paramyxo.names)

paramyxo.PCA <- MFI.mock.corrected.wide.caps %>% 
  dplyr::select(all_of(paramyxo.columns)) %>% prcomp(scale. = TRUE)
fviz_pca_biplot(paramyxo.PCA,habillage=MFI.mock.corrected.wide.caps$Site, 
                addEllipses=TRUE, ellipse.level=0.95, label="var")
#ggsave("FigS1.jpg", width=10,height=6)

## empirical seroprevalence

paramyxo_cutoffs3 <- data.frame(matrix(ncol=5, nrow=1, 
                                       dimnames = list(c(1:1),c("NiV","GhV","HeV","MojV","CedV"))))

paramyxo_cutoffs3[1,1:5] <- c(median(MFI.mock.corrected.wide$NIV), 
                              median(MFI.mock.corrected.wide$GHV),
                              median(MFI.mock.corrected.wide$HEV),
                              median(MFI.mock.corrected.wide$MOJV),
                              median(MFI.mock.corrected.wide$CEDV))

paramyxo_medians <- data.frame()
MFI.mock.corrected$`Sampling month` <- as.Date(MFI.mock.corrected$`Sampling month`)

for (i in 1:5) {
  y <- MFI.mock.corrected %>% drop_na(Sex) %>%
    filter(Analyte==colnames(paramyxo_cutoffs3[i])) %>%
    mutate(Positive = (log.MFI.mc>paramyxo_cutoffs3[1,i])) %>%
    group_by(Site, `Sampling month`, Analyte, Sex) %>% 
    dplyr::summarise(Prevalence = mean(as.integer(Positive)),
                     Positives = sum(as.integer(Positive)),
                     n=n())
  paramyxo_medians <- rbind(paramyxo_medians, y)
}

paramyxo_medians$`Sampling month` <- as.Date(paramyxo_medians$`Sampling month`)

conf <- binconf(paramyxo_medians$Positives, paramyxo_medians$n, alpha=0.05)[,2:3]
paramyxo_medians <- cbind(paramyxo_medians, conf)

NIVmedian <- paramyxo_medians %>% filter(Analyte == "NiV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex, ymin = Lower, ymax = Upper)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_pointrange(aes(size=n), position=position_dodge(width=25)) + scale_size_continuous(range = c(0.1, 1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("NiV seroprevalence")

GHVmedian <- paramyxo_medians %>% filter(Analyte == "GhV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex, ymin = Lower, ymax = Upper)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_pointrange(aes(size=n), position=position_dodge(width=25)) + scale_size_continuous(range = c(0.1, 1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("GhV seroprevalence")

fig2 <- ggarrange(GHVmedian, NIVmedian, labels = c("A", "B"), align="v",ncol=1)
#ggsave("Fig2.png", width = 150, height = 234, units= "mm", bg="white")

HEVmedian <- paramyxo_medians %>% filter(Analyte == "HeV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex, ymin = Lower, ymax = Upper)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_pointrange(aes(size=n), position=position_dodge(width=25)) + scale_size_continuous(range = c(0.1, 1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("HeV seroprevalence")

MOJVmedian <- paramyxo_medians %>% filter(Analyte == "MojV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex, ymin = Lower, ymax = Upper)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_pointrange(aes(size=n), position=position_dodge(width=25)) + scale_size_continuous(range = c(0.1, 1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("MojV seroprevalence")

CEDVmedian <- paramyxo_medians %>% filter(Analyte == "CedV") %>%
  ggplot(aes(x=`Sampling month`, y=Prevalence, col=Sex, ymin = Lower, ymax = Upper)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_pointrange(aes(size=n), position=position_dodge(width=25)) + scale_size_continuous(range = c(0.1, 1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = as.Date("2019-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  geom_vline(aes(xintercept = as.Date("2020-07-01")), 
             data = subset(paramyxo_medians, Site == "Captive colony"),
             colour="black", linetype = "longdash") +
  ylab("CedV seroprevalence")

suppfig2 <- ggarrange(HEVmedian, MOJVmedian, CEDVmedian, 
                      labels = c("A", "B", "C"), ncol=2, nrow=2)
#ggsave("FigS2.png", width = 300, height = 234, units= "mm", bg="white")

####Boxplots of MFI values for all roosts

MFI.mock.corrected.wide$`Sampling month` <- as.Date(MFI.mock.corrected.wide$`Sampling month`)
MFI.mock.corrected.wide$Sex <- as.factor(MFI.mock.corrected.wide$Sex)
MFI.mock.corrected.wide$SexMonth <- interaction(MFI.mock.corrected.wide$`Sampling month`, MFI.mock.corrected.wide$Sex) 

ghvbox <- MFI.mock.corrected.wide %>% drop_na(Sex) %>%
  ggplot(aes(x=`Sampling month`, y=GHV, group=SexMonth)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_boxplot(aes(fill=Sex), width=50) +
  ylab("GhV log MFI values") +
  xlab("Sampling month") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

nivbox <- MFI.mock.corrected.wide %>% drop_na(Sex) %>%
  ggplot(aes(x=`Sampling month`, y=NIV, group=SexMonth)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_boxplot(aes(fill=Sex), width=50) +
  ylab("NiV log MFI values") +
  xlab("Sampling month") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hevbox <- MFI.mock.corrected.wide %>% drop_na(Sex) %>%
  ggplot(aes(x=`Sampling month`, y=HEV, group=SexMonth)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_boxplot(aes(fill=Sex), width=50) +
  ylab("HeV log MFI values") +
  xlab("Sampling month") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mojvbox <- MFI.mock.corrected.wide %>% drop_na(Sex) %>%
  ggplot(aes(x=`Sampling month`, y=MOJV, group=SexMonth)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_boxplot(aes(fill=Sex), width=50) +
  ylab("MojV log MFI values") +
  xlab("Sampling month") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cedvbox <- MFI.mock.corrected.wide %>% drop_na(Sex) %>%
  ggplot(aes(x=`Sampling month`, y=CEDV, group=SexMonth)) +
  theme_bw() + facet_wrap(vars(Site)) +
  geom_boxplot(aes(fill=Sex), width=50) +
  ylab("CedV log MFI values") +
  xlab("Sampling month") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggarrange(ghvbox, nivbox, hevbox, mojvbox, cedvbox, 
          ncol=2, nrow=3, labels = c("A","B","C","D","E"))
#ggsave("FigS3.png", width = 300, height = 400, units= "mm", bg="white")

###combining reproduction + sex variables
MFI.mock.corrected.wide$ReproSex <- MFI.mock.corrected.wide$Repro
MFI.mock.corrected.wide$ReproSex[which(MFI.mock.corrected.wide$Sex == "M")] <- "Male"
MFI.mock.corrected.wide$ReproSex[which(is.na(MFI.mock.corrected.wide$ReproSex))] <- "Non-preg female"
MFI.mock.corrected.wide$ReproSex <- factor(MFI.mock.corrected.wide$ReproSex, 
                                           levels = c("Non-preg female","Pregnant","Lactating","Male"))

#####GAMMs of wild roosts

MFI.mock.corrected.wide$date <- as.numeric(as.Date(MFI.mock.corrected.wide$`Sampling month`))
MFI.mock.corrected.wide$Sex <- as.factor(MFI.mock.corrected.wide$Sex)

summary(gamm(GHV~s(date, bs="cr", k=5, by=Site)+Sex+Age+BodyCondition,
          data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony") %>% drop_na(Sex))$gam)
summary(gamm(NIV~s(date, bs="cr", k=5, by=Site)+Sex+Age+BodyCondition,
             data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony") %>% drop_na(Sex))$gam)
summary(gamm(HEV~s(date, bs="cr", k=5, by=Site)+Sex+Age+BodyCondition,
             data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony") %>% drop_na(Sex))$gam)
summary(gamm(MOJV~s(date, bs="cr", k=5, by=Site)+Sex+Age+BodyCondition,
             data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony") %>% drop_na(Sex))$gam)
summary(gamm(CEDV~s(date, bs="cr", k=5, by=Site)+Sex+Age+BodyCondition,
             data=MFI.mock.corrected.wide %>% filter(Site != "Captive colony") %>% drop_na(Sex))$gam)

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
#ggsave("Fig3.png", width = 300, height = 180, units= "mm", bg="white")

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

###CEDV

#by sex
CEDVfem <- ggplot(AZall %>% filter(Age == "A", Sex == "F"), aes(x=`Sampling month`, y=CEDV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="red") + geom_point(colour="red") + 
  ylim(-0.7,2.5) +
  theme(legend.position = "none") + ylab("CedV log MFI values") + theme_minimal()

CEDVmal <- ggplot(AZall %>% filter(Age == "A", Sex == "M"), aes(x=`Sampling month`, y=CEDV, group=`Bat ID`)) +
  geom_line(linewidth=0.1, colour="blue") + geom_point(colour="blue") + 
  ylim(-0.7,2.5) +
  theme(legend.position = "none") + ylab("CedV log MFI values") + theme_minimal()

##GAMMs

cedv2<-((gamm(CEDV~s(date, bs="cr", k=5, by=Sex)+Sex,
              data=AZall %>% filter(Age == "A"),random=list(ID=~1))))

summary(cedv2$gam)

predCEDV <- bind_cols(newdata, as.data.frame(predict(cedv2$gam, newdata=newdata, type="link", se.fit=TRUE)))
predCEDV$date <- as.POSIXct(as.Date(predCEDV$date, origin = "1970-01-01"))
predCEDV$upr <- predCEDV$fit + (1.96 * predCEDV$se.fit)
predCEDV$lwr <- predCEDV$fit - (1.96 * predCEDV$se.fit)

CEDVgamplot <- AZall %>% filter(!is.na(Sex), Age == "A") %>% ggplot(aes(x=`Sampling month`)) +
  geom_line(data=predCEDV, aes(x=date, y=fit, colour=Sex)) +
  geom_ribbon(data=predCEDV, aes(x=date, ymin=lwr, ymax=upr, fill=Sex), alpha=0.2) +
  ylim(-0.7,2.5) + theme_bw() + ylab("CedV log MFI values") +
  scale_color_manual(values=c("red","blue")) + theme_minimal() +
  geom_segment(aes(x = as.POSIXct("2020-07-01"), y = 1.5, 
                   xend = as.POSIXct("2020-07-01"), yend = 0.7),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = as.POSIXct("2019-07-01"), y = 1.5, 
                   xend = as.POSIXct("2019-07-01"), yend = 0.7),
               lineend = "round", linejoin = "round",
               size=0.8, colour = "black",
               arrow = arrow(length = unit(0.5, "cm")))

ggarrange(HEVfem, HEVmal, HEVgamplot, MOJVfem, MOJVmal, MOJVgamplot,
          CEDVfem, CEDVmal, CEDVgamplot,
          ncol=3, nrow=3, labels = c("A","B","C","D","E","F","G","H","I"))
#ggsave("FigS4.png", width = 300, height = 270, units= "mm", bg="white")

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

#MOJV seroconversions

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

####CEDV seroconversion

AZ5_CEDV <- AZ5 %>% dplyr::select(`Sampling month`, Sex, `Bat ID`, CEDV) %>% 
  pivot_wider(names_from = `Sampling month`, values_from = CEDV) %>%
  dplyr::select(Sex, `Bat ID`, `2019-02-01`, `2019-07-01`, `2020-02-01`,
                `2020-07-01`, `2020-10-01`)

AZ5_CEDV$Jul19 <- NA
AZ5_CEDV$Feb20 <- NA
AZ5_CEDV$Jul20 <- NA
AZ5_CEDV$Oct20 <- NA

for (i in 1:nrow(AZ5_CEDV)) {
  AZ5_CEDV$Jul19[i] <- as.numeric(AZ5_CEDV[i,4]-AZ5_CEDV[i,3])
  AZ5_CEDV$Feb20[i] <- as.numeric(AZ5_CEDV[i,5]-AZ5_CEDV[i,4])
  AZ5_CEDV$Jul20[i] <- as.numeric(AZ5_CEDV[i,6]-AZ5_CEDV[i,5])
  AZ5_CEDV$Oct20[i] <- as.numeric(AZ5_CEDV[i,7]-AZ5_CEDV[i,6])
}

for (i in 1:nrow(AZ5_CEDV)) {
  for (j in 8:11) {
    if (AZ5_CEDV[i,j] > 0.5) {
      AZ5_CEDV[i,j] <- 1
    } else {
      AZ5_CEDV[i,j] <- 0
    }
  }
}

table(AZ5_CEDV$Jul19)
table(AZ5_CEDV$Feb20)
table(AZ5_CEDV$Jul20)
table(AZ5_CEDV$Oct20)

AZ5_CEDV[which(AZ5_CEDV$Jul20 == 1),]

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
  geom_boxplot(varwidth = TRUE, aes(fill=Sex), outlier.shape = NA) +
  geom_point(aes(color=Sex), position=position_jitterdodge(dodge.width = 0.7), size = 0.5) +
  ylab("HeV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

nivage <- ggplot(AZallage, aes(x=Agebin, y=NIV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex), outlier.shape = NA) +
  geom_point(aes(color=Sex), position=position_jitterdodge(dodge.width = 0.7), size = 0.5) +
  ylab("NiV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ghvage <- ggplot(AZallage, aes(x=Agebin, y=GHV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex), outlier.shape = NA) +
  geom_point(aes(color=Sex), position=position_jitterdodge(dodge.width = 0.7), size = 0.5) +
  ylab("GhV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mojvage <- ggplot(AZallage, aes(x=Agebin, y=MOJV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex), outlier.shape = NA) +
  geom_point(aes(color=Sex), position=position_jitterdodge(dodge.width = 0.7), size = 0.5) +
  ylab("MojV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cedvage <- ggplot(AZallage, aes(x=Agebin, y=CEDV)) + 
  theme_minimal() +
  geom_boxplot(varwidth = TRUE, aes(fill=Sex), outlier.shape = NA) +
  geom_point(aes(color=Sex), position=position_jitterdodge(dodge.width = 0.7), size = 0.5) +
  ylab("CedV log MFI values") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fig4 <- ggarrange(ghvage, nivage, labels = c("A", "B"), align="v",ncol=1)
#ggsave("Fig4.png", width = 150, height = 180, units= "mm", bg="white")

figs5 <- ggarrange(hevage, mojvage, cedvage, labels = c("A", "B", "C"), align="v",ncol=2, nrow=2)
#ggsave("FigS5.png", width = 200, height = 180, units= "mm", bg="white")

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

AZagemod$Repro <- as.factor(AZagemod$Repro)
summary(lmer(GHV ~ Agebin + ReproSex + BodyCondition + (1|`Bat ID`), data=AZagemod))
summary(lmer(NIV ~ Agebin + ReproSex + BodyCondition + (1|`Bat ID`), data=AZagemod))
summary(lmer(HEV ~ Agebin + ReproSex + BodyCondition + (1|`Bat ID`), data=AZagemod))
summary(lmer(MOJV ~ Agebin + ReproSex + BodyCondition + (1|`Bat ID`), data=AZagemod))
summary(lmer(CEDV ~ Agebin + ReproSex + BodyCondition + (1|`Bat ID`), data=AZagemod))
