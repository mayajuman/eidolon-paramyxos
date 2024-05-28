## eidolon helvum: paramyxovirus data processing and supp data generation
## maya juman
## updated 02/09/24

## clean up space
rm(list=ls()) 
graphics.off()

## set up

library(readxl)
library(writexl)
library(plyr)
library(dplyr)

setwd("~/Documents/PhD/eidolon")
MFI.data <- read_excel("Eidolon serology PREEMPT 2021-12-07.xlsx",sheet="samples", guess_max = 5000) %>% 
  mutate(Site = factor(Site)) %>%
  dplyr::rename(`Sampling month` = `Date (MM-YY)`) %>%
  filter(Sample != "ZJU20-143 SE")
#removed this sample because the bat was sampled twice in a month, and this sample had weird measurements

####associating data from sample sheets and data cleaning

filenames <- list.files("sample sheets", pattern="*.xlsx", full.names=TRUE)
ldf <- lapply(filenames, read_excel)
ldf[[5]] <- ldf[[5]][-1,]

for (i in 1:18) {
  names(ldf[[i]]) <- unlist(ldf[[i]][6,])
  ldf[[i]] <- ldf[[i]][-c(1:6),]
}

ldf[[19]] <- ldf[[19]] %>%
  dplyr::rename('Sample Code' = sample_ID, "Bat ID" = "Bat ID (Last 7/8 digits)")

#deal with weird spacing in sample ID and extract variables of interest
#clean up Bat ID column

ldf[[11]]$`Bat ID (Last 7/8 digits)` <- NA

for (i in 1:12) {
  names(ldf[[i]])[names(ldf[[i]]) == "Bat ID (Last 7/8 digits)"] <- "Bat ID"
}

for (i in 13:18) {
  names(ldf[[i]])[names(ldf[[i]]) == "Bat ID (Last 7 digits)"] <- "Bat ID"
}

for (i in 1:19) {
  ldf[[i]]$`Sample Code` <- gsub(" ", "", ldf[[i]]$`Sample Code`, fixed = TRUE)
  ldf[[i]] <- dplyr::select(ldf[[i]], 
                            c("Sample Code","Species","Age","Sex",
                              "Repro. Status","Bat weight/g","Forearm/mm",
                              "Bat ID"))
  ldf[[i]]$`Sample Code` <- as.character(ldf[[i]]$`Sample Code`)
  ldf[[i]]$`Species` <- as.factor(ldf[[i]]$`Species`)
  ldf[[i]]$`Age` <- as.character(ldf[[i]]$`Age`)
  ldf[[i]]$`Sex` <- as.factor(ldf[[i]]$`Sex`)
  ldf[[i]]$`Bat ID` <- as.character(ldf[[i]]$`Bat ID`)
  ldf[[i]]$`Repro. Status` <- as.character(ldf[[i]]$`Repro. Status`)
  ldf[[i]]$`Bat weight/g` <- as.numeric(ldf[[i]]$`Bat weight/g`)
  ldf[[i]]$`Forearm/mm` <- as.numeric(ldf[[i]]$`Forearm/mm`)
  ldf[[i]] <- ldf[[i]] %>% dplyr::rename('Sample' = 'Sample Code')
  ldf[[i]] <- ldf[[i]] %>% dplyr::rename('Repro' = 'Repro. Status')
  ldf[[i]] <- ldf[[i]] %>% dplyr::rename('Weight' = 'Bat weight/g')
  ldf[[i]] <- ldf[[i]] %>% dplyr::rename('Forearm' = 'Forearm/mm')
}

demdata <- bind_rows(ldf)

#cleaning up
demdata <- demdata[-which(is.na(demdata$'Sample')),]

demdata$Repro[demdata$Repro == "lactating"] <- "Lactating"
demdata$Repro[demdata$Repro == "pregnant" | demdata$Repro == "Preg"] <- "Pregnant"
demdata$Repro[demdata$Repro == "-"] <- NA
demdata$Weight[demdata$Sample == "ZOC20-076"] <- NA

demdata$Repro[is.na(demdata$Repro)] <- "-"
demdata$Repro <- as.factor(demdata$Repro)

demdata$Age[demdata$Age == "JUV"] <- "J"
demdata$Age[demdata$Age == "Juv"] <- "J"
demdata$Age <- as.factor(demdata$Age)
demdata$`Bat ID`[demdata$`Bat ID` == "-"] <- NA

demdata$Weight[demdata$Weight<0] <- NA
demdata$Forearm[demdata$Forearm > 1000] <- 114

###fixing tag typos

demdata$`Bat ID`[demdata$`Bat ID` == 12062095] <- 12062096
demdata$`Bat ID`[demdata$`Bat ID` == 12067582] <- 12067580
demdata$`Bat ID`[demdata$`Bat ID` == 13063052] <- 12063052
demdata$`Bat ID`[demdata$`Bat ID` == 1607829] <- 1607828
demdata$`Bat ID`[demdata$`Bat ID` == 4895894] <- 5895894
demdata$`Bat ID`[demdata$`Bat ID` == 5896393] <- 5896939
demdata$`Bat ID`[demdata$`Bat ID` == 5896940] <- 5896540
demdata$`Bat ID`[demdata$`Bat ID` == 6898478] <- 5898478
demdata$`Bat ID`[demdata$`Bat ID` == 2902876] <- 5902876

###fixing sexes

demdata$Sex[demdata$`Bat ID` == 12061432] <- "M"
demdata$Sex[demdata$`Bat ID` == 12063180] <- "F"
demdata$Sex[demdata$`Bat ID` == 12063497] <- "M"
demdata$Sex[demdata$`Bat ID` == 12064970] <- "F"
demdata$Sex[demdata$`Bat ID` == 12070425] <- "F"
demdata$Sex[demdata$`Bat ID` == 1933940] <- "F"
demdata$Sex[demdata$`Bat ID` == 5897532] <- "M"
demdata$Sex[demdata$`Bat ID` == 5897884] <- NA
demdata$Sex[demdata$`Bat ID` == 5907743] <- "F"

###fixing ages

demdata$Age[demdata$`Bat ID` == 5903706] <- "A"
demdata$Age[demdata$`Bat ID` == 1591879] <- "A"
demdata$Age[demdata$`Bat ID` == 1607615] <- "A"
demdata$Age[demdata$`Bat ID` == 1925151] <- "A"
demdata$Age[demdata$`Bat ID` == 3461662] <- "A"
demdata$Age[demdata$`Bat ID` == 3470474] <- "A"
demdata$Age[demdata$`Bat ID` == 3476475] <- "A"
demdata$Age[demdata$`Bat ID` == 3594197] <- "A"
demdata$Age[demdata$Sample == "ZFE19-098"] <- "SA"
demdata$Age[demdata$`Bat ID` == 5897203] <- "A"
demdata$Age[demdata$`Bat ID` == 5897532] <- "A"
demdata$Age[demdata$`Bat ID` == 5898663] <- "A"
demdata$Age[demdata$`Bat ID` == 5899054] <- "A"
demdata$Age[demdata$`Bat ID` == 5903706] <- "A"

MFI.data$Sample <- gsub(" SE", "", MFI.data$Sample, fixed = TRUE)
MFI.data <- MFI.data %>% left_join(demdata, by='Sample')

####remove non eidolon
#also remove sample ZOC20-019, which was a repeated sample with the wrong forearm length for same individual
MFI.data <- MFI.data %>% filter(Species == "E.helvum", Sample !="ZOC20-019")
rm(demdata, ldf, filenames, i)

levels(MFI.data$Site) <- c("Accra urban roost","Captive colony",
                           "Akosombo rural roost","Kumasi urban roost")

bead.key <- read_excel("Eidolon serology PREEMPT 2021-12-07.xlsx",sheet="Bead Key")
bead.key <- bead.key %>% filter(Family != "Filoviridae")
bead.key$Genus[which(bead.key$Abbreviation == "MENV")] <- "Pararubulavirus"

controls <- read_excel("Eidolon serology PREEMPT 2021-12-07.xlsx",sheet="controls",na = "NaN")
controls <- controls[-c(1:5, 146:290),]
controls <- controls %>% dplyr::select(Sample, `Analyte 29 MOJV`, `Analyte 35 GHV`,
                                       `Analyte 43 HEV`, `Analyte 46 NIV`,
                                       `Analyte 53 CEDV`, `Analyte 28 MENV`,
                                       `Analyte 26 MOCK`, `Luminex Date`)
controls$`Luminex Date` <- as.character(controls$`Luminex Date`)

ages <- read.csv("ages.csv")
ages <- ages %>% dplyr::rename("Bat ID" = Bat.ID, "Birthyear" = Year.of.birth) %>%
  dplyr::select(`Bat ID`, Birthyear)
ages$`Bat ID` <- as.character(ages$`Bat ID`)
MFI.data <- merge(MFI.data, ages, by="Bat ID", all=T)
MFI.data$`Sampling month` <- as.character(MFI.data$`Sampling month`)

MFI.data$`Age in Years` <- NA
for (i in 1:nrow(MFI.data)) {
  if(!is.na(MFI.data$Birthyear[i])){
  MFI.data$`Age in Years`[i] <- as.numeric(strsplit(MFI.data$`Sampling month`[i], split="-")[[1]][1]) - MFI.data$Birthyear[i]
  }
}

MFI.data$Repro[which(MFI.data$Repro == "-")] <- NA

MFI.data <- MFI.data %>% dplyr::select(`Bat ID`, Site, `Sampling month`, Sample,
                                       `Analyte 29 MOJV`, `Analyte 35 GHV`,
                                       `Analyte 43 HEV`, `Analyte 46 NIV`,
                                       `Analyte 53 CEDV`, `Analyte 28 MENV`,
                                       `Analyte 26 MOCK`, Age, `Age in Years`,
                                       Sex, Repro, Weight, Forearm) %>%
  dplyr::rename(`Reproductive Status` = `Repro`,
                `Weight (g)` = Weight,
                `Forearm (mm)` = Forearm,
                `Age Category` = Age)

suppdata <- list("Key" = bead.key, "Data" = MFI.data, "Controls" = controls) #assume sheet1 and sheet2 are data frames
write_xlsx(suppdata, "Supplementary Data.xlsx")
