## Robustness Check DID ##

# Setup ----------------------------------------------------------

library(tidyverse)
library(sf)
library(spdep)
library(AggregateR)
library(rlist)
library(gtsummary)
library(readxl)
library(data.table)
library(moments)
#for rasterization
library(stars)
library(terra)
library(raster)
library(rmapshaper)
#for connectivity metrics
library(sfnetworks)
library(tidygraph)
library(igraph)
#for estimation
library(AER)
library(moments)
library(lmtest)
library(sandwich)
#for regression output
library(texreg)
#no scientific notation
options(scipen=999) 

# Load functions ----------------------------------------------------------

source("./Functions/Graph_function.R")
source("./Functions/NeighborCon_function.R")

# Load data ---------------------------------------------------------------

#Districts
Districts = read.csv("./Data/GemeindenKtBE.csv")

#Biodiversity promotion sites AES
AES = st_read("./Data/BFF1_BFF2_Vernetzung_Fläche.shp")
# CRS: CH1903+ / LV95

# Sites -------------------------------------------------------------------

#Only inlcude sites in specific areas of intervention
AES = AES[AES$Massnahmen %in% c("Vernetzungsgebiet Hügel / Hang","Erhaltungsgebiet strukturreiche Landschaft","Vernetzungsgebiet Tal / offenes Agrarland", "Vernetzungsgebiet offene Wiesenlandschaft"),]

#Only inlcude sites 611,612,617 smaller than 0.3 ha
AES = AES[AES$CODE %in% c("611","612","617") & AES$VERNETZUNG < 30,]

#Divide into control and treatment
With_rule = c('Vernetzungsgebiet Tal / offenes Agrarland','Vernetzungsgebiet Hügel / Hang','Vernetzungsgebiet offene Wiesenlandschaft')
AES$Rule = as.factor(ifelse(AES$Massnahmen %in% With_rule, 1, 0))

#Test
table(AES$Rule) [[2]] == dim(filter(AES, Massnahmen %in% With_rule)) [1]

#Adding district identification to the dataset
AES = left_join(AES,Districts, by = c('Gemeinde_I' = 'BFS'), keep = FALSE)

#Exclude sites outside the canton of Bern
AES = AES[!AES$Gemeinde_I >= 1000 & !AES$Gemeinde_I <301,]

#Exclude sites where treatment and control overlap
AES$AES_ID = c(1:dim(AES)[1])
AES0 = AES[AES$Rule==0,]
AES1 = AES[AES$Rule==1,]
inter_AES = st_intersection(AES0, AES1, model = "closed")
inter_AES_ID = sort(rbind(inter_AES$AES_ID, inter_AES$AES_ID.1))
AES2 = AES[!AES$AES_ID %in% inter_AES_ID,]

#Add an ID combining municipality and treatment status 
AES$Municipality_rule = paste(AES$Gemeinde_I, AES$Rule)

#Exclude unused rows, rename
AES_A = AES %>% 
  select(Gemeinde_I,Betriebs_I,VK_Nr,CODE,VERNETZUNG,Kulturname,Massnahmen,Rule, Municipality_rule) %>% 
  rename(Municipality_ID = Gemeinde_I, Farm_ID = Betriebs_I, District_ID = VK_Nr, Site_types_ID = CODE, Site_types_name = Kulturname, Area_of_intervention = Massnahmen)

#Municipalities as in the main analysis
AES_A <- AES_A[AES_A$District_ID %in% c(242:246), ]

#Dataset with a network bonus payment NBP
NBP = filter(AES_A, AES_A$VERNETZUNG > 0)

#Dataset with biodiversity promotion areas with no NBP
BPA = filter(AES_A, AES_A$VERNETZUNG == 0)


# Connectivity Metrics BPA ----------------------------------------------------

# Compute Gamma Control
BPA0 <- BPA[BPA$Rule == 0, ]
Municipality_ID_sorted <- sort(unique(BPA0$Municipality_ID))
Gamma <- numeric(length(Municipality_ID_sorted))

for (i in 1:length(Municipality_ID_sorted)) {
  BPA_Com <- BPA0[BPA0$Municipality_ID == Municipality_ID_sorted[i],
  ]
  Gamma[i] <- Connectance(BPA_Com)
}

Gamma0 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 0, Gamma))
Gamma0$ID <- paste(Gamma0$Municipality_ID_sorted, Gamma0$Rule)

# Compute Gamma Treatment
BPA1 <- BPA[BPA$Rule == 1, ]
Municipality_ID_sorted <- sort(unique(BPA1$Municipality_ID))
Gamma <- numeric(length(Municipality_ID_sorted))

for (i in 1:length(Municipality_ID_sorted)) {
  BPA_Com <- BPA1[BPA1$Municipality_ID == Municipality_ID_sorted[i],
  ]
  Gamma[i] <- Connectance(BPA_Com)
}

Gamma1 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 1, Gamma))
Gamma1$ID <- paste(Gamma1$Municipality_ID_sorted, Gamma1$Rule)

Gamma <- rbind(Gamma0, Gamma1)

Missing <- Gamma[!complete.cases(Gamma), ] #check NA values
Gamma[is.na(Gamma)] <- 0  #the NaN values arise due to 0/0 cases
Gamma <- rename(Gamma, Connectance_BPA = Gamma)

# Compute NC Control
NCon <- by(BPA0, BPA0$Municipality_ID, NeighborCon)

Municipality_ID_sorted <- sort(unique(BPA0$Municipality_ID))
NC <- numeric(length(Municipality_ID_sorted))
for (i in 1:length(Municipality_ID_sorted)) {
  if (length(NCon[[i]]) == 1) {
    NC[i] <- 0
  } else {
    NC[i] <- mean(NCon[[i]][, 2])/(mean(NCon[[i]][, 1] + NCon[[i]][,
                                                                   2]))
  }
}

NC0 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 0, NC))
NC0$ID <- paste(NC0$Municipality_ID_sorted, NC0$Rule)

# Compute NC Treatment
NCon <- by(BPA1, BPA1$Municipality_ID, NeighborCon)

Municipality_ID_sorted <- sort(unique(BPA1$Municipality_ID))
NC <- numeric(length(Municipality_ID_sorted))
for (i in 1:length(Municipality_ID_sorted)) {
  if (length(NCon[[i]]) == 1) {
    NC[i] <- 0
  } else {
    NC[i] <- mean(NCon[[i]][, 2])/(mean(NCon[[i]][, 1] + NCon[[i]][,
                                                                   2]))  #mean of external over overall mean
  }
}

NC1 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 1, NC))
NC1$ID <- paste(NC1$Municipality_ID_sorted, NC1$Rule)

NC <- rbind(NC0, NC1) 
Missing <- NC[!complete.cases(NC), ] #check NA values
NC[is.na(NC)] <- 0 #the NA values arise due to 0/0 cases
NC <- rename(NC, NC_BPA = NC)

#load MG_A
load("./Data/MG_A.Rda")

#Merge with MG
MG = left_join(MG_A, Gamma[,c("Connectance_BPA","ID")], by=c("Municipality_rule"="ID"))
MG <- rename(MG, Connectance_NBP = Connectance)
MG = left_join(MG, NC[,c("NC_BPA","ID")], by=c("Municipality_rule"="ID"))
MG <- rename(MG, NC_NBP = NC)

#Check for areas where there are no BPA sites
Missing <- MG[!complete.cases(MG), ]

#Exclude these observations
MG <- MG[complete.cases(MG), ]

#Reshape dataset so that we have connectance and a dummy whether it is BPA or NBP
MG_long <- MG %>%
  pivot_longer(
    cols = c("Connectance_BPA", "Connectance_NBP", "NC_BPA", "NC_NBP"),
    names_to = c(".value", "source"),
    names_pattern = "(.*)_([^_]*)$"
  )
View(MG_long)

MG_diffdiff <- mutate(MG_long,
                      Policy = ifelse(source == "NBP", 1, 0),
)

# Difference-in-Differences Estimation for Connectance---------------------------

hist(MG_diffdiff$Connectance)

MG_diffdiff$ConnectanceL = log(MG_diffdiff$Connectance)
summary_lny = summary(MG_diffdiff$ConnectanceL[is.finite(MG_diffdiff$ConnectanceL)])
min_value = summary_lny[1]
gamma = min_value - 0.000001
MG_diffdiff$ConnectanceL = ifelse(MG_diffdiff$ConnectanceL == -Inf, gamma, MG_diffdiff$ConnectanceL)

MG_diffdiff$ConnectanceL <- MG_diffdiff$ConnectanceL*100

Connectance_min <- min(MG_diffdiff$ConnectanceL)  #LowerLimit

hist(MG_diffdiff$ConnectanceL)

#Invert farms
MG_diffdiff$Farms_nR <- 1/MG_diffdiff$Farms_n

#Tobit
did_model_connectance_1 <- tobit(ConnectanceL ~ Rule + Policy + Rule:Policy + Farm_size + Farms_nR + Forest_pct + Pasture + SHDI +
                                   Ascent_AIV + Soil_sui + Tree_per_ha, left = Connectance_min, data = MG_diffdiff, robust = TRUE)
summary(did_model_connectance_1)

ks.test(resid(did_model_connectance_1), "pnorm", mean=mean(resid(did_model_connectance_1)), sd=sd(resid(did_model_connectance_1)))
skewness(resid(did_model_connectance_1))
kurtosis(resid(did_model_connectance_1))
hist(resid(did_model_connectance_1))
vif(did_model_connectance_1)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_connectance_1)

# Check confidence intervals with bootstrapping: 
set.seed(888)
fit_b <- Boot(did_model_connectance_1, R = 5000)
summary(fit_b)
confint(fit_b, level = .95)

png(filename="did_model_connectance_1_bootstrpped.png", width = 900, height = 1300)
par(mar=c(5.1, 4.1, 6.1, 2.1), oma=c(0, 0, 2, 0), cex=1.8)
hist(fit_b, 12)
dev.off()

#Linear
did_model_connectance_2 <- lm(ConnectanceL ~ Rule + Policy + Rule:Policy + Farm_size + Farms_nR + Forest_pct + Pasture + SHDI +
                                Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_diffdiff)
summary(did_model_connectance_2)
BIC(did_model_connectance_2)

did_model_connectance_2_robust <- coeftest(did_model_connectance_2, vcov = vcovHC(did_model_connectance_2, type = "HC0"))

ks.test(resid(did_model_connectance_2), "pnorm", mean=mean(resid(did_model_connectance_2)), sd=sd(resid(did_model_connectance_2)))

skewness(resid(did_model_connectance_2))
kurtosis(resid(did_model_connectance_2))
hist(resid(did_model_connectance_2))
vif(did_model_connectance_2)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_connectance_2)

# Difference-in-Differences Estimation for Neighbor Connectivity-----------------

MG_diffdiff$NC <- MG_diffdiff$NC*100

hist(MG_diffdiff$NC)


#Tobit
did_model_nc_1 <- tobit(NC ~ Rule + Policy + Rule:Policy + Farm_size + Farms_n + Forest_pct + Pasture + SHDI +
                          Ascent_AIV + Soil_sui + Tree_per_ha, left = Connectance_min, data = MG_diffdiff, robust = TRUE)
summary(did_model_nc_1)
ks.test(resid(did_model_nc_1), "pnorm", mean=mean(resid(did_model_nc_1)), sd=sd(resid(did_model_nc_1)))

skewness(resid(did_model_nc_1))
kurtosis(resid(did_model_nc_1))
hist(resid(did_model_nc_1))
vif(did_model_nc_1)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_nc_1)

set.seed(888)
fit_b <- Boot(did_model_nc_1, R = 5000)
summary(fit_b)
confint(fit_b, level = .95)
png(filename="did_model_nc_1_bootstrpped.png", width = 900, height = 1300)
par(mar=c(5.1, 4.1, 6.1, 2.1), oma=c(0, 0, 2, 0), cex=1.8)
hist(fit_b, 12)
dev.off()


#Linear
did_model_nc_2 <- lm(NC ~ Rule + Policy + Rule:Policy + Farm_size + Farms_n + Forest_pct + Pasture + SHDI +
                       Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_diffdiff)
summary(did_model_nc_2)
BIC(did_model_nc_2)

did_model_nc_2_robust <- coeftest(did_model_nc_2, vcov = vcovHC(did_model_nc_2, type = "HC0"))

ks.test(resid(did_model_nc_2), "pnorm", mean=mean(resid(did_model_nc_2)), sd=sd(resid(did_model_nc_2)))

skewness(resid(did_model_nc_2))
kurtosis(resid(did_model_nc_2))
hist(resid(did_model_nc_2))
vif(did_model_nc_2)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_nc_2)

# DID with the matched sample--------------------------------------------

#load MG_A
load("./Data/MG_matched.Rda")

#Merge with MG
MG = left_join(MG_matched, Gamma[,c("Connectance_BPA","ID")], by=c("Municipality_rule"="ID"))
MG <- rename(MG, Connectance_NBP = Connectance)
MG = left_join(MG, NC[,c("NC_BPA","ID")], by=c("Municipality_rule"="ID"))
MG <- rename(MG, NC_NBP = NC)

#Check for areas where there are no BPA sites
Missing <- MG[!complete.cases(MG), ]

#Exclude these observations
MG <- MG[complete.cases(MG), ]

#Reshape dataset so that we have connectance and a dummy whether it is BPA or NBP
MG_long <- MG %>%
  pivot_longer(
    cols = c("Connectance_BPA", "Connectance_NBP", "NC_BPA", "NC_NBP"),
    names_to = c(".value", "source"),
    names_pattern = "(.*)_([^_]*)$"
  )


MG_diffdiff <- mutate(MG_long,
                      Policy = ifelse(source == "NBP", 1, 0),
)

# Difference-in-Differences Estimation for Connectance---------------------------

hist(MG_diffdiff$Connectance)

MG_diffdiff$ConnectanceL = log(MG_diffdiff$Connectance)
summary_lny = summary(MG_diffdiff$ConnectanceL[is.finite(MG_diffdiff$ConnectanceL)])
min_value = summary_lny[1]
gamma = min_value - 0.000001
MG_diffdiff$ConnectanceL = ifelse(MG_diffdiff$ConnectanceL == -Inf, gamma, MG_diffdiff$ConnectanceL)

MG_diffdiff$ConnectanceL <- MG_diffdiff$ConnectanceL*100

Connectance_min <- min(MG_diffdiff$ConnectanceL)  #LowerLimit

hist(MG_diffdiff$ConnectanceL)
#MG_diffdiff = na.omit(MG_diffdiff)

#Invert farms
MG_diffdiff$Farms_nR <- 1/MG_diffdiff$Farms_n

#Tobit
did_model_connectance_1_matched <- tobit(ConnectanceL ~ Rule + Policy + Rule:Policy + Farm_size + Farms_nR + Forest_pct + Pasture + SHDI +
                                           Ascent_AIV + Soil_sui + Tree_per_ha, left = Connectance_min, data = MG_diffdiff, robust = TRUE, weights = weights)
summary(did_model_connectance_1_matched)
ks.test(resid(did_model_connectance_1_matched), "pnorm", mean=mean(resid(did_model_connectance_1_matched)), sd=sd(resid(did_model_connectance_1_matched)))

skewness(resid(did_model_connectance_1_matched))
kurtosis(resid(did_model_connectance_1_matched))
hist(resid(did_model_connectance_1_matched))
vif(did_model_connectance_1_matched)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_connectance_1_matched)

set.seed(888)
fit_b <- Boot(did_model_connectance_1_matched, R = 5000)
summary(fit_b)
confint(fit_b, level = .95)
png(filename="did_model_connectance_1_matched_bootstrpped.png", width = 900, height = 1300)
par(mar=c(5.1, 4.1, 6.1, 2.1), oma=c(0, 0, 2, 0),
    cex=1.8)
hist(fit_b, 12)
dev.off()

#Linear
did_model_connectance_2_matched <- lm(ConnectanceL ~ Rule + Policy + Rule:Policy + Farm_size + Farms_nR + Forest_pct + Pasture + SHDI +
                                        Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_diffdiff, weights = weights)
summary(did_model_connectance_2_matched)

did_model_connectance_2_matched_robust <- coeftest(did_model_connectance_2_matched, vcov = vcovHC(did_model_connectance_2_matched, type = "HC0"))

summary(did_model_connectance_2_matched_robust)
BIC(did_model_connectance_2_matched_robust)

ks.test(resid(did_model_connectance_2_matched), "pnorm", mean=mean(resid(did_model_connectance_2_matched)), sd=sd(resid(did_model_connectance_2_matched)))

skewness(resid(did_model_connectance_2_matched))
kurtosis(resid(did_model_connectance_2_matched))
hist(resid(did_model_connectance_2_matched))
vif(did_model_connectance_2_matched)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_connectance_2_matched)

# Difference-in-Differences Estimation for Neighbor Connectivity-----------------

MG_diffdiff$NC <- MG_diffdiff$NC*100

hist(MG_diffdiff$NC)


did_model_nc_1_matched <- tobit(NC ~ Rule + Policy + Rule:Policy + Farm_size + Farms_n + Forest_pct + Pasture + SHDI +
                                  Ascent_AIV + Soil_sui + Tree_per_ha, left = Connectance_min, data = MG_diffdiff, robust = TRUE, weights = weights)
summary(did_model_nc_1_matched)
ks.test(resid(did_model_nc_1_matched), "pnorm", mean=mean(resid(did_model_nc_1_matched)), sd=sd(resid(did_model_nc_1_matched)))

skewness(resid(did_model_nc_1_matched))
kurtosis(resid(did_model_nc_1_matched))
hist(resid(did_model_nc_1_matched))
vif(did_model_nc_1_matched)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_nc_1_matched)

set.seed(888)
fit_b <- Boot(did_model_nc_1_matched, R = 5000)
summary(fit_b)
confint(fit_b, level = .95)
png(filename="did_model_nc_1_matched_bootstrpped.png", width = 900, height = 1300)
par(mar=c(5.1, 4.1, 6.1, 2.1), oma=c(0, 0, 2, 0),cex=1.8)
hist(fit_b, 12)
dev.off()

#Linear
did_model_nc_2_matched <- lm(NC ~ Rule + Policy + Rule:Policy + Farm_size + Farms_n + Forest_pct + Pasture + SHDI +
                               Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_diffdiff, weights = weights)
summary(did_model_nc_2_matched)
BIC(did_model_nc_2_matched)

did_model_nc_2_matched_robust <- coeftest(did_model_nc_2_matched, vcov = vcovHC(did_model_nc_2_matched, type = "HC0"))

ks.test(resid(did_model_nc_2_matched), "pnorm", mean=mean(resid(did_model_nc_2_matched)), sd=sd(resid(did_model_nc_2_matched)))

skewness(resid(did_model_nc_2_matched))
kurtosis(resid(did_model_nc_2_matched))
hist(resid(did_model_nc_2_matched))
vif(did_model_nc_2_matched)

pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj,
                                                                  . ~ 1)))
pseudoR2(did_model_nc_2_matched)

# Out----------------------------------------------------------------------------
out <- list(did_model_connectance_1, did_model_connectance_2_robust,
            did_model_nc_1, did_model_nc_2_robust, 
            did_model_connectance_1_matched, did_model_connectance_2_matched_robust, 
            did_model_nc_1_matched, did_model_nc_2_matched_robust)

wordreg(c(out), stars = c(0.001, 0.01, 0.05, 0.1), file = "Diff_in_diffs.doc")
