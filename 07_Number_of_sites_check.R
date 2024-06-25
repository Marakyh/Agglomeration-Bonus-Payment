## Compare the number of sites in treatment and contol area ##

# Setup ----------------------------------------------------------

library(tidyverse)
library(sf)
library(readxl)

options(scipen=999) 

#Districts
Districts = read.csv("./Data/GemeindenKtBE.csv")

#Geo-referenced communal borders
Municipality_georeferenced = st_read("./Data/GRENZ5_G5.shp")

#Biodiversity promotion sites BPS
BPS = st_read("./Data/BFF1_BFF2_Vernetzung_Fläche.shp")
# CRS: CH1903+ / LV95

#Areas of intervention
MG = st_read("./Data/VERNETZ_MG.shp")

# Sites -------------------------------------------------------------------

#Exclude NBP-Sites outside the canton of Bern
BPS = BPS[!BPS$Gemeinde_I >= 1000 & !BPS$Gemeinde_I <301,]

#Only inlcude sites in specific areas of intervention
BPS = BPS[BPS$Massnahmen %in% c("Vernetzungsgebiet Hügel / Hang","Erhaltungsgebiet strukturreiche Landschaft","Vernetzungsgebiet Tal / offenes Agrarland", "Vernetzungsgebiet offene Wiesenlandschaft"),]

#Divide into control and treatment
With_rule = c('Vernetzungsgebiet Tal / offenes Agrarland','Vernetzungsgebiet Hügel / Hang','Vernetzungsgebiet offene Wiesenlandschaft')
BPS$Rule = as.numeric(ifelse(BPS$Massnahmen %in% With_rule, 1, 0))

#Add identifier
BPS$Municipality_rule <- paste(BPS$Gemeinde_I, BPS$Rule)

#Only include observations from the final analysis 
load("./Data/MG_A.Rda")
Municipalities_analysis <- unique(MG_A$Municipality_rule)
BPS <- BPS[BPS$Municipality_rule %in% Municipalities_analysis,]

#Add the number of farms
BPS <- left_join(BPS, MG_A[,c("Municipality_rule","Farms_n")], by = c("Municipality_rule" = "Municipality_rule"), keep = FALSE)

BPS <- st_drop_geometry(BPS)

save(BPS, file = "./Data/BPS.Rda")

# Pasture and meadows smaller than 0.3 ha----------------------------------------
BPS_pasture_meadows_NBP_small= BPS[BPS$CODE %in% c("611","612","617") & BPS$VERNETZUNG < 30 & BPS$VERNETZUNG > 0,]

#Sum up per Municipality_rule
Check0 <- BPS_pasture_meadows_NBP_small %>%
  group_by(Municipality_rule) %>%
  summarise(
    across(c(Farms_n, Rule), ~mean(., na.rm = TRUE)),
    Sites_number = n(),
    Farms_specific = n_distinct(Betriebs_I)
  )

Check0$sites_per_farm <- Check0$Sites_number / Check0$Farms_specific

#as the variables are clearly not normally distributed, we use wilcox.test is thus appropriate
wilcox.test(Check0[Check0$Rule == 0, ]$sites_per_farm, Check0[Check0$Rule == 1, ]$sites_per_farm,
            alternative = "two.sided")

# Pasture and meadows larger than 0.3 ha----------------------------------------
BPS_pasture_meadows_NBP_large = BPS[BPS$CODE %in% c("611","612","617") & BPS$VERNETZUNG > 30,]

#Sum up per Municipality_rule
Check1 <- BPS_pasture_meadows_NBP_large %>%
  group_by(Municipality_rule) %>%
  summarise(
    across(c(Farms_n, Rule), ~mean(., na.rm = TRUE)),
    Sites_number = n(),
    Farms_specific = n_distinct(Betriebs_I)
  )

Check1$sites_per_farm <- Check1$Sites_number / Check1$Farms_specific

#as the variables are clearly not normally distributed, we use wilcox.test is thus appropriate
wilcox.test(Check1[Check1$Rule == 0, ]$sites_per_farm, Check1[Check1$Rule == 1, ]$sites_per_farm,
            alternative = "two.sided")









