## Prepare Outcome Variables ##

# Setup ----------------------------------------------------------

library(tidyverse)
library(sf)
library(spdep)
library(AggregateR)
library(rlist)
library(gtsummary)
library(readxl)
library(data.table)
# for rasterization
library(stars)
library(terra)
library(raster)
library(rmapshaper)
# for connectivity metrics
library(sfnetworks)
library(tidygraph)
library(igraph)
# no scientific notation
options(scipen = 999)

# Load functions
# ----------------------------------------------------------

source("./Functions/NeighborCon_function.R")

# Abbreviations
# -----------------------------------------------------------

# NBP = network bonus payment NBP_A = network bonus payment data set
# for analysis BPS = biodiversity promotion sites ABP = agglomeration
# bonus payment MG = areas of intervention (dt. Massnahmengebiet, MG)

# Load data
# ---------------------------------------------------------------

# Districts
Districts <- read.csv("./Data/GemeindenKtBE.csv")

# Geo-referenced communal borders
Municipality_georeferenced <- st_read("./Data/GRENZ5_G5.shp")

# Biodiversity promotion sites BPS
BPS <- st_read("./Data/BFF1_BFF2_Vernetzung_Fläche.shp")
# CRS: CH1903+ / LV95

# Areas of intervention
MG <- st_read("./Data/VERNETZ_MG.shp")


# Sites
# -------------------------------------------------------------------

# Only inlcude sites with a network bonus payment NBP
NBP <- filter(BPS, BPS$VERNETZUNG > 0)

# Only inlcude sites in specific areas of intervention
NBP <- NBP[NBP$Massnahmen %in% c("Vernetzungsgebiet Hügel / Hang", "Erhaltungsgebiet strukturreiche Landschaft",
                                 "Vernetzungsgebiet Tal / offenes Agrarland", "Vernetzungsgebiet offene Wiesenlandschaft"),
]

# Only inlcude sites 611,612,617 smaller than 0.3 ha
NBP <- NBP[NBP$CODE %in% c("611", "612", "617") & NBP$VERNETZUNG < 30,
]

# Divide into control and treatment
With_rule <- c("Vernetzungsgebiet Tal / offenes Agrarland", "Vernetzungsgebiet Hügel / Hang",
               "Vernetzungsgebiet offene Wiesenlandschaft")
NBP$Rule <- as.factor(ifelse(NBP$Massnahmen %in% With_rule, 1, 0))

# Test
table(NBP$Rule)[[2]] == dim(filter(NBP, Massnahmen %in% With_rule))[1]

# Adding information on districts to the dataset
NBP <- left_join(NBP, Districts, by = c(Gemeinde_I = "BFS"), keep = FALSE)

# Exclude NBP-Sites outside the canton of Bern
NBP <- NBP[!NBP$Gemeinde_I >= 1000 & !NBP$Gemeinde_I < 301, ]

# Exclude sites where treatment and control overlap
NBP$NBP_ID <- c(1:dim(NBP)[1])
NBP0 <- NBP[NBP$Rule == 0, ]
NBP1 <- NBP[NBP$Rule == 1, ]
inter_NBP <- st_intersection(NBP0, NBP1, model = "closed")
inter_NBP_ID <- sort(rbind(inter_NBP$NBP_ID, inter_NBP$NBP_ID.1))
NBP2 <- NBP[!NBP$NBP_ID %in% inter_NBP_ID, ]

# Add an ID combining municipality and treatment status
NBP$Municipality_rule <- paste(NBP$Gemeinde_I, NBP$Rule)

# Exclude unused rows, rename
NBP_A <- NBP |>
  select(Gemeinde_I, Betriebs_I, VK_Nr, CODE, Kulturname, Massnahmen,
         Rule, Municipality_rule) |>
  rename(Municipality_ID = Gemeinde_I, Farm_ID = Betriebs_I, District_ID = VK_Nr,
         Site_types_ID = CODE, Site_types_name = Kulturname, Area_of_intervention = Massnahmen)

# Save the dataset
save(NBP_A, file = "./Data/NBP_A.Rda")

#------------Areas of intervention (dt. Massnahmengebiete MG)

# Subset AIVs
MG <- MG[MG$MG_TYP %in% c(2, 3, 6, 9), ]

# Subset into treatment and control areas of intervention (treatment
# = type 2 (Vernetzungsgebiet Tal/offenes Agrarland), 3
# (Vernetzungsgebiet Huegel/Hang), 6 (Vernetzungsgebiet offene
# Wiesenlandschaft)
With_rule <- c(2, 3, 6)
MG$Rule <- as.factor(ifelse(MG$MG_TYP %in% With_rule, 1, 0))

# Add information on municipality and district to the data set
Municipality_georeferenced <- Municipality_georeferenced[, c(2:3)]
MG <- st_join(MG, Municipality_georeferenced, join = st_intersects, left = TRUE,
              largest = TRUE)
MG <- left_join(MG, Districts, by = c(BFSNR = "BFS"), keep = FALSE)

# Exclude lakes
MG <- MG[!MG$BFSNR >= 1000, ]

# Exclude MGs outside of Bern (no matching BFSNR)
MG <- MG[!is.na(MG$VK), ]

# Check for missing values
missing <- MG[rowSums(is.na(MG)) > 0, ]

# Add an ID combining municipality and treatment status
MG$Municipality_rule <- paste(MG$BFSNR, MG$Rule)

# Exclude unused rows, rename
MG <- MG |>
  select(MG_TYP, Shape_Area, Rule, VK_Nr, BFSNR, Municipality_rule) |>
  rename(AIV_Typ = MG_TYP, Municipality_ID = BFSNR, District_ID = VK_Nr)


# Connectivity Metrics
# ----------------------------------------------------

# Compute Gamma Control
NBP_A0 <- NBP_A[NBP_A$Rule == 0, ]
Municipality_ID_sorted <- sort(unique(NBP_A0$Municipality_ID))
Gamma <- numeric(length(Municipality_ID_sorted))

for (i in 1:length(Municipality_ID_sorted)) {
  NBP_Com <- NBP_A0[NBP_A0$Municipality_ID == Municipality_ID_sorted[i],
  ]
  Gamma[i] <- Connectance(NBP_Com)
}

Gamma0 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 0, Gamma))
Gamma0$ID <- paste(Gamma0$Municipality_ID_sorted, Gamma0$Rule)

# Compute Gamma Treatment
NBP_A1 <- NBP_A[NBP_A$Rule == 1, ]
Municipality_ID_sorted <- sort(unique(NBP_A1$Municipality_ID))
Gamma <- numeric(length(Municipality_ID_sorted))

for (i in 1:length(Municipality_ID_sorted)) {
  NBP_Com <- NBP_A1[NBP_A1$Municipality_ID == Municipality_ID_sorted[i],
  ]
  Gamma[i] <- Connectance(NBP_Com)
}

Gamma1 <- as.data.frame(cbind(Municipality_ID_sorted, Rule = 1, Gamma))
Gamma1$ID <- paste(Gamma1$Municipality_ID_sorted, Gamma1$Rule)

Gamma <- rbind(Gamma0, Gamma1)

# Check with several municipalities by hand that gamma is correct
# (382 0, 385 0, 438 0, 590 0, 492 1)

# Check NA values
Missing <- Gamma[!complete.cases(Gamma), ]
Gamma[is.na(Gamma)] <- 0  #the NaN values arise due to 0/0 cases

# Compute Neighbor Connectivity Control
NCon <- by(NBP_A0, NBP_A0$Municipality_ID, NeighborCon)

Municipality_ID_sorted <- sort(unique(NBP_A0$Municipality_ID))
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

# Compute Neighbor Connectivity Treatment
NCon <- by(NBP_A1, NBP_A1$Municipality_ID, NeighborCon)

Municipality_ID_sorted <- sort(unique(NBP_A1$Municipality_ID))
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
Missing <- NC[!complete.cases(NC), ]
NC[is.na(NC)] <- 0 #the NaN values arise due to 0/0 cases

# Merge with MG
MG <- left_join(MG, Gamma[, c("Gamma", "ID")], by = c(Municipality_rule = "ID")) |>
  left_join(NC[, c("NC", "ID")], by = c(Municipality_rule = "ID"))

MG_df <- st_drop_geometry(MG)
Missing <- MG_df[!complete.cases(MG_df), ]

Missing_ID <- unique(Missing$Municipality_rule)

# Subset the NBP data frame to only include rows where Municipality_rule
# matches Missing_ID -> are there no sites that cause the NaN value?
sub_df <- NBP_A[NBP_A$Municipality_rule %in% Missing_ID, ]

nrow(sub_df) == 0

# NA values because they have no sites for this area, thus we exclude them
MG <- MG[!MG$Municipality_rule %in% Missing_ID, ]

# Subset canton of Bern to the districts under investigation
MG <- MG[MG$District_ID %in% c(242:246), ]

#Rename Gamma to Connectance
MG <- rename(MG, Connectance = Gamma)

# Save the dataset
save(MG, file = "./Data/MG.Rda")

# -> go to 02_Control_variables

# Cite packages -----------------------------------------------------------

library(NCmisc)
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)
library(report)
cite_packages()
