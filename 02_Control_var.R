## Prepare Control Variables ##

# Setup ----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(readxl)
# for spatial data
library(sf)
library(stars)
library(rmapshaper)
# for rasterization
library(terra)
library(raster)
# for SHDI
library(landscapemetrics)

 
# Load data ---------------------------------------------------------------

load("./Data/MG.Rda")
load("./Data/NBP_A.Rda")

# Full dataset of BPS
BPS = st_read("./Data/BFF1_BFF2_Vernetzung_Fläche.shp")

# BFS Agricultural data
Agric <- read_excel("./Data/LN_Kennzahlen.xlsx", sheet = "Landwirtschaft")

# Slope
dhm25_epsg21781 <- raster("./Data/dhm25/dhm25_grid_raster")

# Trees Einzelbaum_West =
# st_read('./Data/TLM_BB/swissTLM3D_TLM_EINZELBAUM_GEBUESCH_WEST.shp')
# save(Einzelbaum_West, file = './Data/Einzelbaum_West.Rda')
load("./Data/Einzelbaum_West.Rda")

# Forest patches Bodenbedeckung_West =
# st_read('./Data/TLM_BB/swissTLM3D_TLM_BODENBEDECKUNG_WEST.shp')
# save(Bodenbedeckung_West, file = './Data/Bodenbedeckung_West.Rda')
load("./Data/Bodenbedeckung_West.Rda")

# Arealstatistik
AS <- read.csv("./Data/ag-b-00.03-37-area-csv.csv", header = TRUE, sep = ";")

# Swiss Soil Suitability Map
Soil <- st_read("./Data/Bodeneignungskarte/Bodeneignungskarte_LV95/Bodeneignungskarte_LV95.shp")

# Subset canton of Bern to the districts under investigation

#!remove, is now in 01
MG <- MG[MG$District_ID %in% c(242:246), ]


# Farms------------------------------------------------------------------

Farms <- NBP_A |>
  group_by(Municipality_rule) |>
  summarise(Farms_n = n_distinct(Farm_ID)) |>
  st_drop_geometry()

MG <- left_join(MG, Farms, by = "Municipality_rule", keep = FALSE)

# Farm Size---------------------------------------------------------

#Subset to specific areas of intervention
BPS = BPS[BPS$Massnahmen %in% c("Vernetzungsgebiet Hügel / Hang","Erhaltungsgebiet strukturreiche Landschaft","Vernetzungsgebiet Tal / offenes Agrarland", "Vernetzungsgebiet offene Wiesenlandschaft"),]

#Divide into control and treatment
With_rule = c('Vernetzungsgebiet Tal / offenes Agrarland','Vernetzungsgebiet Hügel / Hang','Vernetzungsgebiet offene Wiesenlandschaft')
BPS$Rule = as.factor(ifelse(BPS$Massnahmen %in% With_rule, 1, 0))

#Add ID
BPS$Municipality_rule <- paste(BPS$Gemeinde_I, BPS$Rule)

Farms_BPS <- BPS |>
  group_by(Municipality_rule) |>
  summarise(Farms_BPS_n = n_distinct(Betriebs_I)) |>
  st_drop_geometry()

MG <- left_join(MG, Farms_BPS, by = "Municipality_rule", keep = FALSE)
#divide area to farms after adding up the data set (below)

# Area---------------------------------------------------------------

# Check if Shape_Area is the correct area of the AIVs
all(round(as.numeric(st_area(MG$geometry)) - MG$Shape_Area, 2) == 0)

# Size of AIV is in m2, we need to divide this number by 10'000 to get hectares
MG$Area_ha <- MG$Shape_Area/10000


# Acre ---------------------------------------------------------

Agric$Municipality_ID <- as.numeric(gsub(".*?([0-9]+).*", "\\1", Agric$BFS))

# Add control variables to MG data set:
MG <- left_join(MG, Agric[, -1], by = "Municipality_ID", keep = FALSE)


# Pasture-----------------------------------------------------------------

Pasture <- NBP_A |>
  st_drop_geometry() |>
  dplyr::group_by(NBP_A$Municipality_ID, NBP_A$Rule) |>
  dplyr::summarise(share = (sum(Site_types_ID == "617")/n()))

# Test if the share is correct
a <- NBP_A[NBP_A$Municipality_ID == "301" & NBP_A$Site_types_ID == "617" &
             NBP_A$Rule == 0, ]
b <- NBP_A[NBP_A$Municipality_ID == "301" & NBP_A$Rule == 0, ]
Pasture[Pasture$`NBP_A$Municipality_ID` == "301" & Pasture$`NBP_A$Rule` ==
          0, ]$share == dim(a)[1]/dim(b)[1]

Pasture$Municipality_rule <- paste(Pasture$"NBP_A$Municipality_ID", Pasture$"NBP_A$Rule")
Pasture <- rename(Pasture, Pasture = share)
MG <- left_join(MG, Pasture[c("Municipality_rule", "Pasture")], by = "Municipality_rule",
                keep = FALSE)


# Slope-------------------------------------------------------------------
# Slope doesn't have any crs in the metadata. we know that it should
# be +init=epsg:21781

proj4string(dhm25_epsg21781) <- crs("+init=epsg:21781")
dhm25_epsg21781

# Check if CRS is correct
Bern_extent <- st_bbox(MG, crs = crs("+init=epsg:2056"))
Bern_extent_epsg21781 <- Bern_extent |>
  st_as_sfc() |>
  st_transform(crs = crs("+init=epsg:21781")) |>
  st_bbox()

# Crop the map to the canton of Bern (more efficient)
dhm25_Bern <- crop(dhm25_epsg21781, Bern_extent_epsg21781)
dhm25 <- projectRaster(dhm25_Bern, crs = crs("+init=epsg:2056"))

compareCRS(MG, dhm25)

slope <- terrain(dhm25, opt = c("slope"), unit = "degrees")

MG$Ascent_AIV <- extract(slope, MG, fun = mean, na.rm = TRUE)


# Trees-------------------------------------------------------------------

# Get rid of 3D and adjust CRS accordingly
Einzelbaum <- Einzelbaum_West |>
  st_zm() |>
  st_set_crs("EPSG:2056")

# Calculate the number of trees in each area of intervention
MG$Tree_n <- lengths(st_intersects(MG, Einzelbaum))


# Forest patches-------------------------------------------------------------------

# Get rid of umlauts
Bodenbedeckung_West$OBJEKTART <- iconv(Bodenbedeckung_West$OBJEKTART, from = "UTF-8",
                                       to = "ASCII//TRANSLIT")

# Select all 'Wald' and 'Gehoelzflaeche' objects from the
# 'Bodenbedeckung_CH' data set and convert to a simple feature (sf)
# object and set coordinate reference system (CRS)
Wald = Bodenbedeckung_West[Bodenbedeckung_West$OBJEKTART %in% c("Wald","Wald offen","Gehoelzflaeche"),] %>% st_zm() %>% st_set_crs(., "EPSG:2056")

# Optionally, check output in QGIS st_write(Wald, './Wald.shp')

# Add ID for summing up intersection area
MG$ID <- c(1:dim(MG)[1])

# Find the intersection between the 'MG' and 'Wald' data sets
# create new column with area of the intersection, Reduce data set to
# the relevant columns, drop geometry as we don't need it

intersect <- st_intersection(MG, Wald) %>% 
  dplyr::mutate(intersect_area = st_area(.)) %>%  
  dplyr::select(ID, intersect_area) %>%   
  st_drop_geometry() 

# Group the intersection data set by feature ID, Calculate the total
# intersection area for each feature ID

intersect <- intersect %>%
  dplyr::group_by(ID) %>% 
  dplyr::summarize(MG_Forest = sum(intersect_area)) 

MG <- left_join(MG, intersect, by = "ID")

MG$Forest_pct <- as.numeric(MG$MG_Forest/MG$Shape_Area)
MG$Forest_pct <- replace(MG$Forest_pct, is.na(MG$Forest_pct), 0)  #if there is no forest, value should be 0 and not NA


# Shannon Index----------------------------------------

# Subset 'Arealstatistik' to the canton of Bern
AS_Bern <- AS[AS$GMDE %in% unique(MG$Municipality_ID), ]
# We use AS97_72 (72 Grundkategorien gemäss Standardnomenklatur der
# Arealstatistik 1992/97)
AS97_72 <- AS_Bern[, c("E", "N", "AS97_72")]

# Rasterize
AS97_72_raster <- rasterFromXYZ(AS97_72, crs = crs(MG))

# Check if extent is comparable
extent(AS97_72_raster)
extent(MG)

# Enlarge map with the surrounding cantons (for a better results at
# the border) (Cantons Luzern, Freiburg, Solothurn, Neuchatel - see
# https://www.bfs.admin.ch/bfs/de/home/grundlagen/agvch.assetdetail.23886073.html
# for municipality ID)
Municipalities_bern_surroundings <- c(301:995, 1001:1151, 2008:2338, 2401:2622,
                                      6404:6512)

AS97_72_raster <- AS[AS$GMDE %in% Municipalities_bern_surroundings, c("E",
                                                                      "N", "AS97_72")] |>
  rasterFromXYZ()

# Set CRS
crs(AS97_72_raster) <- crs(MG)

# Check in QGIS if CRS is set correctly writeRaster(AS97_72_raster,
# filename = './Test/AS97_72_raster.tif', format = 'GTiff',
# overwrite=TRUE)

# Calculate SHDI for the different areas of intervention Ensure if
# input data is valid
check_landscape(AS97_72_raster)

buffered_polygons <- st_buffer(MG, dist = 25)

MG_shdis <- sample_lsm(AS97_72_raster, buffered_polygons, what = "lsm_l_shdi")

MG <- cbind(MG, MG_shdis[, c("value")])
MG <- rename(MG, SHDI = value)

# Soil Suitability ------------------------------

# Replace 99 values with 0
Soil$KU_CODE[Soil$KU_CODE == 99] <- 0

# Perform a spatial join between MG and Soil datasets to associate Soil attributes with MG observations
MG_with_Soil <- st_join(MG, Soil["KU_CODE"])

# Aggregate the K values by each MG observation and calculate the mean
MG_mean_K <- MG_with_Soil %>%
  group_by(ID) %>% 
  summarise(Soil_sui = mean(KU_CODE, na.rm = TRUE)) %>%
  st_drop_geometry()

MG <- left_join(MG, MG_mean_K, by = "ID")

# Sites per AIV------------------------------------

Sites <- NBP_A %>%
  group_by(Municipality_rule) %>%
  summarise(count = n()) %>%
  st_drop_geometry()

MG <- left_join(MG, Sites, by = "Municipality_rule", keep = FALSE)
MG <- rename(MG, Sites_n = count)

# Wrap up-----------------------------------------------------------------

# Test: Compare shapefile MG_out_02 with MG from the beginning in
# QGIS st_write(MG, './Test/MG.shp')

MG <- st_drop_geometry(MG)

# Convert cat. variables to numeric variables:
head(MG$Rule)
MG$Rule <- as.numeric(MG$Rule) - 1
head(MG$Rule)

# Calculate mean or sum up
MG_A <- MG %>%
  group_by(Municipality_rule) %>%
  summarise(
    across(-c(Shape_Area, Area_ha, Tree_n), ~mean(., na.rm = TRUE)),
    across(c(Shape_Area, Area_ha, Tree_n), ~sum(., na.rm = TRUE))
  )

# Calculate the relative variables
MG_A$Farm_size <- MG_A$Area_ha/MG_A$Farms_BPS_n
MG_A$Tree_per_ha <- MG_A$Tree_n/MG_A$Area_ha


# Check dataset
head(MG_A)

# Remove variables that will not be used anymore
columns_to_remove <- c("District_ID", "ID", "MG_Forest")
MG_A <- MG_A[, !(names(MG_A) %in% columns_to_remove)]

# Round certain columns
# Define columns to exclude from rounding
columns_to_exclude <- c("Municipality_rule", "Connectance", "NC")

# Round all columns except those in columns_to_exclude
MG_A[, !(names(MG_A) %in% columns_to_exclude)] <- round(MG_A[, !(names(MG_A) %in% columns_to_exclude)], 5)

# Check for missing values
Missing <- MG_A[!complete.cases(MG_A), ]

# Exclude observations with only 1 farm
Farms_one <- unique(MG_A[MG_A$Farms_n==1,]$Municipality_rule)
length(Farms_one)
MG_A <- MG_A[!(MG_A$Municipality_rule %in% Farms_one), ]

# Save the dataset
save(MG_A, file = "./Data/MG_A.Rda")

# -> go to 03_Descriptive

# Cite packages -----------------------------------------------------------

library(NCmisc)
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)
library(report)
cite_packages()