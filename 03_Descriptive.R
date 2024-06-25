## Descriptive Analysis ##

# Setup-------------------------------------------------------------------

library(tidyverse)
# for summary statistics
library(AggregateR)
# for histograms
library(Hmisc)
# for corr.plot
library(corrplot)
# for summary statistics
library(gtsummary)
# for geocomputation
library(sf)
# for correlation
library(heatmaply)

# Load data
# ---------------------------------------------------------------

load("./Data/NBP_A.Rda")
load("./Data/MG_A.Rda")

# Number of municipalities, number of control and treatment---------------
length(MG_A$Municipality_ID)
length(unique(MG_A$Municipality_ID))
table(MG_A$Rule)

a <- as.data.frame.matrix(table(MG_A$Municipality_ID, MG_A$Rule))
a$Sum <- a$`0` + a$`1`
dim(a[a$Sum == 2, ])

# Site types per preservation (control) and networking (treatment) area -------- 
#Table 1

Municipalities_analysis <- unique(MG_A$Municipality_ID)

NBP_A <- NBP_A[NBP_A$Municipality_ID %in% Municipalities_analysis, ]
tab <- table(NBP_A$Site_types_name, NBP_A$Rule)
tab
prop.table(tab, margin = 2)

# Summary statistics of the data set-------------------------------------
# Table 2

MG_A_descr <- subset(MG_A, select = c(Rule, Connectance, NC, Tree_per_ha, Forest_pct,
                                      SHDI, Ascent_AIV, Farms_n, Farm_size, Soil_sui,
                                      Area_ha, Pasture))

cr_table_bM <- MG_A_descr |>
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd}) {min} {max}"),
  )


# Summary statistics of the data set by control and treatment
MG_A_descr <- subset(MG_A, select = c(Rule, Connectance, NC, Tree_per_ha, Forest_pct,
                                      SHDI, Ascent_AIV, Farms_n, Farm_size, Soil_sui,
                                      Pasture))

cr_table_bM <- MG_A_descr |>
  tbl_summary(by = Rule, statistic = list(all_continuous() ~ "{mean} ({sd}) {min} {max}"))

print(cr_table_bM)
summary(MG_A$Rule)

# Histogram Connectance------------------------------------------------

tiff("Histogram_Connectance.tiff", units = "in", width = 5, height = 5,
     res = 300)

Hist.Con <- MG_A |>
  mutate(Rule = factor(Rule, labels = c("Control", "Treatment"))) |>
  ggplot(aes(x = Connectance, fill = Rule)) + geom_histogram(alpha = 0.55,
                                                             position = "identity", bins = 15) + labs(title = "Histogram Connectance",
                                                                                                      x = "Connectance", y = "Count") + labs(fill = "Rule") + scale_fill_manual(values = c(Control = "#bebada",
                                                                                                                                                                                           Treatment = "#8dd3c7"))

print(Hist.Con)

dev.off()

wilcox.test(MG_A[MG_A$Rule == 0, ]$Connectance, MG_A[MG_A$Rule == 1, ]$Connectance,
            alternative = "two.sided")

# Histogram NC----------------------------------------------------------
tiff("Histogram_NC.tiff", units = "in", width = 5, height = 5, res = 300)

Hist.NC <- MG_A |>
  mutate(Rule = factor(Rule, labels = c("Control", "Treatment"))) |>
  ggplot(aes(x = NC, fill = Rule)) + geom_histogram(alpha = 0.55, position = "identity",
                                                    bins = 15) + labs(title = "Histogram Neighbor Connectivity", x = "Neighbor Connectivity",
                                                                      y = "Count") + labs(fill = "Rule") + scale_fill_manual(values = c(Control = "#bebada",
                                                                                                                                        Treatment = "#8dd3c7"))

print(Hist.NC)

dev.off()

t.test(MG_A[MG_A$Rule == 0, ]$NC, MG_A[MG_A$Rule == 1, ]$NC)

# Outliers----------------------------------------------------------------

pdf("Boxplots.pdf", width = 8, height = 6)
# Boxplot
for (col in names(MG_A_descr)) {
  boxplot(MG_A_descr[[col]], main = col)
}
dev.off()

# Correlation-------------------------------------------------------------

Model_variables <- subset(MG_A, select = c(Connectance, NC, Tree_per_ha, Forest_pct,
                                           SHDI, Ascent_AIV, Farms_n, Farm_size, Soil_sui,
                                           Area_ha, Pasture))

cor_matrix <- round(cor(Model_variables, use = "complete.obs"), 2)

heatmap <- heatmaply_cor(x = cor_matrix, xlab = "Features", ylab = "Features",
                         k_col = 2, k_row = 2, dendrogram = "none", cellnote = cor_matrix, cellnote_size = 6)

print(heatmap)

# -> go to 04_Matching

# Cite packages----------------------------------------------------------

library(NCmisc)
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)
library(report)
cite_packages()(heatmap)