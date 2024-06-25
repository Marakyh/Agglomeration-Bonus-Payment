## Matching ##

# Setup ----------------------------------------------------------
library(dplyr)
library(report)

# for spatial data
library(sf)

# for hetcor
library(polycor)

# for the matching
library(MatchIt)
library(optmatch)
library(Rglpk)
library(gurobi)

# for love plot
library(cobalt)

# for modeling
library(car)

# for likelihood ratio test
library(lmtest)

# for density plot
library(cobalt)
library(gridExtra)
library(ggplot2)

load("./Data/MG_A.Rda")


# Check the model---------------------------------------------------------

glm <- glm(Rule ~ Farm_size + Farms_n + Forest_pct + Pasture + SHDI +
             Ascent_AIV + Soil_sui + Tree_per_ha, family = binomial(), data = MG_A)

summary(glm)
hist(resid(glm))
vif(glm)


# Variable names for plots----------------------------------------------------

full_names <- c(Forest_pct = "Forest (pct)", SHDI = "Shannon's diversity index (absolute)",
                Tree_per_ha = "Individual trees (per ha)", Ascent_AIV = "Slope (degrees)",
                Farms_n = "Farms (absolute)", Farm_size = "Farm size (ha)",
                Pasture = "Pasture (pct)", Soil_sui = "Soil suitability (absolute)")


# Check different kinds of matching----------------------------------------------------------------

# A) Optimal Full matching (method = 'full')

Match1 <- matchit(Rule ~ Farm_size + Farms_n + Forest_pct + Pasture +
                    SHDI + Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_A, method = "full",
                  distance = "glm", estimand = "ATE")
summary(Match1)
love.plot(Match1, grid = TRUE)

Match2 <- matchit(Rule ~ Farm_size + Farms_n + Forest_pct + Pasture +
                    SHDI + Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_A, method = "full",
                  link = "probit", estimand = "ATE")
summary(Match2)
love.plot(Match2, grid = TRUE)

# B) Generalized full matching (?method_quick)

Match3 <- matchit(Rule ~ Farm_size + Farms_n + Forest_pct + Pasture +
                    SHDI + Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_A, method = "quick",
                  distance = "glm", estimand = "ATE")
summary(Match3)
love.plot(Match3)

# C) Subclassification (method = 'subclass') See ?method_subclass for
# the documentation for matchit() with method = 'subclass'

Match4 <- matchit(Rule ~ Farm_size + Farms_n + Forest_pct + Pasture +
                    SHDI + Ascent_AIV + Soil_sui + Tree_per_ha, data = MG_A, method = "subclass",
                  distance = "glm", estimand = "ATE")
summary(Match4)
love.plot(Match4)
abline(v = 0.1, col = "red", lty = 2)


# D) Profile Matching (method = 'cardinality', solver = gurobi)

Match5 <- matchit(Rule ~ Farms_n + Farm_size + Forest_pct + Tree_per_ha + Pasture +
                    SHDI + Ascent_AIV + Soil_sui, data = MG_A, method = "cardinality",
                  estimand = "ATE", ratio = NA, verbose = TRUE, tols = 0.1, time = 120,
                  solver = "gurobi")
summary(Match5)

tiff("Love_plot.tiff", units = "in", width = 10, height = 5, res = 300)
love_plot <- love.plot(Match5, thresholds = c(m = 0.1), var.order = "unadjusted",
                       var.names = full_names, colors = c("#8dd3c7", "#bebada"), sample.names = c("Unmatched",
                                                                                                  "Matched"))
print(love_plot)
dev.off()

# Density plot Vector
variable_names <- c("Farms_n","Farm_size","Forest_pct","Tree_per_ha","Pasture",
                      "SHDI","Ascent_AIV","Soil_sui")

# List to store the ggplot objects for each balance plot
plot_list <- lapply(variable_names, function(var_name) {
  p <- bal.plot(Match5, var.name = var_name, which = "both", sample.names = c("Unmatched",
                                                                              "Matched"))
  p <- p + theme(plot.title = element_blank()) + labs(x = full_names[var_name]) +
    scale_fill_discrete(name = "Treatment \nstatus") + scale_fill_manual(values = c("#bebada",
                                                                                    "#8dd3c7"))
  return(p)
})

# Number of columns in the grid arrangement
num_cols <- 2

# Convert the ggplot objects to 'grobs' using ggplotGrob()
grob_list <- lapply(plot_list, ggplotGrob)

pdf("plot_density_all_portrait.pdf", width = 8.5, height = 11)

grid_arranged_plots <- do.call(grid.arrange, c(grob_list, ncol = num_cols))
print(grid_arranged_plots)

dev.off()


# Save data
MG_matched <- match.data(Match5)
save(MG_matched, file = "./Data/MG_matched.Rda")

# -> go to 05_Estimation

# Cite packages -----------------------------------------------------------

library(NCmisc)
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)
library(report)
cite_packages()