## Estimation Model ##

# Setup ----------------------------------------------------------
library(tidyverse)
library(moments)
# library(sf)
library(rcompanion)
# for tobit
library(AER)
# for regression output
library(texreg)

load("./Data/MG_matched.Rda")


# Transform connectance-----------------------------------------------------
hist(MG_matched$Connectance)
MG_matched$ConnectanceL <- log(MG_matched$Connectance)
hist(MG_matched$ConnectanceL)
#Set the NA values according to Cameron and Trivedi (2022) 
summary_lny <- summary(MG_matched$ConnectanceL[is.finite(MG_matched$ConnectanceL)])
min_value <- summary_lny[1]
gamma <- min_value - 1e-06
MG_matched$ConnectanceL <- ifelse(MG_matched$ConnectanceL == -Inf, gamma,
                                  MG_matched$ConnectanceL)

Connectance_min <- min(MG_matched$ConnectanceL)  #LowerLimit

# Check normality
skewness(MG_matched$ConnectanceL)  # (normal: 0)
kurtosis(MG_matched$ConnectanceL)  # (normal: 3)
ks.test(MG_matched$ConnectanceL, "pnorm", mean = mean(MG_matched$ConnectanceL),
        sd = sd(MG_matched$ConnectanceL))

plotNormalHistogram(MG_matched$ConnectanceL, prob = FALSE, main = "Normal Distribution overlay on Histogram")

# Transform other variables-----------------------------------------------
MG_matched$Farms_nR <- 1/MG_matched$Farms_n

head(MG_matched$Rule)
MG_matched$Rule <- factor(MG_matched$Rule, levels = 0:1, labels = c("Control",
                                                                    "Treatment"))
head(MG_matched$Rule)

# Dataset without outliers------------------------------------------------------
MG_matched_subset <- subset(MG_matched, select = c(ConnectanceL, NC, Tree_per_ha,
                                                   Forest_pct, SHDI, Ascent_AIV, 
                                                   Farms_n, Farm_size, Pasture,
                                                   Soil_sui))

for (col in names(MG_matched_subset)) {
  boxplot(MG_matched_subset[[col]], main = col)
}

top_1_Farms_n <- tail(order(MG_matched$Farms_n), 1)
top_2_Forest_pct <- tail(order(MG_matched$Forest_pct), 2)
top_2_Tree_per_ha <- tail(order(MG_matched$Tree_per_ha), 2)
top_1_Pasture <- tail(order(MG_matched$Pasture), 2)

top_indices <- Reduce(union, list(top_1_Farms_n, top_2_Forest_pct, top_2_Tree_per_ha, top_1_Pasture))

MG_matched_without_outliers <- MG_matched[-top_indices, ]


# Tobit regressions-------------------------------------------------------

# Connectance-------------------------------------------------------------

tobit1 <- tobit(ConnectanceL ~ Rule + 
                  Forest_pct + SHDI + Tree_per_ha + Ascent_AIV + Farms_nR + 
                  Farm_size + Pasture + Soil_sui, left = Connectance_min,
                data = MG_matched, robust = TRUE, weights = weights, x=TRUE)
tobit2 <- tobit(ConnectanceL ~ Rule + Forest_pct + SHDI + Tree_per_ha +
                  Farms_nR + Farm_size + Pasture + Soil_sui, left = Connectance_min,
                data = MG_matched, robust = TRUE, weights = weights, x=TRUE)
tobit3 <- tobit(ConnectanceL ~ Rule + Forest_pct + SHDI + Tree_per_ha +
                  Farms_nR + Farm_size + Pasture + Soil_sui, left = Connectance_min,
                data = MG_matched_without_outliers, robust = TRUE, weights = weights, x=TRUE)
#No log
MG_matched$Connectance <- MG_matched$Connectance*100
tobit4 <- tobit(Connectance ~ Rule + 
                  Forest_pct + SHDI + Tree_per_ha + Ascent_AIV + Farms_nR + 
                  Farm_size + Pasture + Soil_sui, left = 0,
                data = MG_matched, robust = TRUE, weights = weights)
#No inverted farms
tobit5 <- tobit(ConnectanceL ~ Rule + 
                  Forest_pct + SHDI + Tree_per_ha + Ascent_AIV + Farms_n + 
                  Farm_size + Pasture + Soil_sui, left = Connectance_min,
                data = MG_matched, robust = TRUE, weights = weights)



# Check model output
tobits_Con <- list(tobit1, tobit2, tobit3, tobit4, tobit5)

screenreg(tobits_Con)

lapply(tobits_Con, vif)
lapply(tobits_Con, resid) |>
  lapply(mean)
lapply(tobits_Con, resid) |>
  lapply(skewness)
lapply(tobits_Con, resid) |>
  lapply(kurtosis)

ks.test(resid(tobit1), "pnorm", mean = mean(resid(tobit1)), sd = sd(resid(tobit1)))
ks.test(resid(tobit2), "pnorm", mean = mean(resid(tobit2)), sd = sd(resid(tobit2)))
ks.test(resid(tobit3), "pnorm", mean = mean(resid(tobit3)), sd = sd(resid(tobit3)))

plotNormalHistogram(resid(tobit1), prob = FALSE, main = "Normal Distribution overlay on Histogram")

fit_tobit_1 <- predict(tobit1, type = "response")
plot(fit_tobit_1, resid(tobit1))

fit_tobit_1 <- predict(tobit1, type = "response")
residuals <- resid(tobit1)

plot(fit_tobit_1, residuals, xlab = "Fitted values", ylab = "Standardized residuals",
     main = "Tukey-Anscombe Plot")
abline(h = 0, col = "red")  # Add a horizontal line at y = 0

#Censored mean in levels, following Cameron and Trivedi, keeping all variables constant except for "rule"
MG_matched_without_outliers$Rule2 <- as.numeric(MG_matched_without_outliers$Rule)-1

mu <- sum(apply(tobit3$x,2,FUN=median)[-2] * tobit3$coef[-2]) + tobit3$coef[2]*MG_matched_without_outliers$Rule2
sigma <- tobit3$scale
ey0 <- exp(mu + sigma^2/2)
p0 <- 1 - (pnorm((gamma-mu-sigma^2)/sigma))
ey <- ey0*p0

e <- as.data.frame(cbind(MG_matched_without_outliers$Rule, ey))

conditional_mean <- e %>%
  group_by(V1) %>%
  summarise(mean_outcome = mean(ey, na.rm = TRUE))

conditional_mean

#Adjusted estimates
mu <- sum(apply(tobit3$x,2,FUN=median) * tobit3$coef)
p0 <- 1 - (pnorm((gamma-mu-sigma^2)/sigma))
p0 * tobit3$coef[-1]


# Neighbor Connectivity---------------------------------------------------

hist(MG_matched$NC)
skewness(MG_matched$NC)  
kurtosis(MG_matched$NC)  
ks.test(MG_matched$NC, "pnorm", mean = mean(MG_matched$NC),sd = sd(MG_matched$NC))

mean(MG_matched$NC)

tobit1 <- tobit(NC ~ Rule + Forest_pct + SHDI + Tree_per_ha + Ascent_AIV +
                  Farms_n + Farm_size + Pasture + Soil_sui, left = 0, right = 100,
                data = MG_matched, robust = TRUE, weights = weights, x=TRUE)
tobit2 <- tobit(NC ~ Rule + Forest_pct + SHDI + Tree_per_ha + 
                  Farms_n + Farm_size + Pasture + Soil_sui, left = 0, right = 100, data = MG_matched,
                robust = TRUE, weights = weights, x=TRUE)
tobit3 <- tobit(NC ~ Rule + Forest_pct + SHDI + Tree_per_ha + 
                  Farms_n + Farm_size + Pasture + Soil_sui, left = 0, right = 100, data = MG_matched_without_outliers,
                robust = TRUE, weights = weights, x=TRUE)
tobit4 <- tobit(NC ~ Rule + SHDI + Tree_per_ha + 
                  Farms_n + Farm_size + Pasture + Soil_sui, left = 0, right = 1, data = MG_matched_without_outliers,
                robust = TRUE, weights = weights, x=TRUE)


# Check model output
tobits_NC <- list(tobit1, tobit2, tobit3, tobit4)

screenreg(tobits_NC)

lapply(tobits_NC, resid) |>
  lapply(mean)
lapply(tobits_NC, resid) |>
  lapply(skewness)
lapply(tobits_NC, resid) |>
  lapply(kurtosis)
lapply(tobits_NC, vif)

ks.test(resid(tobit1), "pnorm", mean = mean(resid(tobit1)), sd = sd(resid(tobit1)))
ks.test(resid(tobit2), "pnorm", mean = mean(resid(tobit2)), sd = sd(resid(tobit2)))
ks.test(resid(tobit3), "pnorm", mean = mean(resid(tobit3)), sd = sd(resid(tobit3)))
ks.test(resid(tobit3), "pnorm", mean = mean(resid(tobit4)), sd = sd(resid(tobit4)))

plotNormalHistogram(resid(tobit1), prob = FALSE, main = "Normal Distribution overlay on Histogram")


#Marginal effects on y 
mu <- sum(apply(tobit4$x,2,FUN=median) * tobit4$coef)
sigma <- tobit4$scale

p0 <- pnorm((1-mu)/sigma) - pnorm((0-mu)/sigma)

p0 * tobit4$coef[-1]


# Export results--------------------------------------------------------------

tobits <- c(tobits_Con, tobits_NC)

wordreg(tobits, stars = c(0.001, 0.01, 0.05, 0.1), file = "Tobits.doc")

# Cite packages-----------------------------------------------------------

library(NCmisc)
list.functions.in.file(rstudioapi::getSourceEditorContext()$path, alphabetic = TRUE)
library(report)
cite_packages()