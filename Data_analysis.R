### Global -------------------------------------------------------
source("Data_prep.R")

library(glmmTMB)
library(visreg)
library(performance)
library(corrplot)
library(censReg)
library(car)
library(sp)
library(spdep)

### Raw data --------------------------------------------------------------
cm <- cor(df_mod[, c(2,3,5,6, 11:13, 20, 21)], use = "pairwise.complete.obs")
corrplot(cm)

### Models ----------------------------------------------------------------
colnames(df_mod)[c(2:4)] <- c("Pb210", "Cs137", "Am241")

# Censor data Am-241
df_mod$Am241_c <- df_mod$Am241
df_mod$Am241_c[c(1,7,8,9,11,14,15,18)] <- 0

# log-scaling
df_mod$Size_kmq_l <- log(df_mod$Size_kmq)
df_mod$prec_CHB_l <- log(df_mod$prec_CHB)
df_mod$prec_GLF_l <- log(df_mod$prec_GLF)

### R-2 data for Fig. 2
m1_r <- lm(log(Cs137) ~ log(Pb210), data = df_mod); check_model(m1_r)
m2_r <- lm(log(Am_241_h) ~ log(Cs137), data = df_mod); check_model(m2_r)
m3_r <- lm(log(Am_241_h) ~ log(Pb210), data = df_mod); check_model(m2_r)

R_r2 <- data.frame(Pair = c("Cs-137~Pb-210", "Am-241~Cs-137", "Am-241~Pb-210"), 
                   "R.sq" = round(c(summary(m1_r)$adj.r.squared, summary(m2_r)$adj.r.squared, summary(m3_r)$adj.r.squared),3))

### Lead-210 
mP.0 <- lm(log(Pb210) ~ 1, 
           data = df_mod); summary(mP.0)
mP.1 <- lm(log(Pb210) ~ OM, 
           data = df_mod); summary(mP.1)
mP.2 <- lm(log(Pb210) ~ OM + Average_altitude, 
           data = df_mod); summary(mP.2)
mP.3 <- lm(log(Pb210) ~ OM + Size_kmq_l, 
           data = df_mod); summary(mP.3)
mP.4 <- lm(log(Pb210) ~ OM + Size_kmq_l + Average_altitude, 
           data = df_mod); summary(mP.4)

check_model(mP.1); check_outliers(mP.1)# Two detected outliers
check_model(mP.2); check_outliers(mP.2)# One detected outlier
check_model(mP.3); check_outliers(mP.3)# One detected outlier
check_model(mP.4); check_outliers(mP.4)# No outliers detected 

anova(mP.0, mP.1, mP.2, mP.4, mP.3)
AIC(mP.0, mP.1, mP.2, mP.3, mP.4)

summary(mP.4); Anova(mP.4)

### Cs-137
mC.0 <- lm(log(Cs137) ~ 1,
           data = df_mod); summary(mC.0)
mC.1 <- lm(log(Cs137) ~ OM,
           data = df_mod); summary(mC.1)
mC.2 <- lm(log(Cs137) ~ OM + Size_kmq_l,
           data = df_mod); summary(mC.2)
mC.3 <- lm(log(Cs137) ~ OM + prec_GLF_l + prec_CHB_l,
           data = df_mod); summary(mC.3)
mC.4 <- lm(log(Cs137) ~ OM + Size_kmq_l + prec_GLF_l + prec_CHB_l,
           data = df_mod); summary(mC.4)

check_model(mC.1); check_outliers(mC.1)# No outliers detected 
check_model(mC.2); check_outliers(mC.2)# No outliers detected 
check_model(mC.3); check_outliers(mC.3)# No outliers detected 
check_model(mC.4); check_outliers(mC.4)# No outliers detected 

AIC(mC.0, mC.1, mC.2, mC.3, mC.4)
anova(mC.0, mC.1, mC.2, mC.3, mC.4)

summary(mC.2); Anova(mC.2)

### Am-241 / Model - censored regression 
mA.1 <- censReg(log(Am241_c+0.001) ~ OM + Size_kmq_l + prec_GLF_l, 
                data = df_mod); summary(mA.1)
Anova(mA.1, test.statistic = "F")

### Spatial autocorrelation -----------------------------------------------
df_spat <- df_mod

# Create Spatial Points Data Frame
coordinates(df_spat) <- c("Lat", "Lon")
proj4string(df_spat) <- CRS("+proj=longlat +ellps=WGS84")

# Spatial weights using K-nearest neighbors
spatial_weights <- knn2nb(knearneigh(coordinates(df_spat), k = 3))
summary(spatial_weights)

spatial_weights <- nb2listw(spatial_weights, style = "W", zero.policy = TRUE)

# Moran's I for each isotope
SM_P <- moran.test(df_spat$Pb210, listw = spatial_weights, zero.policy = TRUE, alternative = "two.sided")
SM_C <- moran.test(df_spat$Cs137, listw = spatial_weights, zero.policy = TRUE, alternative = "two.sided")
SM_A <- moran.test(df_spat$Am_241_h, listw = spatial_weights, zero.policy = TRUE, alternative = "two.sided")

# Calculate Z-scores for moran's I
Z_M_P <- SM_P$estimate[1] - SM_P$estimate[2]/sqrt(SM_P$estimate[3])
Z_M_C <- SM_C$estimate[1] - SM_C$estimate[2]/sqrt(SM_C$estimate[3])
Z_M_A <- SM_A$estimate[1] - SM_A$estimate[2]/sqrt(SM_A$estimate[3])

Morans.ix <- data.frame(Isotope = c("Lead-210", "Caesium-137", "Americum-241"), 
                        Moran = c(SM_P$estimate[1], SM_C$estimate[1], SM_A$estimate[1] ), 
                        Z = c(Z_M_P, Z_M_C, Z_M_A), 
                        P = c(SM_P$p.value, SM_C$p.value, SM_A$p.value))
colnames(Morans.ix)[2] <- "Moran's I statistic"

write.csv(Morans.ix, "Output/Spatial AC - results.csv")

# Compute lagged variable
df_spat$lag_Pb <- lag.listw(spatial_weights, df_spat$Pb210)
df_spat$lag_Cs <- lag.listw(spatial_weights, df_spat$Cs137)
df_spat$lag_Am <- lag.listw(spatial_weights, df_spat$Am_241_h)