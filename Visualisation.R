### Global ----------------------------------------------------------------
source("Data_analysis.R")

library(ggplot2)
library(ggcorrplot)
library(ggpubr)
library(visreg)

### Rename variables -----------------------------------------------------------
rownames(cm) <- c("Pb-210","Cs-137","Am-241^","Organic matter","Avg. altitude","Max. altitude", "Surface area", "Prec. GF", "Prec. CHB")
colnames(cm) <- rownames(cm)

levels(df_l$Isotope) <- c("Pb-210","Cs-137","Am-241")
levels(df_ms$Isotope) <- c("Pb-210","Cs-137","Am-241")

### Raw data plots --------------------------------------------------------
## Fig. 2a - Activity concentrations variation between glaciers
df_box <- df_l 
df_box$Glacier <- factor(df_box$Glacier, levels = df_mod[order(df_mod$Lon), ]$Glacier)

g1 <- ggplot(df_box, aes(Glacier, value, color = Isotope, fill = Isotope))

g1.o <- g1 + geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter() + facet_wrap(.~Isotope, ncol = 1, nrow = 3, scales = "free") + 
  geom_hline(data = df_ms, aes(yintercept = value)) + 
  theme_classic() + ylab(expression(Activity~concentration~(Bq~kg^{-1}))) + 
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none");g1.o

## Fig. 2b - Variables correlation
g2 <- ggcorrplot(cm, hc.order = FALSE, type = "lower", lab = TRUE, 
                 tl.cex = 9, lab_size = 3,
                 outline.col = "white",
                 colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation", ggtheme = ggplot2::theme_classic) +
  theme(axis.text.x = element_text(margin = margin(-2,0,0,0), angle = 30, vjust = 0.9, hjust = 1),
        axis.text.y = element_text(margin = margin(0,-2,0,0), angle = 45, vjust = -1, hjust = 1),
        panel.grid.minor = element_line(size = 10), legend.key.size = unit(0.4, 'cm'), 
        legend.text = element_text(size = 7, vjust = 0.8), legend.title = element_text(angle = 0, size = 8))

## Fig. 2c - table with explained variance 
gt <- ggtexttable(R_r2, rows = NULL, cols = c("Linear regression", "adj. R-squared"), 
                  theme = ttheme(base_size = 8, base_style = "classic"))

## Fig. 2d - correlation of isotopes 
gx <- ggplot(data = df_mod, aes(x = log(Cs137), y =  log(Pb210)))
gx.o <- gx + geom_point(aes(color = Am241), size = 2.5) + geom_smooth(method = "lm", color = "grey", alpha = 0.2) + 
  theme_classic() + theme(legend.key.size = unit(0.4, 'cm'), legend.title = element_text(angle = 0, size = 7), legend.text = element_text(size = 8)) +
  scale_color_gradient(name = expression(paste(""^"241", "Am (Bq kg" ^-1,")*")), low =  "#ffbc51",
                       high = "black") + 
  xlab(expression(paste(""^"137", "Cs (Bq kg" ^-1,")*"))) + ylab(expression(paste(""^"210", "Pb (Bq kg" ^-1,")*")));gx.o

## Arrange and save Fig. 2  
a1 <- ggarrange(g1.o, 
                ggarrange(g2, gt, gx.o, ncol = 1, heights = c(1.5,0.5,1), labels = c("B", "C", "D")), 
                ncol = 2, widths = c(1.1, 1), labels = c("A")); a1 

ggsave("Output/Fig. 2 - Panneled.png", scale = 1, width = 2800, height = 2000, units = "px", bg = "white")

### Spatial autocorrelation 
## Lead-210
S.g1 <- ggplot(data.frame(df_spat), aes(Pb210, lag_Pb)) + 
  geom_point(color = "#a1a391", size = 3.5, shape = 16) + 
  geom_smooth(method = "lm", color = "#91b1a0", se = FALSE, size = 2) + 
  theme() + theme_classic(base_size = 18) + 
  xlab(expression(paste(""^"210", "Pb [Bq kg" ^-1,"]"))) + 
  ylab(expression(paste("Spatially lagged "^"210", "Pb [Bq kg" ^-1,"]"))) + 
  geom_hline(yintercept = mean(df_spat$lag_Pb), linetype = "dashed", 
             color = "grey50", linewidth = 0.5) + 
  geom_vline(xintercept = mean(df_spat$Pb210), linetype = "dashed", 
             color = "grey50", linewidth = 0.5); S.g1
ggsave("Output/Paper_3a_Pb_lag.png", height = 12, width = 16, dpi = 600, units = "cm")

summary(lm(lag_Pb ~ Pb210, data = df_spat))

S.g2 <- ggplot(data.frame(df_spat), aes(Cs137, lag_Cs)) + 
  geom_point(color = "#a1a391", size = 3.5, shape = 17) + 
  geom_smooth(method = "lm", color = "#91b1a0", se = FALSE, size = 2) + 
  theme() + theme_classic(base_size = 18) + 
  xlab(expression(paste(""^"137", "Cs [Bq kg" ^-1,"]"))) + 
  ylab(expression(paste("Spatially lagged "^"137", "Cs [Bq kg" ^-1,"]"))) + 
  geom_hline(yintercept = mean(df_spat$lag_Cs), linetype = "dashed", 
             color = "grey50", linewidth = 0.5) + 
  geom_vline(xintercept = mean(df_spat$Cs137), linetype = "dashed", 
             color = "grey50", linewidth = 0.5); S.g2
ggsave("Output/Paper_3b_Cs_lag.png", height = 12, width = 16, dpi = 600, units = "cm")

summary(lm(lag_Cs ~ Cs137, data = df_spat))

S.g3 <- ggplot(data.frame(df_spat), aes(Am_241_h, lag_Am)) + 
  geom_point(color = "#a1a391", size = 3.5, shape = 18) + 
  geom_smooth(method = "lm", color = "#91b1a0", se = FALSE, size = 2) + 
  theme() + theme_classic(base_size = 18) + 
  xlab(expression(paste(""^"241", "Am [Bq kg" ^-1,"]"))) + 
  ylab(expression(paste("Spatially lagged "^"241", "Am [Bq kg" ^-1,"]"))) + 
  geom_hline(yintercept = mean(df_spat$lag_Am), linetype = "dashed", 
             color = "grey50", linewidth = 0.5) + 
  geom_vline(xintercept = mean(df_spat$Am_241_h), linetype = "dashed", 
             color = "grey50", linewidth = 0.5); S.g3
ggsave("Output/Paper_3c_Am_lag.png", height = 12, width = 16, dpi = 600, units = "cm")

summary(lm(lag_Am ~ Am_241_h, data = df_spat))

### Models results - partial residuals-------------------------------------
## Models_coefficients
mP.4_coe <- summary(mP.4)$coefficients
mC.2_coe <- summary(mC.2)$coefficients
mA.1_coe <- data.frame(summary(mA.1)$estimate)

## Fig. 3_P2 - Models results (Lead)
g3.1 <- visreg(mP.4, "OM", gg = TRUE, line.par = list(col = "#e78f06"))
g3.1 + theme_classic(base_size = 18) + geom_point(size = 3.5, col = "#e78f06", alpha = 0.5) +
  theme() +
  xlab("Average organic matter content \n in cryoconite per glacier [%]") + ylab(expression(paste("log("^"210", "Pb [Bq kg" ^-1,"])")))
ggsave("Output/Paper_4a_Pb_OM.png", height = 12, width = 16, dpi = 600, units = "cm")

g3.2 <- visreg(mP.4, "Average_altitude", gg = TRUE, line.par = list(col = "#2c5b95"))
g3.2 + theme_classic(base_size = 18) + geom_point(size = 3.5, col = "#2c5b95", alpha = 0.5) +
  theme() +
  xlab("Average altitude of glaciers [m a.s.l.]") + ylab(expression(paste("log("^"210", "Pb [Bq kg" ^-1,"])")))
ggsave("Output/Paper_4b_Pb_altitude.png", height = 12, width = 16, dpi = 600, units = "cm")

g3.3 <- visreg(mP.4, "Size_kmq_l", gg = TRUE, line.par = list(col = "#3F002C"))
g3.3 + theme_classic(base_size = 18) + geom_point(size = 3.5, col = "#3F002C", alpha = 0.5) +
  theme() +
  xlab(expression(paste("log(surface area of glaciers", " [km" ^2,"])"))) + ylab(expression(paste("log("^"210", "Pb [Bq kg" ^-1,"])")))
ggsave("Output/Paper_4c_Pb_size.png", height = 12, width = 16, dpi = 600, units = "cm")

## Fig. 4_C2 - Models results (Ceasium) 
g4.1 <- visreg(mC.2, "OM", gg = TRUE, line.par = list(col = "#e78f06"))
g4.1 + theme_classic(base_size = 18) + geom_point(size = 3.5, col = "#e78f06", shape = 17, alpha = 0.5) +
  theme() +
  xlab("Average organic matter content \n in cryoconite per glacier [%]") + ylab(expression(paste("log("^"137", "Cs [Bq kg" ^-1,"])")))
ggsave("Output/Paper_5a_Cs_OM.png", height = 12, width = 16, dpi = 600, units = "cm")

g4.2 <- visreg(mC.2, "Size_kmq_l", gg = TRUE, line.par = list(col = "#3F002C"))
g4.2 + theme_classic(base_size = 18) + geom_point(size = 3.5, col = "#3F002C", shape = 17, alpha = 0.5) +
  theme() +
  xlab(expression(paste("log(surface area of glaciers", " [km" ^2,"])"))) + ylab(expression(paste("log("^"137", "Cs [Bq kg" ^-1,"])")))
ggsave("Output/Paper_5b_Cs_size.png", height = 12, width = 16, dpi = 600, units = "cm")