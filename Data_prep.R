### Global ----------------------------------------------------------------
source("Climate_data.R")

library(data.table)
library(dplyr)
library(scales)

### Wide format - all data 
df_all <- read_xlsx("Input/spec_cryo.xlsx")

### Long format - isotopes inter- and intra variation data
df_l <- melt(as.data.table(df_all[,c(2,3,5:7)]), id.vars = c("Sample_ID", "Glacier"), variable.name = "Isotope")
str(df_l)

### Overall mean for isotopes
df_ms <- aggregate(value ~ Isotope, data = df_l[,c(3:4)], FUN = mean)

### Mean for glaciers - all data
df_mod <- df_all[,c(3,5:9)] %>% 
  group_by(Glacier) %>% 
  summarise_all("mean", na.rm = TRUE)

### Mass of samples
df_mass <- df_all[,c(3,4)] %>% 
  group_by(Glacier) %>% 
  summarise_all("mean", na.rm = TRUE)

### Add meta and precipitation for df_mod
df_meta <- read_xlsx("Input/meta_glaciers.xlsx")
df_meta <- merge(df_meta, df_prec[, c(1, 4:8)], by = "ID"); rm(df_prec)
df_mod <- merge(df_mod, df_meta, by = "Glacier"); rm(df_meta)

### Standard deviations for each glacier
df_sd_out <- df_all[,c(3,5:9)] %>% 
  group_by(Glacier) %>% 
  summarise_all("sd", na.rm = TRUE)

write.csv(df_mod, "Output/df_mod.csv")
write.csv(df_sd_out, "Output/df_sd_out.csv")
