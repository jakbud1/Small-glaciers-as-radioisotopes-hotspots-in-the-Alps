### global ----------------------------------------------------------------
library(readxl)

### data input
file_list <- list.files("Input/Precipitation/1962", full.names = TRUE)

### dfs for each year
## 62
df_62 <- read_xlsx("Input/Precipitation/1962/1962_1.xlsx")

for (i in 2:length(file_list)){
  df_62 <- cbind(df_62, read_xlsx(file_list[i])[4])
}; df_62 
df_62$sum_62 <- apply(as.matrix(df_62[,-c(1:3)]), 1, FUN = sum)

## 63
df_63 <- read_xlsx("Input/Precipitation/1963/1963_1.xlsx")

for (i in 2:length(file_list)){
  df_63 <- cbind(df_63, read_xlsx(file_list[i])[4])
}; df_63
df_63$sum_63 <- apply(as.matrix(df_63[,-c(1:3)]), 1, FUN = sum)

## 64
df_64 <- read_xlsx("Input/Precipitation/1964/1964_1.xlsx")

for (i in 2:length(file_list)){
  df_64 <- cbind(df_64, read_xlsx(file_list[i])[4])
}; df_64
df_64$sum_64 <- apply(as.matrix(df_64[,-c(1:3)]), 1, FUN = sum)

### Overall sum -----------------------------------------------------------
df_prec <- cbind(df_62[,c(1:3, 16)], df_63[,16], df_64[,16])
df_prec$prec_GLF <- rowSums(df_prec[,c(4:6)])

### Add Chernobyl (May 1986) precipitation --------------------------------
df_CHB <- read_xlsx("Input/Precipitation/1986_05.xlsx")
df_prec$prec_CHB <- df_CHB$prec_may86
rm(file_list, i, df_62, df_63, df_64)

colnames(df_prec)[c(5,6)] <- c("sum_63", "sum_64")  

### Export data -----------------------------------------------------------
write.csv(df_prec, "Output/precipitation.csv")
writexl::write_xlsx(df_prec, "Output/precipitation.xlsx")

rm(df_CHB)