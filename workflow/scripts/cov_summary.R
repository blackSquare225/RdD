library(data.table)
library(tidyverse)
library(writexl)

rm(list = ls())

files <- list.files(path = "alignments", pattern=".mosdepth.regions.bed.gz$", full.names=TRUE)
names <- gsub(".mosdepth.regions.bed.gz$","",basename(files))

cov <- lapply(files, function(file) {
	df <- fread(file)
	df$sample <- gsub(".mosdepth.regions.bed.gz$","",basename(file))
	return(df)
	})

df <- do.call(rbind, cov)

# Perform Summary of coverage from all runs and samples 
summary(df$V4)

# Perform coverage per run
df_ <- df %>%
	mutate(RUN = gsub("_barcode.*","",sample)) %>%
	group_by(RUN) %>%
  summarise(
    Mean = mean(V4, na.rm = TRUE),
    Median = median(V4, na.rm = TRUE),
    Min = min(V4, na.rm = TRUE),
    Max = max(V4, na.rm = TRUE),
    Q2 = quantile(V4, 0.25, na.rm = TRUE), # Second quartile (25th percentile)
    Q3 = quantile(V4, 0.75, na.rm = TRUE)  # Third quartile (75th percentile)
  ) %>%
	ungroup() %>%
	as.data.frame() 


summary(df$V4)



write_xlsx(x=df_, path="docs/CovForRuns.xlsx")
