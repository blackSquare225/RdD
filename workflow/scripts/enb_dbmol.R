library(tidyverse)
library(readxl)

files <- list.files(path = "data", pattern = "_ENB.xlsx", full.names = TRUE)

list <- lapply(files, function(file) {

	print(file)
	sample <- gsub(".xlsx$", "", basename(file))

	var <- as.data.frame(read_excel(file))
	var <- var[,c(1,3,4,5,6,9,10,11)]
	var <- separate(var, Locus, into = c("Chromosome", "Position"), sep = ":", remove = FALSE)
	var <- var[,c("Chromosome","Position","Ref","Observed Allele","Type","Genes","Allele Frequency %","Amino Acid Change", "Coverage")]
	colnames(var) <- c("CHR","POS","REF","ALT","TYPE","SYMBOL","AF","AminoAcidChange", "Coverage")
	var$AF <- as.numeric(var$AF)
	var$AF <- trunc(var$AF / 100 * 100) / 100

	var$sample <- rep(sample, nrow(var))

	var$SYMBOL <- gsub("MFSD11,MIR636,", "", var$SYMBOL)
	var$SYMBOL <- gsub(",TET2-AS1", "", var$SYMBOL)
	var$SYMBOL <- gsub(",U2AF1L5", "", var$SYMBOL)

	return(var)

	})

df <- do.call(rbind, list)

df_wider <- df %>%
	select(sample, SYMBOL) %>%
	pivot_wider(names_from = SYMBOL, values_from = SYMBOL, 
		values_fn = list(SYMBOL = ~ as.numeric(length(.) > 0)),
		values_fill = list(SYMBOL = 0)) %>%
	as.data.frame()

rownames(df_wider) <- df_wider$sample

# Remove the 'sample' column
df_wider <- df_wider[, -which(names(df_wider) == "sample")]

# Convert to matrix
matrix <- as.matrix(df_wider)

ptsEnb <- readLines("resources/pts_enb.txt") 

matEnb <- matrix[rownames(matrix) %in% ptsEnb,]
dfEnb <- df %>%
	filter(sample %in% ptsEnb) 
