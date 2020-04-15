########################################################################################################
# Calculating gene set lengths and overlaps for figures
########################################################################################################

library(tidyverse)
library(VennDiagram)

# set wd

setwd("~/PhD/Thesis/Results/LTP_seq_chapter/Figure_prep")

# import files

colNames <- c("entrezgene", "geneset")

Total_120 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_120/Ready_for_MAGMA/LTP_total_001.txt", 
                        delim = " ", col_names = colNames)

Trap_120 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_120/Ready_for_MAGMA/LTP_TRAP_001.txt", 
                        delim = " ", col_names = colNames)

Total_60 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_60/Ready_for_MAGMA/LTP_total_001_60.txt", 
                        delim = " ", col_names = colNames)

Trap_60 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_60/Ready_for_MAGMA/LTP_TRAP_001_60.txt", 
                       delim = " ", col_names = colNames)

Total_30 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_30/Ready_for_MAGMA/LTP_total_001_30.txt", 
                       delim = " ", col_names = colNames)

Trap_30 <- read_delim("../../../../LTP_seq/prep_for_disease_association/output/LTP_30/Ready_for_MAGMA/LTP_TRAP_001_30.txt", 
                       delim = " ", col_names = colNames)

# generate lists and calculate overlap

LTP_120 <- list(Total_120$entrezgene,Trap_120$entrezgene)

LTP_120 <- calculate.overlap(LTP_120)

names(LTP_120) <- c("Total_120", "TRAP_120", "Overlap")

LTP_60 <- list(Total_60$entrezgene, Trap_60$entrezgene)

LTP_60 <- calculate.overlap(LTP_60)

names(LTP_60) <- c("Total_60", "TRAP_60", "Overlap")

LTP_30 <- list(Total_30$entrezgene, Trap_30$entrezgene)

LTP_30 <- calculate.overlap(LTP_30)

names(LTP_30) <- c("Total_30", "TRAP_30", "Overlap")

# export files

LTP_120_df <- as.data.frame(sapply(LTP_120, '[', seq(max(sapply(LTP_120, length)))))

LTP_60_df <- as.data.frame(sapply(LTP_60, '[', seq(max(sapply(LTP_60, length)))))

LTP_30_df <- as.data.frame(sapply(LTP_30, '[', seq(max(sapply(LTP_30, length)))))

write_delim(LTP_120_df, "LTP_120.csv", delim = ",", na = "")

write_delim(LTP_60_df, "LTP_60.csv", delim = ",", na = "")

write_delim(LTP_30_df, "LTP_30.csv", delim = ",", na = "")
