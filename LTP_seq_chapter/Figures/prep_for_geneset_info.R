########################################################################################################
# Calculating gene set lengths and overlaps for figures
########################################################################################################

library(tidyverse)
library(VennDiagram)
library(eulerr)

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

LTP_120 <- list(Total_120 = Total_120$entrezgene, Trap_120 = Trap_120$entrezgene)

LTP_120_overlap <- calculate.overlap(LTP_120)

names(LTP_120_overlap) <- c("Total_120", "TRAP_120", "Overlap")

LTP_60 <- list(Total_60 = Total_60$entrezgene, Trap_60 = Trap_60$entrezgene)

LTP_60_overlap <- calculate.overlap(LTP_60)

names(LTP_60_overlap) <- c("Total_60", "TRAP_60", "Overlap")

LTP_30 <- list(Total_30 = Total_30$entrezgene, Trap_30 = Trap_30$entrezgene)

LTP_30_overlap <- calculate.overlap(LTP_30)

names(LTP_30_overlap) <- c("Total_30", "TRAP_30", "Overlap")

# export files

LTP_120_df <- as.data.frame(sapply(LTP_120_overlap, '[', seq(max(sapply(LTP_120_overlap, length)))))

LTP_120_df <- mutate_all(LTP_120_df, as.character)

LTP_60_df <- as.data.frame(sapply(LTP_60_overlap, '[', seq(max(sapply(LTP_60_overlap, length)))))

LTP_60_df <- mutate_all(LTP_60_df, as.character)

LTP_30_df <- as.data.frame(sapply(LTP_30_overlap, '[', seq(max(sapply(LTP_30_overlap, length)))))

LTP_30_df <- mutate_all(LTP_30_df, as.character)

write_delim(LTP_120_df, "LTP_120.txt", delim = "\t", na = "")

write_delim(LTP_60_df, "LTP_60.txt", delim = "\t", na = "")

write_delim(LTP_30_df, "LTP_30.txt", delim = "\t", na = "")

# re format for BioVinci programme and write

write.table(LTP_120_venn, "LTP_120_venn.csv", sep = ",", row.names = F)

write.table(LTP_60_df, "LTP_60_venn.csv", sep = ",", row.names = F)

write.table(LTP_30_df, "LTP_30_venn.csv", sep = ",", row.names = F)

# play around with euler diagram

LTP_120_euler <- euler(LTP_120)

LTP_120_plot <- plot(LTP_120_euler, fills = c("#FE7B0D", "#FEC392"), legend = TRUE, quantities = TRUE)
 
print(LTP_120_plot)
 
LTP_60_euler <- euler(LTP_60)

LTP_60_plot <- plot(LTP_60_euler, fills = c("#7A0CB5", "#B580D2"), legend = TRUE, quantities = TRUE)

print(LTP_60_plot)

LTP_30_euler <- euler(LTP_30)

LTP_30_plot <- plot(LTP_30_euler, fills = c("#1C7E1C", "#5DBF5D"), legend = TRUE, quantities = TRUE)

print(LTP_30_plot) 
 
# total

Total <- list(Total_120 = Total_120$entrezgene, Total_60 = Total_60$entrezgene, Total_30 = Total_30$entrezgene)

Total_overlap <- calculate.overlap(Total)

Total_euler <- euler(Total)

plot(euler(Total, shape = "ellipse"), quantities = TRUE) 

names(Total_overlap) <- c("Total_120_Total_60_Total_30", "Total_120_Total_60", "Total_120_Total_30",
                          "Total_60_Total_30","Total_120", "Total_60", "Total_30") # this is manually calibrated

LTP_Total_df <- as.data.frame(sapply(Total_overlap, '[', seq(max(sapply(Total_overlap, length)))))

LTP_Total_df <- mutate_all(LTP_Total_df, as.character)

write_delim(LTP_Total_df, "LTP_total.txt", delim = "\t", na = "")

# Trap

Trap <- list(Trap_120 = Trap_120$entrezgene, Trap_60 = Trap_60$entrezgene, Trap_30 = Trap_30$entrezgene)
 
Trap_overlap <- calculate.overlap(Trap) 

Trap_euler <- euler(Trap) 

plot(euler(Trap, shape = "ellipse"), quantities = TRUE)   

names(Trap_overlap) <- c("Trap_120_Trap_60_Trap_30", "Trap_120_Trap_60", "Trap_120_Trap_30",
                         "Trap_60_Trap_30", "Trap_120", "Trap_60", "Trap_30") # manually calibrated

LTP_Trap_df <- as.data.frame(sapply(Trap_overlap, '[', seq(max(sapply(Trap_overlap, length)))))

LTP_Trap_df <- mutate_all(LTP_Trap_df, as.character)

write_delim(LTP_Trap_df, "LTP_trap.txt", delim = "\t", na = "")
