#######################################################################################
# Volcano plots for 2 and 5 hour
######################################################################################

library(tidyverse)

######################
# variables
######################

recall_path_2hr <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_2hr/Limma/results/rm_outliers/RUVg/prep/output_files/recall_k_21.txt"

extinction_path_2hr <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_2hr/Limma/results/rm_outliers/RUVg/prep/output_files/extinction_k_21.txt"

recall_path_5hr <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_5hr/results/RUVg/output_files/recall_k_29.txt"

extinction_path_5hr <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_5hr/results/RUVg/output_files/recall_k_29.txt"

recallExtinction_path_5hr <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_5hr/results/RUVg/output_files/recall_extinction_k_29.txt"
#####################
# read in files
####################

cols_want <- c("ensembl", "symbol", "entrez_gene", "logFC","adj.P.Val")

recall_2hr <- read_delim(recall_path_2hr, delim = "\t") |>
  select(cols_want) |>
  add_column(group = "Recall") |>
  mutate(deg = ifelse(adj.P.Val <= 0.05, "yes", "no"))

extinction_2hr <- read_delim(extinction_path_2hr, delim = "\t")|>
  select(cols_want) |>
  add_column(group = "Extinction")|>
  mutate(deg = ifelse(adj.P.Val <= 0.05, "yes", "no"))

recall_5hr <- read_delim(recall_path_5hr, delim = "\t")|>
  add_column(group = "Recall")

extinction_5hr <- read_delim(extinction_path_5hr, delim = "\t")|>
  add_column(group = "Extinction")

recallExtinction_5hr <- read_delim(recallExtinction_path_5hr, delim = "\t") |>
  add_column(group = "Recall_Extinction")

####################
# bind files
####################

hr_2 <- left_join(recall_2hr, extinction_2hr, by = c("ensembl", "symbol", "entrez_gene"))

hr_2 <- mutate(hr_2, deg_group = case_when(group.x == "Recall" & deg.x == "yes" &
                                                group.y == "Extinction" & deg.y == 
                                                "no" ~ "recall_only",
                  group.x == "Recall" & deg.x == "yes" & group.y == "Extinction" & 
                    deg.y == "yes" ~ "both",
                  group.x == "Recall" & deg.x == "no"& group.y == "Extinction" & 
                    deg.y == "yes" ~ "extinction_only",
                  group.x == "Recall" & deg.x == "no" &
                    group.y == "Extinction" & deg.y == "no" ~ "non_sig"))

recall_hr_2 <- select(hr_2, "ensembl", "symbol", "logFC.x", "adj.P.Val.x","deg_group") 
# need to replace extinction only with non_sig
