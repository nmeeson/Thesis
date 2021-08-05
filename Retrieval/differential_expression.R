# CFC retrieval differential expression analysis 2 hour

# load packages

library(tidyverse)
library(biomaRt)
library(limma)
library(edgeR)
library(ggExtra)
library(gtools)
library(RColorBrewer)
library(ggrepel)

# get date for log file

start_time <- date() 

date <- str_split(start_time, pattern = " ", simplify = TRUE)

date <- paste(date[3], date[2], date[5], sep = "_")

# set working directory

setwd("~/Documents/PhD/Retrieval_2hr/Feature_count_files/featurecounts/")

# set variables

n_group <- 11 # how many samples per group

output_path <- "~/Documents/PhD/Retrieval_2hr/DGE_results/"

design_formula <- "~ 0 + group + RIN + pool"
p_cutoff <- 0.1

log_file_name <- str_glue("logs/", "log_{date}", ".txt")
log_description <- "Differential gene expression 2 hours excluding outliers, limma voom"


# read in files, firstly make a vector with all the file names in, then read in with readDGE from edgeR

samples <- as.character(1:36) # number of samples

outliers <- as.character(c("6", "14", "21"))#as.character(c("6", "14","21"))

samples <- samples[!samples %in% outliers]

files <- c(paste0("CA1_CFC_", samples, ".markdup.featurecount")) # file names

x <- readDGE(files, columns = c(1,7), comment.char = "#")

n_samples <- length(samples)

## set WD for rest of script, create directories

setwd("~/Documents/PhD/Retrieval_2hr/DGE_results/")


# initiate log file

cat(str_glue("\n","Log file for : {log_description}", "Design: {design_formula}", .sep = "\n"), file = log_file_name, append = TRUE)

#########################################################################################################
# set up sample information
#########################################################################################################

# sample names

sample_names <- str_extract(x$samples$files, "CA1_CFC_[1-9][0-9]?")

sample_names_short <- str_extract(sample_names, "[1-9][0-9]?$")

colnames(x) <- sample_names_short

# read in a file that has the info you want

sample_info <- read_delim("../sample_condition_summary.txt", delim = "\t") %>%
  filter(!sampleID %in% outliers) 


# change col types
sample_info$rat <- as.factor(sample_info$sampleID)
sample_info$condition <- as.factor(sample_info$condition)
sample_info$pool <- as.factor(sample_info$pool)

#sample_info$isolation_batch <- as.factor(sample_info$isolation_batch)
#sample_info$library_batch <- as.factor(sample_info$library_batch)

# add variables to x

x$samples$rat <- sample_info$rat
x$samples$group <- sample_info$condition
x$samples$pool <- sample_info$pool
x$samples$RIN <- sample_info$RIN
x$samples$freeze <- sample_info$percent_freezing
#x$samples$isolation_batch <- sample_info$isolation_batch
#x$samples$library_batch <- sample_info$library_batch
#x$samples$library_bp <- sample_info$library_bp
#x$samples$library_conc <- sample_info$library_conc
x$samples$library_nM <- sample_info$library_nM


########################################################################################################
# get gene annotations using biomaRt, other host arguments are asia.ensembl.org and useast.ensembl.org
########################################################################################################
geneid <- as.data.frame(rownames(x))

mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl") # host arg. host = "asia.ensembl.org"

genes_biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "wikigene_description", "gene_biotype" ),
                       
                       filters = "ensembl_gene_id",
                       
                       values = geneid$`rownames(x)`,
                       
                       mart = mart)

genes <- left_join(geneid, genes_biomart, by = c("rownames(x)" = "ensembl_gene_id"))

colnames(genes) <- c("ensembl", "symbol", "entrez_gene", "gene_description", "gene_type")

# remove duplicate IDs 

genes <- genes[!duplicated(genes$ensembl),]

# add to sample info

x$genes <- genes

names(x$genes) <- c("ensembl", "symbol", "entrez_gene", "gene_description", "gene_type")

# filter for protein coding genes

#x <- x[x$genes$gene_type == "protein_coding", ]

#########################################################################################################
# data preprocessing 
#########################################################################################################
# generating CPM and log CPM (counts per million and log counts per million), can also use RPKM or FPKM if
# have gene lengths

cpm <- cpm(x)

lcpm <- cpm(x, log = TRUE)

# generate mean (L) and median (M) library size and show in e notation 

L <- mean(x$samples$lib.size) * 1e-6 
M <- median(x$samples$lib.size) * 1e-6 
summary (lcpm)

# add to log file

cat(str_glue("\n", "Mean library size : {L}", "Median library size: {M}", .sep = "\n"), file = log_file_name, append = TRUE)

## remove genes that are lowly expressed 

no_expression <- table(rowSums(x$counts ==0)==length(samples)) 

cat(str_glue("\n", "Genes with no expression: {no_expression[2]}", "Genes remaining: {no_expression[1]}", .sep = "\n"), file = log_file_name, append = TRUE)

keep <- filterByExpr(x, group = x$samples$condition)

x_all <- x

x <- x[keep, , keep.lib.sizes = FALSE]

cat(str_glue("\n", "Average library size: {mean(x$samples$lib.size)}", 
             "Min library size: {min(x$samples$lib.size)}",
             "Max library size: {max(x$samples$lib.size)}", "\n",
             .sep = "\n"), file = log_file_name, append = TRUE)

# calculate normalisation factors 

x <- calcNormFactors(x, method = "TMM")

cat(str_glue("\n", "Min norm factor: {min(x$samples$norm.factors)}",
             "Max norm factor: {max(x$samples$norm.factors)}", 
             .sep = "\n"), file = log_file_name, append = TRUE)

# create model matrix

design <- model.matrix(as.formula(design_formula), data = x$samples)

contr.matrix <- makeContrasts(
  norecallvsRecall = grouprecall - groupnoRecall,
  norecallvsExtinction = groupextinction - groupnoRecall,
  recallvsExtinction = groupextinction - grouprecall,
  levels = colnames(design))


# limma trend

# logCPM <- cpm(x, log = TRUE, prior.count = 3)
# tfit <- lmFit(logCPM, design)
# tfit <- (contrasts.fit(tfit, contrasts = contr.matrix))
# efit <- eBayes(tfit, trend = TRUE)


# fit the model voom with quality weights

# v <- voomWithQualityWeights(x, design, plot = FALSE)
# 
# vfit <- lmFit(v, design)
# 
# vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
# 
# efit <- eBayes(tfit)

# voom

v <- voom(x, design)

vfit <- lmFit(v, design)

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

efit <- eBayes(vfit)


# Save files using topTable

recall <-topTable(efit, coef = 1, n= Inf) %>%
  arrange(adj.P.Val)

extinction <- topTable(efit, coef = 2, n = Inf) %>%
  arrange(adj.P.Val)

recall_extinction <- topTable(efit, coef = 3, n = Inf) %>%
  arrange(adj.P.Val)

# write_delim(recall, file = str_glue(output_path, "files/",  
#                                     "recall.txt"), delim = "\t")
# 
# write_delim(extinction, file = str_glue(output_path, "files/", 
#                                         "extinction.txt"), delim = "\t")
# 
# write_delim(recall_extinction, file = str_glue(output_path, "files/", 
#                                                "recall_extinction.txt"), delim = "\t")
# save differentially expressed genes

recall_deg <- filter(recall, adj.P.Val <= p_cutoff)

extinction_deg <- filter(extinction, adj.P.Val <= p_cutoff)

recall_extinction_deg <- filter(recall_extinction, adj.P.Val <= p_cutoff)

# write_delim(recall_deg, file = str_glue(output_path, "files/",  
#                                     "recall_deg.txt"), delim = "\t")
# 
# write_delim(extinction_deg, file = str_glue(output_path, "files/", 
#                                         "extinction_deg.txt"), delim = "\t")
# 
# write_delim(recall_extinction_deg, file = str_glue(output_path, "files/", 
#                                                "recall_extinction_deg.txt"), delim = "\t")

# explore differentially expressed genes

de_summary <- summary(decideTests(efit, p.value = p_cutoff))

cat(str_glue("\n", "Differential expression for design {design_formula}, with p-value cut off of {p_cutoff}", "\n",
             .sep = "\n"), file = log_file_name, append = TRUE)

invisible(capture.output(de_summary, file = log_file_name, append = TRUE))

# plot venn diagram

de.genes <- decideTests(efit, p.value = p_cutoff)

# tiff(str_glue(output_path,"figures/",
#               "Venn_diagram__{p_cutoff}.tiff"),
#      width = 21, height = 21, units = "cm", res = 300)
# 
# vennDiagram(de.genes[,1:2], circle.col= c("#D95F02", "#7570B3"), 
#             names = c("Recall", "Extinction"))
# 
# mtext(str_glue("FDR = {p_cutoff}"),cex = 1.5, line = 1 )
# dev.off()

# volcano plots

cols_want <- c("ensembl", "symbol", "entrez_gene", "logFC","adj.P.Val")

recall_vol <- select(recall, cols_want) |>
  add_column(group = "Recall") |>
  mutate(deg = ifelse(adj.P.Val <= p_cutoff, "yes", "no"),
         shared = ifelse(ensembl %in% extinction_deg$ensembl & adj.P.Val <= p_cutoff , "shared", "not_shared"),
         log_10_p = -log10(adj.P.Val),
         status = ifelse(deg == "yes" & shared =="shared", "DEG_shared",
                         ifelse(deg == "yes" & shared == "not_shared", "DEG_not_shared", "not_DEG")))

extinction_vol <- select(extinction, cols_want) |>
  add_column(group = "Extinction")|>
  mutate(deg = ifelse(adj.P.Val <= p_cutoff, "yes", "no"),
         shared = ifelse(ensembl %in% recall_deg$ensembl & adj.P.Val <= p_cutoff , "shared", "not_shared"),
         log_10_p = -log10(adj.P.Val),
         status = ifelse(deg == "yes" & shared =="shared", "DEG_shared",
                         ifelse(deg == "yes" & shared == "not_shared", "DEG_not_shared", "not_DEG")))


extinction_not_shared <- filter(extinction_vol, shared =="not_shared") %>%
  arrange(adj.P.Val)

# save dfs

recall_only <- filter(recall_vol, status == "DEG_not_shared") |>
  write_delim("group_rin_pool/voom/files/recall_only.txt", delim = "\t")

shared <- filter(recall_vol, status == "DEG_shared") |>
  write_delim("group_rin_pool/voom/files/shared_recall.txt", delim = "\t")

extinction_only <- filter(extinction_vol, status == "DEG_not_shared") |>
  write_delim("group_rin_pool/voom/files/extinction_only.txt", delim = "\t")

shared_e <- filter(extinction_vol, status == "DEG_shared") |>
  write_delim("group_rin_pool/voom/files/shared_extinction.txt", delim = "\t")

# recall plot

recall_plot <- ggplot(recall_vol, aes(x = logFC, y = log_10_p, col = status)) +
  geom_point(show.legend =  FALSE, size = 6)+
  geom_text_repel(aes(label =ifelse(ensembl %in% head(ensembl, 5) | (shared == "not_shared" & adj.P.Val <= p_cutoff) , as.character(symbol), "")), 
                  max.overlaps = 18, colour = ifelse(recall_vol$shared == "shared","#1B9E77","#D95F02"),
                  size = 8)+
  scale_color_manual(values = c("#D95F02","#1B9E77", "black")) +
  xlab("Log fold change") +
  ylab("-log10(adjusted FDR)") +
  ggtitle("Volcano plot: Recall vs Control") +
  scale_x_continuous(expand = c(0,0), limits = c(-1.5, 2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 22),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), size = 22),
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.text = element_text(size = 18, face = "bold"))+
  removeGrid()

ggsave(filename = str_glue(output_path,"figures/", "recall_volcano_plot.tiff"), 
       plot = recall_plot, width = 30, height = 25, units = "cm", dpi = 300)

# Extinction plot

extinction_plot <- ggplot(extinction_vol, aes(x = logFC, y = log_10_p, col = status)) +
  geom_point(show.legend =  FALSE, size = 6)+
  geom_text_repel(aes(label =ifelse(ensembl %in% head(ensembl, 5) | (ensembl %in% head(extinction_not_shared$ensembl, 5) & adj.P.Val <= p_cutoff) , as.character(symbol), "")), 
                  max.overlaps = 30, colour = ifelse(extinction_vol$shared == "shared","#1B9E77","#7570B3"),
                  size = 8)+
  scale_color_manual(values = c("#7570B3","#1B9E77", "black")) +
  xlab("Log fold change") +
  ylab("-log10(adjusted FDR)") +
  ggtitle("Volcano plot: Extinction vs Control") +
  scale_x_continuous(expand = c(0,0), limits = c(-1.5, 2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 22),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), size = 22),
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.text = element_text(size = 18, face = "bold"))+
  removeGrid()

ggsave(filename = str_glue(output_path,"figures/", "extinction_volcano_plot.tiff"), 
       plot = extinction_plot, width = 30, height = 25, units = "cm", dpi = 300)

