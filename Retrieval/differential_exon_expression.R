# CFC retrieval differential exon analysis 2 hour

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

setwd("~/Documents/PhD/Retrieval_2hr/exon_featurecount_files/")

# set variables

n_group <- 11 # how many samples per group

output_path <- "~/Documents/PhD/Retrieval_2hr/DEE_results/"

design_formula <- "~ 0 + group + RIN + pool"
p_cutoff <- 0.1

log_file_name <- str_glue("logs/", "log_{date}", ".txt")
log_description <- "Differential exon expression 2 hours excluding outliers, limma voom"


# read in files, firstly make a vector with all the file names in, then read in with readDGE from edgeR

samples <- as.character(1:36) # number of samples

outliers <- as.character(c("6", "14", "21"))#as.character(c("6", "14","21"))

samples <- samples[!samples %in% outliers]

files <- c(paste0("CA1_CFC_", samples, ".markdup.featurecount")) # file names

# remove duplicate exons before readDGE

for (f in files){
  file_f <- read_delim(f, delim = "\t", comment = "#") %>%
    filter(!duplicated(Geneid))
  write_delim(file_f, str_glue("no_dups/{f}"), delim = "\t")
}

setwd("~/Documents/PhD/Retrieval_2hr/exon_featurecount_files/no_dups/")

# read in files with readDGE

x <- readDGE(files, columns = c(1,7), comment.char = "#")

n_samples <- length(samples)

## set WD for rest of script, create directories

setwd("~/Documents/PhD/Retrieval_2hr/DEE_results/")


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

# Get gene names that link to exon ID

exonid <- as.data.frame(rownames(x))

mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl") # host arg. host = "asia.ensembl.org"

genes_biomart <- getBM(attributes = c("ensembl_exon_id", "ensembl_gene_id","external_gene_name"),
                       
                       filters = "ensembl_exon_id",
                       
                       values = exonid$`rownames(x)`,
                       
                       mart = mart)

genes <- left_join(exonid, genes_biomart, by = c("rownames(x)" = "ensembl_exon_id"))

colnames(genes) <- c("ensembl_exon","ensembl_gene", "symbol")

# add to sample info

x$genes <- genes

names(x$genes) <- c("ensembl_exon","ensembl_gene", "symbol")

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

# voom

v <- voom(x, design)

vfit <- lmFit(v, design)

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

efit <- eBayes(vfit) #not sure whether I need this

ex <- diffSplice(efit, geneid = "ensembl_gene")

invisible(capture.output(diffSplice(efit, geneid = "ensembl_gene"), file = log_file_name, append = TRUE))

# differential expression
recall <-topSplice(ex, coef = 1, test = "simes", n= Inf) %>%
  arrange(FDR)

recall_t <- topSplice(ex, coef = 1, test = "t", n = Inf) %>%
  arrange(FDR)

extinction <-topSplice(ex, coef = 2, test = "simes", n= Inf) %>%
  arrange(FDR)

recall_extinction <-topSplice(ex, coef = 3, test = "simes", n= Inf) %>%
  arrange(FDR)

write_delim(recall, file = str_glue(output_path, "files/",
                                    "recall_exon.txt"), delim = "\t")

write_delim(extinction, file = str_glue(output_path, "files/",
                                        "extinction_exon.txt"), delim = "\t")

write_delim(recall_extinction, file = str_glue(output_path, "files/",
                                               "recall_extinction_exon.txt"), delim = "\t")

# plotSplice

plotSplice(ex, coef = 1, genecolname = "symbol")

