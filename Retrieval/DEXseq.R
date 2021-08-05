##############################################################################################################
# Pre processing for DEXSeq from featureCount files
##############################################################################################################
library(DEXSeq)
library(edgeR)
library(tidyverse)

# set wd and variables

setwd("~/Documents/PhD/Retrieval_2hr/exon_featurecount_files/no_dups/")

output_path <- "~/Documents/PhD/Retrieval_2hr/DEXSeq_results/"

# get samples

samples <- as.character(1:36) # number of samples

outliers <- as.character(c("6", "14", "21"))#as.character(c("6", "14","21"))

samples <- samples[!samples %in% outliers]

files <- c(paste0("CA1_CFC_", samples, ".markdup.featurecount"))

# read in as DGE list, convert to matrix for DEXSeq

x <- readDGE(files, columns = c(1,7), comment.char = "#")

sample_names <- str_extract(x$samples$files, "CA1_CFC_[1-9][0-9]?")

sample_names_short <- str_extract(sample_names, "[1-9][0-9]?$")

colnames(x) <- sample_names_short

sample_info <- read_delim("../../sample_condition_summary.txt", delim = "\t") %>%
  filter(!sampleID %in% outliers) 

row.names(sample_info) <- sample_info$sampleID

x_matrix <- as.matrix(x$counts)

# biomaRt
# Get gene names that link to exon ID

exonid <- as.data.frame(rownames(x_matrix))

mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl") # host arg. host = "asia.ensembl.org"

genes_biomart <- getBM(attributes = c("ensembl_exon_id", "ensembl_gene_id","external_gene_name"),
                       
                       filters = "ensembl_exon_id",
                       
                       values = exonid$`rownames(x_matrix)`,
                       
                       mart = mart)

genes <- left_join(exonid, genes_biomart, by = c("rownames(x_matrix)" = "ensembl_exon_id"))

colnames(genes) <- c("ensembl_exon","ensembl_gene", "symbol")

# create DEXSeq object

dxd <- DEXSeqDataSet(x_matrix, sample_info, design = ~sample + exon + condition:exon,
                     featureID = genes$ensembl_exon, groupID = genes$ensembl_gene)

# run DEXSeq

dxd <- estimateSizeFactors(dxd)

dxd <- estimateDispersions(dxd)
