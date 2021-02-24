#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
tree_file <- args[1]
gene_expression_file <- args[2]
out_dir <- args[3] 
fc_out <- as.numeric(args[4])
abs_leakyness <- as.numeric(args[5])
abs_leakynessbg <- as.numeric(args[6])

source("/users/asebe/aelek/proj/scRNAseq_spis/Stylophora_single_cell_atlas/metacell_downstream_functions/Tree_functions.R")

#wdir <- "/users/asebe/aelek/proj/scRNAseq_spis/Stylophora_single_cell_atlas/cross_stages/tree"
#setwd(wdir)

#tree <- readRDS(file.path(wdir,"cell_type_tree_vargenes1.8_coocaveragep0.75hclusth0.75,0.95_politomies0.10.RDS"))
tree <- readRDS(tree_file)

#mc_fp <- read.table(file.path(wdir,"Spis_all_stages_gene_FC.tsv"))
mc_fp <- read.table(gene_expression_file)
colnames(mc_fp) <- str_replace(colnames(mc_fp),"\\.","-")

dir.create(out_dir,showWarnings=FALSE)

library(foreach)
library(doParallel)
ncores=64
registerDoParallel(ncores)
fcrange_in <- seq(from=1,to=2,by=0.1)
#fcrange_out <- seq(from=1,to=2,by=0.1)
datadir <- file.path(out_dir,sprintf("abs%s_abs%s",abs_leakyness,abs_leakynessbg))
dir.create(datadir,showWarnings=FALSE)
foreach (i=fcrange_in) %dopar% {
    exist_files <- list.files(datadir,pattern="*.RDS")
    out_file <- sprintf("gap_genes_spis_infc_%s_outfc_%s.RDS",i,fc_out)
    if (!out_file %in% exist_files) {
          print(sprintf("Will calculate %s/%s",datadir,out_file))
        gap_genes <- tree_gap_genes(
          tree = tree, feature_matrix = mc_fp,
          branch_length_thrs = 0.01,
          feature_in_thrs = i, feature_out_thrs = fc_out,
	  method="absolute", methodbg="absolute",
	  abs_leakyness=abs_leakyness, abs_leakynessbg=abs_leakynessbg,
          verbose = FALSE, ncores=1
        )
        saveRDS(gap_genes, file.path(datadir,out_file))
    }
}
stopImplicitCluster()
