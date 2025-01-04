library(data.table)
library(stringr)

# load gene markers object
dat = readRDS("clustering_coral/scdb/Spis_mc_coral_markers_list.RDS")

# make markers table
dt = data.table(gene=dat$genes)

# load gene expression data in metacells
mc = read.table("clustering_coral/scdb/Spis_coral_metacell_gene_FC")
mt = as.matrix(mc)
colnames(mt) <- str_remove(colnames(mt), "X") # a fix because read.table adds 'X' in front of column names that start with number

# get cell type with max fc for every gene
max_mc <- structure(
    colnames(mt)[apply(mt, 1, which.max)],
    names = rownames(mt)
)
stopifnot(all(dt$gene %in% names(max_ct)))

# add to table
dt[, metacell := max_mc[gene]]
dt[, metacell := as.integer(metacell)]
stopifnot(!any(is.na(dt$metacell)))

# get metacell to cell type assignment
ann = fread("clustering_coral/scdb/Spis_coral_metacell_annotation")

# add cell type to table
dt <- merge.data.table(dt, ann, by="metacell", all.x=TRUE, sort=FALSE)

# save
fwrite(dt, "clustering_coral/scdb/Spis_mc_coral_markers_list_cell_type.tsv", sep = "\t")