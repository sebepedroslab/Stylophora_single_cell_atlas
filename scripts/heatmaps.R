library(data.table)
library(stringr)

# get metacell to cell type assignment
ann = fread("clustering_coral/scdb/Spis_coral_metacell_annotation")
ann[, metacell := as.character(metacell)]

# get heatmap data
fn <- "Spis_coral_gene_expression_sc.RDS"
hmp = readRDS(file.path("clustering_coral/scdb/figs/", fn))

# get matrix
mt = hmp@matrix

# maka data table
dt = data.table(gene = rownames(mt))

# get metacell with max fc for every gene
max_mc = structure(
    colnames(mt)[apply(mt, 1, which.max)],
    names = rownames(mt)
)
stopifnot(all(dt$gene %in% names(max_mc)))
dt[, metacell := max_mc[gene]]

# add cell type to table
dt <- merge.data.table(dt, ann, by="metacell", all.x=TRUE, sort=FALSE)

# save
fwrite(dt, file.path("clustering_coral/scdb/figs/", str_replace(fn, ".RDS", ".tsv")), sep="\t")

