library(data.table)
library(stringr)

# get metacell to cell type assignment
st = "coral" # coral, polyp or larva
ann = fread(sprintf("clustering_%s/scdb/Spis_%s_metacell_annotation", st, st))
ann[, metacell := as.character(metacell)]

# get heatmap data
fn <- sprintf("Spis_%s_gene_expression.RDS", st)
hmp = readRDS(file.path(sprintf("clustering_%s/scdb/figs/", st), fn))

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
dt <- dt[.N:1]

# save
fwrite(dt, file.path(
    sprintf("clustering_%s/scdb/figs/", st),
    str_replace(fn, ".RDS", ".tsv")
), sep="\t")
