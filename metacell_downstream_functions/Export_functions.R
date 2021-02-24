require(metacell)
require(Seurat)
require(data.table)
require(stringr)
require(ggplot2)

#' Export Metacell solution as Seurat object
#' 
#' @param mat_id character, id of single cell matrix object in scdb
#' @param mc_id character, id of metacell object in scdb
#' @param mc_ann character, path to file containin metacell annotation, this should be a tsv file
#'   with the three columns containing metacell, cell type and color annotaions
#' @param gset_id id of markers gene set in scdb (used for building Metacell solution)
#'   if this is not specified, `nfeatures` will be selected by Seurat
#' @param nfeatures integer, number of variable features to select (using Seurat), 
#'   this doesn't have to be specified if `gset_id` is used (default: NULL)
#' @param dims dimensions of reduction to use as input, passed to `Seurat::FindNeighbors` 
#'   and `Seurat::RunUMAP`
#' @param out_dir character, where to save the object and plots (default: working directory)
#' @param verbose logical, print detailed log messages (default: TRUE)
#' @param ... Other arguments passed to Seurat methods
#' 
#' @return Seurat object of class `Assay`. Also saves 2D projection plots to `out_dir`.
#' 
mc_export_seurat <- function(
  mat_id, mc_id, mc_ann, gset_id, nfeatues=NULL, dims=1:10, out_dir=".", verbose=TRUE
) {
  
  # TBA check that scdb is initialized...
  # TBA create object w/o cell type annotations
  
  dir.create(out_dir)
  
  # seurat object
  if (verbose) message("Loading matrix data")
  mat <- scdb_mat(mat_id)@mat
  if (verbose) message("Gene names must not contain '_', replacing them with '-'.")
  rownames(mat) <- stringr::str_replace_all(rownames(mat),"_","-")
  if (verbose) message("Loading metacell data")
  mc <- scdb_mc(mc_id)
  mc_cells <- intersect(colnames(mat),names(mc@mc))
  if (verbose) message("Creating Seurat object with ",length(mc_cells)," cells.")
  seurat_obj <- CreateSeuratObject(counts=mat[,mc_cells])
  
  # metadata
  seurat_obj@meta.data$orig.ident=mc_id
  if (!is.null(mc_ann)) {
    if (verbose) message("Loading metacell annotations")
    ann <- fread(mc_ann,select=1:3,header=TRUE)
    setnames(ann,c("metacell","cell_type","color"))
    cell_asgn <- data.table(cell=names(mc@mc),metacell=mc@mc)
    cell_ann <- merge.data.table(cell_asgn,ann,by="metacell",all.x=TRUE)
    cell_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    meta_dt <- data.table(seurat_obj@meta.data,keep.rownames="cell")
    meta_dt_ann <- merge.data.table(meta_dt,cell_ann,by="cell",all.x=TRUE)
    meta_dt_ann[is.na(cell_type),c("cell_type","color"):=list("unassigned","gray")]
    meta <- copy(meta_dt_ann); class(meta) <- "data.frame"; rownames(meta) <- meta_dt$cell
    seurat_obj@meta.data <- meta
  } else {
    stop("You have to specify metacell annotation file!")
  }
  
  # variable features and scaling
  if (!is.null(gset_id)) {
    marker_gset <- scdb_gset(gset_id)
    marker_genes <- names(marker_gset@gene_set)
    if (verbose) message("Gene names must not contain '_', replacing them with '-'.")
    marker_genes <- stringr::str_replace_all(marker_genes,"_","-")
    variable_genes <- intersect(marker_genes,rownames(mat))
    if (verbose) message("Supplied ",length(marker_genes)," marker genes; using ",length(variable_genes)," variable genes.")
    VariableFeatures(seurat_obj) <- marker_genes
  } else if (!is.null(nfeatures)) {
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=3000)
    if (verbose) message("Selecting ",nfeatures," variable genes.")
  } else {
    stop("You have to specify either gset_id or nfeatures!")
  }
  if (verbose) message("Scaling data.")
  seurat_obj <- ScaleData(seurat_obj)
  
  # clustering
  if (verbose) message("Running PCA.")
  seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(object=seurat_obj), verbose=verbose)
  if (verbose) message("Finding clusters.")
  seurat_obj <- FindNeighbors(seurat_obj, dims=dims, verbose=verbose)
  seurat_obj <- FindClusters(seurat_obj, resolution=0.5, verbose=verbose)
  if (verbose) message("Running UMAP.")
  seurat_obj <- RunUMAP(seurat_obj, dims=dims, verbose=verbose)
  
  # saving outputs
  if (verbose) message("Plotting 2d projection.")
  gp_col_dt <- unique(meta_dt_ann[,.(cell_type,color)])
  if (verbose) print(gp_col_dt)
  gp_col <- gp_col_dt$color
  names(gp_col) <- gp_col_dt$cell_type
  plots <- lapply(c("pca","umap"), function(x) 
    DimPlot(seurat_obj, reduction=x, group.by="cell_type", combine=TRUE) + 
      theme(legend.position="top") + labs(title = x) +
      scale_color_manual(values=gp_col) +
      guides(color=guide_legend(nrow=3, byrow=TRUE, override.aes=list(size=3))) 
  )
  if (verbose) message("Saving plots.")
  seurat_pdf <- file.path(out_dir,sprintf("seurat_%s.pdf",mc_id))
  pdf(seurat_pdf,height=6,width=12)
  print(plots)
  dev.off()
  seurat_rds <- file.path(out_dir,sprintf("seurat_%s.RDS",mc_id))
  if (verbose) message("Saving Seurat object to ",seurat_rds)
  saveRDS(seurat_obj, file=seurat_rds)
  if (verbose) message("Done.")
  
  return(seurat_obj)
  
}
