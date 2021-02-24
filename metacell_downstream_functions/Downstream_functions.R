require(metacell)
require(tgconfig)
require(tgstat)
require(tglkmeans)
require(zoo)
require(data.table)
require(dplyr)
require(stringr)
require(ggplot2)
require(scales)
require(RColorBrewer)
require(ComplexHeatmap)
require(rasterpdf)

#' Compute cell type footprint based on a standard metacell->cell type 
#' definition table, using same strategy as mc_fp.
#' 
#' @param input_table data.frame with three columns: metacell, cell_type and color,
#'    or a path to tsv file with this table
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#'  
#' @return analogous mc_object for cell types
#'   
sca_cell_type_fp <- function( input_table, mc_object, mat_object) {
  cell_type_table=read.table(input_table,h=TRUE,sep="\t",comment.char="")
  rownames(cell_type_table)=cell_type_table$metacell
  
  sc_ct_label=as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
  names(sc_ct_label)=names(mc_object@mc)
  
  cells_cols=cell_type_table[as.character(mc_object@mc),"color"]
  cells_cols=as.character(cells_cols)
  names(cells_cols)=names(mc_object@mc)
  
  # filter low expression genes not included in the mc_fp
  umis=as.matrix(mat_object@mat)
  umis=umis[rownames(mc_object@mc_fp),names(mc_object@mc)]
  
  #ct_counts=t(apply(umis,1,function(x) tapply(x, sc_ct_label, sum)))
  #ct_size=colSums(ct_counts)
  #ct_umifrac=t(apply(ct_counts,1,function(x) x*1000/ct_size))
  #ct_umifrac_n=(0.1 + ct_umifrac)/apply(0.1 + ct_umifrac, 1, median)
  
  ct_geomean=t(apply(umis, 1,  function(x) tapply(x, sc_ct_label,function(y) exp(mean(log(1+y)))-1)))
  ct_meansize=tapply(colSums(umis), sc_ct_label, mean)
  ideal_cell_size=pmin(1000,median(ct_meansize))
  g_fp=t(ideal_cell_size * t(ct_geomean)/as.vector(ct_meansize))
  fp_reg=0.05
  g_fp_n=(fp_reg + g_fp)/apply(fp_reg + g_fp, 1, median)
  
  # for compatibility with other functions, return as a MC-like object
  ct_table=mc_object
  ct_table@mc_fp=g_fp_n
  ct_table@mc=sc_ct_label
  ct_table@colors=cells_cols
  ct_table@cell_names=names(sc_ct_label)
  return(ct_table)
}

#' Transfer cell annotations to metacells
#'
#' @param cell_df data.frame with two columns containing cell and cell type IDs;
#'    cell IDs should be the same as in `mc_df`
#' @param mc_df data.frame with two columns: cell and metacell IDs; 
#'    cell IDs should be the same as in `cell_df`
#' @param NA_cell_type logical, should the missing cell type assignment (NA) be kept
#'    among final assigned cell types (default: FALSE)
#' @param niche_order character, metacell ids in the order in which they should be 
#'    plotted (default: NULL)
#' @param output_file png file to save the plot to
#' @param interactive_plot logical, whether to output plotly interactive graph
#'    (default: FALSE). Note that this only works in an interactive session i.e. when
#'    `interactive()==TRUE`
#'    
#' @return data.frame with metacells' cell_type annotations. Also saves the plot to disk.
#' 
sca_transfer_annotations <- function(
  cell_df,  mc_df, NA_cell_type = FALSE,
  niche_order = NULL, color_df = NULL,
  output_file = NULL, width = 4800, height = 2100, res = NA, mc_font_size = 4,
  interactive_plot = FALSE
) {
  
  require(RColorBrewer)
  require(ggplot2)
  require(data.table)
  
  # input data frames as data.table
  dt <- copy(cell_df)[,1:2]
  mc_dt <- copy(mc_df)[,1:2]
  setnames(dt, c("cell", "cell_type"))
  setnames(mc_dt, c("cell", "mc"))
  dt[,cell:=as.character(cell)]
  mc_dt[,cell:=as.character(cell)]
  
  # add mc annotation to cells
  dt[mc_dt, on="cell", mc:=i.mc] 
  # remove cells that are not assigned to any mc
  missing_cells <- is.na(dt$mc)
  dt <- dt[!missing_cells]
  message(sprintf(
    "Removed %s cells (out of %s) that are not in any metacell!",
    sum(missing_cells),length(missing_cells)
  ))
  # for each mc, count how many cells in it are assigned to each ct
  dt[,cells_from_mc_in_ct:=.N,by=.(mc,cell_type)]
  
  # find most occuring ct for each mc
  if (NA_cell_type == TRUE) {
    dt[,max_cells_ct:=max(cells_from_mc_in_ct),mc] 
    dt[,max:={cells_from_mc_in_ct==max_cells_ct}]
    max_ct <- unique(dt[max==TRUE,.(cell_type,mc)])
  } else if (NA_cell_type == FALSE) {
    dtfilt <- dt[!is.na(cell_type)][,max_cells_ct:=max(cells_from_mc_in_ct),mc] 
    dtfilt[,max:={cells_from_mc_in_ct==max_cells_ct}]
    max_ct <- unique(dtfilt[max==TRUE,.(cell_type,mc)])
  }
  dt[max_ct, on="mc", assigned_cell_type:=i.cell_type]
  setorder(dt, assigned_cell_type, mc, na.last = TRUE)
  
  # cell assigned to cell types
  cell_annotation <- unique(dt[,.(cell,mc,cell_type,assigned_cell_type)])
  setnames(cell_annotation, "cell_type", "original_cell_type")
  
  # mc assigned to cell types
  mc_annotation <- unique(dt[,.(mc,assigned_cell_type)])
  setnames(mc_annotation, "assigned_cell_type", "cell_type")
  
  # order metacells by niche order
  if (is.null(niche_order)) 
    niche_order <- mc_annotation$mc
  dt[,mc:=factor(mc, levels=niche_order)]
  cell_annotation[,mc:=factor(mc, levels=niche_order)]
  mc_annotation[,mc:=factor(mc, levels=niche_order)]
  
  if (!is.null(col_df)) {
    col_dt <- copy(col_df)[,1:2]
    setnames(col_dt, c("cell_type","color"))
    cols <- col_dt$color
    names(cols) <- col_dt$cell_type
  } else {
    cpalette <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", " cyan"))
    cts <- unique(dt$cell_type)
    cols <- cpalette(length(cts))
    names(cols) <- cts
  }
  
  # barplot
  
  # plotly doesn't like factors, but this means the ordering of bars will change...  
  if (interactive_plot==TRUE) {
    dt[,cell_type:=as.character(cell_type)]
    dt[,mc:=as.character(mc)]
  }
  
  gp <- ggplot(dt, aes(x=mc, fill=cell_type), color="black") + 
    geom_bar(position = "fill") + 
    scale_x_discrete(breaks=dt$mc, labels=dt$mc) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = cols) +
    theme_minimal() + theme(panel.grid = element_blank()) +
    labs(x="metacells",y="percent of cells in mc") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = mc_font_size))
  
  if (is.null(output_file))
    output_file <- file.path(getwd(), "transfer_annot.png")
  #ggsave(filename = name, plot = gp, width = width)
  png(output_file, width = width, height = height, res = res)
  print(gp)
  dev.off()
  
  if (interactive() & interactive_plot) {
    require(plotly)
    ggplotly(gp)
  }
  
  col_dt2 <- copy(col_dt)
  setnames(col_dt2,"cell_type","assigned_cell_type")
  dt_col <- dt[col_dt, on="cell_type", color:=i.color]
  out_dt <- dt_col[col_dt2, on="assigned_cell_type", assigned_color:=i.color]
  #out_dt <- cell_annotation
  out_df <- copy(out_dt)
  class(out_df) <- "data.frame"
  return(out_df)
  
}


#' Plot barplot of metacell batch distribution
#' 
#' Plots barplot of batch diistribution of cells in each metacell 
#' 

#' @param batch_field charter, name of the column in `mat_object@cell_metadata` 
#'   containg batch information
#' @param clust_order
#' @param update_color_file NULL
#' @param output_file output file name, default `NULL` will save 
#'   'Batch_distribution.pdf' file to initialized scdb figures directory
#' @param plot_format character, format of the output file, either 
#'    "png" (default) or "pdf"
#' @param height numeric, height of the plot saved to output file
#' @param width numeric, width of the plot saved to output file
#' 
scp_batch_distribution <- function(
  mc_id, mat_id, 
  batch_field="dataset", clust_order, update_color_file=NULL,
  plot_format="png", height=600, width=1200
){
  if (!scdb_is_valid()) {
    stop("scdb not initialized!")
  }
  mc_object <- scdb_mc(mc_id)
  mat_object <- scdb_mat(mat_id)
  x=sapply(clust_order,function(x) table(factor(mat_object@cell_metadata[names(which(mc_object@mc==x)),batch_field],levels=unique(mat_object@cell_metadata[,batch_field]))))
  colnames(x)=as.character(clust_order)
  # get batches
  categories=as.character(unique(mat_object@cell_metadata[,batch_field]))
  # check if there's a defined batch color        
  if(is.null(mat_object@cell_metadata$color)){
    n_colors=length(categories)
    color=RColorBrewer::brewer.pal(n_colors,"Dark2")
    names(color)=categories
  } else {
    color=sort(sapply(categories,function(x) as.character(unique(mat_object@cell_metadata[which(mat_object@cell_metadata[,batch_field]==x),"color"]))))
    color=color[rownames(x)]
  }
  if(!is.null(update_color_file)){
    update_color_table=read.table(update_color_file,h=T,row.names=1,sep="\t")
    update_color=as.character(unique(update_color_table$color))
    names(update_color)=unique(update_color_table[,batch_field])
    color=update_color[names(color)]
  }
  # batch frequencies
  x_colSums=colSums(x)
  x_freq=t(apply(x,1,function(y) y/x_colSums))
  bckg_freq=table(mat_object@cell_metadata[,batch_field])/ncol(mat_object@mat)
  FC=apply(x_freq,2,function(y) y/bckg_freq)
  logFC=log2(FC+1)
  # save
  if (is.null(.scfigs_base) | !file.exists(.scfigs_base)) {
    stop("figs directory at ", .scfigs_base, " is missing")
  }
  output_file_base <- sprintf("%s/%s.%s", .scfigs_base, mc_id, "batch_distribution")
  if (plot_format=="png") {
    output_file <- paste0(output_file_base,".png")
    png(output_file,h=height,w=width)
  } else if (plot_format=="pdf") {
    output_file <- paste0(output_file_base,".pdf")
    pdf(output_file,h=height,w=width,useDingbats=FALSE)
  } else {
    stop("plot_format should be either 'png' or 'pdf'")
  }
  barplot(x,las=2,col=color,cex.names=1,cex.axis=2)
  legend("topleft",legend=names(color),fill=color,box.lty=0,horiz=F,ncol=round(length(color)/3,0),cex=1.5)
  dev.off()
  
  message("Saved to", output_file)
}


#' Plot heatmap of gene expression
#' 
#' Plots heatmap of gene expression fold change for metacells (and optionally, 
#' for single cells).
#' 
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param black_list character, blacklisted genes
#' @param output_file output file name without extension, all output files 
#'    will have this prefix; if NULL (default), will not save heatmap to file
#' @param height numeric, height of the heatmap image saved to output file
#' @param width numeric, width of the heatmap image saved to output file
#' @param res numeric, resolution of the output file
#' @param plot_format character, format of the output image file, either 
#'    "png" (default) or "pdf"
#' @param print_heatmap logical, whether to draw heatmap to console
#' @param sub_list_mc charcter, subset of metacells to plot (default: NULL)
#' @param plot_sc logical, whether to also plot a heatmap of single cells 
#'    (default: FALSE)
#' @param plot_sc_width_cex numeric, factor by which to scale single cell plot 
#'    width (default: 2)
#' @param plot_sc_height_cex numeric, factor by which to scale single cell plot 
#'    height (default: 1)
#' @param gene_list character, list of genes to plot, if NULL (default) ...
#' @param order_genes logical, whether to cluster genes (default: TRUE)
#' @param highlight_genes logical or character, whether to highlight top expressed 
#'    and annotted genes on the heatmap; if character, a set of genes which to 
#'    highlight on the heatmap
#' @param gene_annot_file charcter, file path to the file containig gene annotations, 
#'    it should have three tab separated columns containing gene ID, pfam architecture 
#'    and any additional annotation in the last column
#' @param tfs_annot_file character, file path to the file containig TFs annotations, 
#'    minimum required is one column with the gene IDs; if this file is provided, TF  
#'    genes annotations will be highlighted in red on the heatmap
#' @param annot_header logical, gene annotation file has column names?
#' @param show_gene_names logical, show gene names as heatmap row names
#'    (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames 
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_ord character, metacells in the order in which they should be 
#'    plotted; if cluster order is not specified (default: NULL), it is 
#'    determined by hierarchical clustering
#' @param clust_col character, colours to asssign to clusters, should be either 
#'    a named vecor where names are names of metacells, or an unamed vector 
#'    ordered in the same order as specified by clust_ord; if NULL (default), 
#'    cluster colour annotation bar is not printed
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of 
#'    heatmap columns (default:  NULL)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(15,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names 
#'    (default: TRUE) ~NOT IMPLEMENTED~
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param mc_label_cex numeric, scalling of the metacell names on sc plot (default: 2)
#' @param per_clust_genes integer, how many genes per cluster to aim to show in the heatmap 
#'    (default: 20)
#' @param gene_min_fold numeric, minimum fold change for a gene to be considered for plotting
#'    (default: 2)
#' @param smoothen integer, width of the window used to calculate rolling mean 
#'    of expression in single cells (default: 5)
#' @param max_expression_fc,max_expression_fc_sc numeric, max expression value 
#'    to scale metacell and single cell heatmap coloring to (default: 5) 
#' @param transverality_N integer, number of metacells in which a gene can be highly expressed (>1.4) 
#'    to be considered for plotting, by default this is the total number of metacells
#' @param transv_excluded_mc character, metacells to be excluded in transversality calculation 
#'    (default: NULL)
#' 
#' 
scp_plot_cmod_markers <- function(
  mc_object, mat_object, black_list=c(),
  output_file = NULL, plot_format = "png", height = 8000, width = 3000, res = NA, 
  print_heatmap = FALSE, show_heatmap_legend = TRUE,
  sub_list_mc = NULL, plot_sc = FALSE, plot_sc_width_cex = 2, plot_sc_height_cex = 1,
  
  gene_list = NULL, order_genes = TRUE, highlight_genes = NULL,
  gene_annot_file = NULL, tfs_annot_file=NULL, annot_header = FALSE,
  show_gene_names = FALSE, gene_font_size = 4, gene_chr_limit = 100,
  
  clust_ord = NULL, clust_col = NULL, clust_bars = NULL, show_clust_borders=TRUE,
  clust_anno_size = unit(15,"mm"),
  show_mc_names = TRUE, mc_font_size = 4, mc_label_cex = 2,
  column_anno_padding=unit(10,"mm"), row_anno_padding=unit(5,"mm"), dimname_padding=unit(5,"mm"),
  
  per_clust_genes = 20, gene_min_fold = 2, smoothen = 5, max_expression_fc = 5, max_expression_fc_sc = 5,
  transverality_N = ncol(mc_object@mc_fp), transv_excluded_mc = NULL,
  verbose=FALSE
)
{ 
  # load gene annotations
  if(!is.null(gene_annot_file))
    annot <- read.table(gene_annot_file,header=annot_header,sep="\t",fill=TRUE,quote="",row.names=1)
  
  # directory to save output files to
  if (is.null(output_file)) {
    outdir <- getwd()
  } else {
    outdir <- dirname(output_file)
  }
  
  # expression matrix
  if(is.null(sub_list_mc)){
    niche_geomean_n= mc_object@mc_fp
  }else{
    niche_geomean_n= mc_object@mc_fp[,sub_list_mc]
    clust_ord=sub_list_mc
  }
  
  # select genes for plotting
  if(is.null(gene_list)){
    # exclude genes with fc < gene_min_fold
    genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes))))))
    # exclude blacklisted genes
    genes=setdiff(genes, black_list)
    # exclude genes with transversality > transverality_N
    transversal_genes=names(which(
      apply(
        niche_geomean_n[,setdiff(as.character(colnames(niche_geomean_n)),transv_excluded_mc)], 
        1, 
        function(x) sort(x,decreasing=T)[transverality_N]>1.4
      )
    ))
    genes=setdiff(genes, transversal_genes)
  } else {
    # plot only genes in gene list, if it is specified
    genes=gene_list
  }
  
  genes=intersect(genes,rownames(niche_geomean_n))
  message("Will use ",length(genes)," genes")
  
  mat_niche <- niche_geomean_n[genes,]	
  
  # if cluster order is not specified, do hierarchical clustering
  if(is.null(clust_ord)) {
    message("Recomputing cell ord")
    hc1 = hclust(dist(cor(mat_niche,method="pearson")), "ward.D2")
    clust_ord = as.character(hc1$order)
    # write.table(
    #   clust_ord,
    #   file=file.path(outdir, "tmp_cell_clusts_ordered_by_scp_markers_plot.txt"),
    #   quote=FALSE,col.names=FALSE,row.names=FALSE
    # )
    # png(
    #   file.path(outdir,"tmp_cell_clusts_ordered_by_scp_markes_tree.png"),
    #   height=500,width=1000
    # )
    # plot(hc1,xlab="",xaxt='n',hang=-1,ylab="",main="",cex=1)
    # dev.off()  
    scr_tmp_niche_order <- as.character(hc1$order)
  }
  
  # barplot annotation
  if (!is.null(clust_bars)) {
    if (is.null(names(clust_bars)))
      names(clust_bars) <- as.character(clust_ord)
    anno_bar <- clust_bars[as.character(clust_ord)]
  }
  
  # if gene order is TRUE, order genes
  if (order_genes){
    message("Ordering genes")
    gene_ord = order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,1))))
  } else {
    gene_ord= 1:length(genes)
  }
  gene_ord <- rev(gene_ord)
  # write.table(
  #   genes[gene_ord],
  #   file=file.path(outdir,"tmp_markers_ordered.txt"),
  #   quote=FALSE,col.names = FALSE,row.names=FALSE
  # )
  
  # gene labels
  if(!is.null(gene_annot_file)) {
    
    gene_labels_0 <- genes[gene_ord]
    
    message("Genes: ", head(gene_labels_0), "...")
    
    gene_labels_1 <- as.character(annot[genes[gene_ord],2])
    message("Gene labels: ", head(gene_labels_0), "...")
    bad_labels <- gene_labels_1 %in% c("","-"," ") | is.na(gene_labels_1)
    message(sum(bad_labels), " bad gene labels")
    gene_labels_1[bad_labels] <- gene_labels_0[bad_labels]
    long_labels <- nchar(gene_labels_1)>gene_chr_limit
    message(sum(long_labels), " long gene labels")
    gene_labels_1[long_labels] <- paste0(substr(gene_labels_1[long_labels],1,gene_chr_limit-3),"...")
    names(gene_labels_1) <- genes[gene_ord]
    
    gene_labels_3 <- as.character(annot[genes[gene_ord],1])
    bad_labels <- gene_labels_3 %in% c("","-"," ") | is.na(gene_labels_3)
    gene_labels_3[bad_labels] <- genes[gene_ord][bad_labels]
    long_labels <- nchar(gene_labels_3)>gene_chr_limit
    gene_labels_3[long_labels] <- paste0(substr(gene_labels_3[long_labels],1,gene_chr_limit-3),"...")
    gene_labels_2 <- ifelse(gene_labels_0==gene_labels_3, gene_labels_0, paste(gene_labels_0,gene_labels_3,sep=" "))
    names(gene_labels_2) <- genes[gene_ord]
    
  } else {
    
    gene_labels_0 <-genes[gene_ord]
    gene_labels_1 <-genes[gene_ord]
    gene_labels_2 <-genes[gene_ord]
    
  }
  
  gene_font_col <- rep("black",length(gene_labels_0))
  if (!is.null(tfs_annot_file)) {
    tfs_ann <- fread(tfs_annot_file,header=FALSE)
    gene_font_col[gene_labels_0 %in% tfs_ann[[1]]] <- "red"
  }
  
  if (!is.null(highlight_genes)) {
    if (class(highlight_genes)=="logical") {
      if (highlight_genes==TRUE) {
        highlight_genes <- gene_labels_0[gene_labels_1 != gene_labels_0]
      } 
    }
    highlight_ids <- match(highlight_genes,gene_labels_0)
    highlight_genes <- highlight_genes[gene_labels_1[highlight_ids] != gene_labels_0[highlight_ids]]
    
    gids <- which(gene_labels_0 %in% highlight_genes)
    if (length(gids>0)) {
      gene_labels_1[!gene_labels_0 %in% highlight_genes] <- ""
      gene_labels_2[!gene_labels_0 %in% highlight_genes] <- ""
    } else {
      gids <- which(gene_labels_0)
    }
    
  }
  
  
  
  ########################### PLOT METACELL PROFILE ###########################
  
  message("Plotting metacell expression")
  
  # heatmap
  ht_opt(
    RESET = TRUE, 
    COLUMN_ANNO_PADDING=column_anno_padding, 
    ROW_ANNO_PADDING=row_anno_padding, 
    DIMNAME_PADDING=dimname_padding
  )
  mat1 <- pmin(log2(niche_geomean_n[genes[gene_ord],as.character(clust_ord)]+1),max_expression_fc)
  if (show_gene_names) {
    if (length(highlight_genes)>1) {
      message("Gene annots highlights")
      row_ha_right <- HeatmapAnnotation(
        which = "row", simple_anno_size = unit(10,"mm"),
        gene = anno_mark(
          which="row", side="right", at=gids, labels=gene_labels_1[gids], 
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), 
          extend=unit(10, "mm")
        )
      )
      row_ha_left <- HeatmapAnnotation(
        which = "row", simple_anno_size = unit(10,"mm"), 
        gene = anno_mark(
          which="row", side="left", at=gids, labels=gene_labels_2[gids], 
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), 
          extend=unit(10, "mm")
        )
      )
    } else {
      message("Gene annots")
      if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
      row_ha_right <- HeatmapAnnotation(
        which = "row", simple_anno_size = unit(10,"mm"),
        gene = anno_text(which = "row", gene_labels_1, location = 0.5, just = "left", gp = gpar(
          fontsize = gene_font_size, col = gene_font_col
        ))
      )
      if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
      row_ha_left <- HeatmapAnnotation(
        which = "row", simple_anno_size = unit(10,"mm"),
        gene = anno_text(which = "row", gene_labels_2, location = 0.5, just = "right", gp = gpar(
          fontsize = gene_font_size, col = gene_font_col
        ))
      )
    }
  } else {
    row_ha_left <- HeatmapAnnotation(
      which = "row", empty = anno_empty(which = "row", border = FALSE)
    )
    row_ha_right <- row_ha_left
  }
  
  message("Expreession colors...")
  col_fun <- colorRampPalette(c("white","white","orange","red","purple","black"))
  shades <- col_fun(1000)
  
  # mc labels
  message("Metacell labels...")
  collabs <- colnames(mat1) 
  if (!show_mc_names) collabs <- rep("",length(collabs))
  column_lab_ha <- HeatmapAnnotation(
    which = "column",
    LAB = anno_text(which = "column", collabs, gp = gpar(fontsize = mc_font_size, rot=90))
  )
  top_column_ha <- c(column_lab_ha)
  bottom_column_ha <- c(column_lab_ha)
  empty_ha <- HeatmapAnnotation(
    empty = anno_empty(which = "column", border = FALSE, height = unit(5,"mm"))
  )
  
  # cluster colours
  if(!is.null(clust_col)){
    message("Columns...")
    if (is.null(names(clust_col))) {
      names(clust_col) <- clust_ord
    } else {
      if (!all(names(clust_col) %in% clust_ord))
        stop("Colour and cluster names do not match!")
      clust_col <- clust_col[clust_ord]
    }
    column_col_ha <- HeatmapAnnotation(
      which = "column",
      MC = colnames(mat1), col = list(MC = clust_col),
      border = c(TRUE), simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
    )
    top_column_ha <- c(column_col_ha,empty_ha,top_column_ha)
    bottom_column_ha <- c(bottom_column_ha,empty_ha,column_col_ha)
    
  }
  
  # barplots
  if (!is.null(clust_bars)) {
    
    column_bar_ha <- HeatmapAnnotation( 
      which = "column",
      BAR = anno_barplot(
        anno_bar, height = 3*clust_anno_size, bar_width = 0.9, 
        gp = gpar(fill = clust_col, fontsize = mc_font_size),
        axis_param = list(gp = gpar(fontsize = mc_font_size))
      ),
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
    )
    top_column_ha <- c(column_bar_ha,empty_ha,top_column_ha)
    
  } 
  
  # expression heatmap
  h1 <- Heatmap(
    mat1, name = "expression", col = shades, #use_raster = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE, 
    show_column_names = FALSE, show_row_names = FALSE, 
    #row_names_side = "left", row_labels = gene_labels_1, row_names_gp = gpar(fontsize = gene_font_size),
    right_annotation = row_ha_right, left_annotation = row_ha_left,
    column_names_gp = gpar(fontsize = mc_font_size),
    top_annotation = top_column_ha, bottom_annotation = bottom_column_ha, 
    border = TRUE, show_heatmap_legend = show_heatmap_legend
  )
  
  hlist <- h1
  
  if(!is.null(clust_col) & show_clust_borders) {
    mat2 <- rbind(clust_col[match(clust_ord, names(clust_col))])
  } else {
    mat2 <- rbind(clust_ord)
  }
  
  # save figure
  if (!is.null(output_file)) {
    output_file <- stringr::str_remove(output_file,"\\.png$")
    output_file_rds <- sprintf("%s.RDS",output_file)
    output_file_png <- sprintf("%s.png",output_file)
    output_file_pdf <- sprintf("%s.pdf",output_file)
    output_file_legend <- sprintf("%s_legend.png",output_file)
    saveRDS(hlist, output_file_rds)
    
    # # ComplexHeatmap version
    if (plot_format=="png") {
      png(output_file_png, h=height, w=width, res=res)
    } else if (plot_format=="pdf") {
      pdf(output_file_pdf, h=height, w=width, useDingbats=TRUE)
    }
    draw(hlist, padding=unit(c(10, 100, 10, 50), "mm"))
    if (show_clust_borders) {
      change_clust <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
      decorate_heatmap_body("expression", {
        for (i in change_clust) {
          grid.lines(x = i/ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.25))
        }
      })
    }
    dev.off()
    
  }
  
  ######################### PLOT SINGLE-CELL PROFILE ########################## 
  if(plot_sc==TRUE)
  {
    message("Plotting single cell expression")
    cell_order=c()  
    for (niche in clust_ord){
      cells=names(mc_object@mc[which(mc_object@mc==niche)])    
      cell_order=c(cell_order,cells)
    }
    cluster_cell_count=as.matrix(table(mc_object@mc))
    n_cells_cluster=cluster_cell_count[clust_ord,1]
    cells_clusts <- unlist(mapply(rep, clust_ord, n_cells_cluster, USE.NAMES=FALSE))
    
    umis=as.matrix(mat_object@mat[,names(mc_object@mc)])
    mat = umis[genes, cell_order]
    totu = colSums(umis[, cell_order])
    mat = t(t(mat)/totu)*800
    
    lus_1 = log2(1+7*mat[genes[gene_ord], cell_order])
    lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
    lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x, smoothen, fill=0)))
    
    # heatmap
    mat1sc <- pmin(lus_smoo,max_expression_fc_sc)
    colnames(mat1sc) <- colnames(lus_smoo)
    shades <- colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
    col_fun <- colorRampPalette(c("white","white","orange","red","purple","black"))
    shades <- col_fun(1000)
    
    # mc labels
    top_column_ha <- HeatmapAnnotation(
      mclabstop = anno_empty(which = "column", border = FALSE, height = unit(15,"mm"))
    )
    bottom_column_ha <- HeatmapAnnotation(
      mclabsbottom = anno_empty(which = "column", border = FALSE, height = unit(15,"mm"))
    )
    # top_column_ha <- c(sc_mc_labels_ha)
    # bottom_column_ha <- c(sc_mc_labels_ha)
    
    if(!is.null(clust_col)){
      
      column_col_ha <- HeatmapAnnotation(
        which = "column",
        MC = cells_clusts, col = list(MC = clust_col),
        border = c(TRUE), simple_anno_size = 2*clust_anno_size,
        show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
      )
      top_column_ha <- c(column_col_ha,empty_ha,top_column_ha)
      bottom_column_ha <- c(bottom_column_ha,empty_ha,column_col_ha)
      
    }
    
    # expression heatmap
    h1sc <- Heatmap(
      mat1sc, name = "sc_expression", col = shades, use_raster = TRUE, 
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_column_names = FALSE, show_row_names = FALSE, 
      #row_names_side = "left", row_labels = gene_labels_1, row_names_gp = gpar(fontsize = gene_font_size),      show_column_names = FALSE,
      right_annotation = row_ha_right, left_annotation = row_ha_left,
      top_annotation = top_column_ha, bottom_annotation = bottom_column_ha,
      show_heatmap_legend = FALSE, border = TRUE
    )
    hlistsc <- h1sc
    
    # save figure
    if (!is.null(output_file)) {
      output_file <- stringr::str_remove(output_file,"\\.png$")
      output_file_sc_rds <- sprintf("%s_sc.RDS",output_file)
      output_file_sc_png <- sprintf("%s_sc.png",output_file)
      output_file_sc_pdf <- sprintf("%s_sc.pdf",output_file)
      saveRDS(hlistsc, output_file_sc_rds)
      if (plot_format=="png") {
        png(output_file_sc_png, h=height*plot_sc_height_cex, w=width*plot_sc_width_cex, res=res)
      } else if (plot_format=="pdf") {
        pdf(output_file_sc_pdf, h=height*plot_sc_height_cex, w=width*plot_sc_width_cex, useDingbats=TRUE)
      }
      # heatmap
      draw(hlistsc, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
      
      # add grid
      if(!is.null(clust_col)) {
        mat2sc <- rbind(clust_col[match(cells_clusts, names(clust_col))])
      } else {
        mat2sc <- rbind(cells_clusts)
        if (is.null(colnames(mat2sc)))
          colnames(mat2sc) <- cells_clusts
      }
      change_clust_sc <- which(sapply(2:ncol(mat2sc), function(i) mat2sc[,i]!=mat2sc[,i-1]))
      change_mc_sc <- c(
        which(sapply(2:ncol(mat2sc), function(i) colnames(mat2sc)[i]!=colnames(mat2sc)[i-1])),
        ncol(mat2sc)
      )
      if (show_clust_borders) {
        decorate_heatmap_body("sc_expression", {
          for (i in change_clust_sc) {
            grid.lines(x = i/ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 1))
          }
        })
      }
      .add_mc_labels <- function(pos,labs) {
        for (j in 1:length(pos)) {
          i <- pos[j]
          iprev <- ifelse(j==1,0,pos[j-1])
          nt <- ncol(mat2sc)
          grid.text(label = labs[j], x = i/nt-(i/nt-iprev/nt)/2, y = 0.5, gp = gpar(fontsize = mc_label_cex*mc_font_size), rot=90)
        }
      }
      if (show_mc_names) {
        decorate_annotation("mclabstop", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
        decorate_annotation("mclabsbottom", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
      }
      dev.off()
      ht_opt(RESET = TRUE)
    }
    
  }
  ########################## RETURN METACELL PROFILE ##########################
  if (print_heatmap==TRUE) {
    draw(hlist)
  }
  
  message("Done.")
}


#' Comparison of gene expression of two clusterig solutions
#' 
#' @param matrix1,matrix2 gene fold change matrices, with genes in rows and 
#'    clusters in columns; shold have both row names and column names
#' @param markers character, optional list of marker genes to use for 
#'    correlation calculation
#' @param master_gene character
#' @param cor_method character, comparison metric to use, should be one of the following
#'    `c("overlap","jaccard","pearson","kendall","spearman")`
#' @param fc_thrs numeric, lfc threshold to use for binarizing gene expression in clusters
#'    (default: 2); used for calculation when `cor_method` is one of `c("overlap",jaccard")`, 
#'    and for reporting overlapping genes
#' @param output_file output heatmap file name
#' @param width numeric, figure height (in px)
#' @param height numeric, figure width (in px)
#' @param res numeric (default: NA)
#' @param original_ordering_1st
#' @param annotation_file_1
#' @param reorder_by_ann1,reorder_by_ann2 logical, whether to plot metacells in order 
#'   in which they appear in annotation files (default: FALSE)
#' @param cor_color character, colors to use for plotting
#' @param plot_type character, either `"heatmap"` (default) or `"dotplot"`
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size numeric, size of annotation labels
#' @param cor_max numeric, color scaling max value (default: 1)
#' @param cex_dot numeric, dot size scaling factor (default: 1), 
#'    only used if `plot_type=="dotplot"`
#' @param grid logical
#' @param annotation_grid_1,annotation_grid_2 logical
#' 
#' @return a list with following elements: 
#'   1) `heatmap` complex heatmap object
#'   2) `cor_matrix` correlation matrix used for plotting
#'   3) `overlap_matrix` matrix with number of overlapping genes
#'   4) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
sca_clustering_comparison <- function(
  matrix1, matrix2, 
  markers=NULL, master_gene=NULL, 
  cor_method="jaccard", fc_thrs=1.2, 
  output_file=NULL, width=3000, height=3000, res=NA,
  name1=NULL, name2=NULL,
  original_ordering_1st=FALSE, original_ordering_2nd=FALSE, 
  annotation_file_1=NULL, annotation_file_2=NULL,
  reorder_sps2=FALSE, reorder_by_ann1=FALSE, reorder_by_ann2=FALSE,
  cor_color = NULL, plot_type="heatmap", annotation_size = 10, label_font_size = 12, cor_max=1, cex_dot=1,
  grid = TRUE,  annotation_grid_1 = TRUE, annotation_grid_2 = TRUE
) {
  # select marker genes
  if(!is.null(markers)){
    markers=intersect(markers,rownames(matrix1))
    markers=intersect(markers,rownames(matrix2))
    mat1=matrix1[markers,]
    mat2=matrix2[markers,]
  } else {
    mat1=matrix1
    mat2=matrix2
  }
  # binarize expression
  intg <- intersect(rownames(mat1),rownames(mat2))
  mat1[mat1<fc_thrs] <- 0
  mat2[mat2<fc_thrs] <- 0
  mat1[!(mat1<fc_thrs)] <- 1
  mat2[!(mat2<fc_thrs)] <- 1
  mat1 <- data.matrix(mat1)[intg,]
  mat2 <- data.matrix(mat2)[intg,]
  # genes in common
  out_list <- vector("list",length = ncol(mat1))
  for (i in 1:ncol(mat1)) {
    intglist <- lapply(1:ncol(mat2), function(j) {
      ints <- which(mat1[,i] > 0 & mat2[,j] > 0)
      if (length(ints)==0) { NA } else { 
        rownames(mat1)[ints]
      }
    })
    names(intglist) <- colnames(mat2)
    out_list[[i]] <- intglist
  }
  names(out_list) <- colnames(mat1)
  
  # compare clusters
  if (cor_method %in% c("jaccard","overlap")) {
    cl_comp=get(cor_method)(mat1,mat2)
  } else if (cor_method %in% c("pearson","spearman","kendall")) {
    cl_comp=cor(mat1,mat2,method=cor_method)
  }
  rownames(cl_comp)=as.character(colnames(mat1))
  colnames(cl_comp)=as.character(colnames(mat2))
  
  # ordering (if order is not retained)
  if(is.null(master_gene) & !original_ordering_1st & !reorder_by_ann1){
    hc=hclust(dist(cor(t(cl_comp))),method="ward.D2")
    tmp_cor=cl_comp[hc$order,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
    order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
  } else if(original_ordering_1st & !original_ordering_2nd & !reorder_by_ann1){
    hc=list()
    hc$order=colnames(mat1)
    tmp_cor=cl_comp[,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
    order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
  } else if (!is.null(master_gene) & !reorder_by_ann1 & !reorder_by_ann2){
    order_cols=names(sort(matrix2[master_gene,],decreasing=T))
    hc=hclust(dist(cor(t(cl_comp[,order_cols]))),method="ward.D2")
  } else {
    hc=list()
    hc$order=colnames(mat1)
    order_cols=colnames(mat2)
  }
  # plotting matrix
  cor_mat <- cl_comp[hc$order,order_cols]
  cor_name <- sprintf("%s \n",cor_method)
  col_ord <- colnames(cor_mat)
  row_ord <- rownames(cor_mat)
  
  # metacell annotations
  if(!is.null(annotation_file_1)){
    clust_anno_size  <- unit(annotation_size,"mm")
    annr <- read.table(annotation_file_1,header=TRUE)
    if (reorder_by_ann1==TRUE) {
      row_ord <- as.character(annr$metacell)
      cor_mat <- cor_mat[row_ord,]
    } 
    rid <- unlist(lapply(row_ord,match,table=annr$metacell))
    row_clusts <- annr[rid,]$cell_type
    row_clust_col <- annr[rid,]$color
    names(row_clust_col) <- row_clusts
    right_row_col_ha <- HeatmapAnnotation(
      which = "row",
      MC = row_clusts,
      lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
      col = list(MC = row_clust_col),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
    left_row_col_ha <- HeatmapAnnotation(
      which = "row",
      lab = anno_text(which = "row", row_ord, just = "right", location = unit(1, 'npc'), gp = gpar(fontsize = label_font_size)),
      MC = row_clusts,
      col = list(MC = row_clust_col),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
  } else {
    right_row_col_ha <- HeatmapAnnotation(
      which = "row",
      lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
    left_row_col_ha <- HeatmapAnnotation(
      which = "row",
      lab = anno_text(which = "row", row_ord, just = "right", location = unit(1, 'npc'), gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
  }
  
  if(!is.null(annotation_file_2)){
    clust_anno_size  <- unit(annotation_size,"mm")
    annc <- read.table(annotation_file_2, header=TRUE)
    if (reorder_by_ann2==TRUE) {
      col_ord <- as.character(annc$metacell)
      cor_mat <- cor_mat[,col_ord]
    } 
    cid <- unlist(lapply(col_ord,match,table=annc$metacell))
    col_clusts <- annc[cid,]$cell_type
    col_clust_col <- annc[cid,]$color
    names(col_clust_col) <- col_clusts
    top_column_col_ha <- HeatmapAnnotation(
      which = "column",
      lab = anno_text(which = "column", col_ord, just = "left", location = unit(0, 'npc'), gp = gpar(fontsize = label_font_size)),
      MC = col_clusts,
      col = list(MC = col_clust_col),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
    bottom_column_col_ha <- HeatmapAnnotation(
      which = "column",
      MC = col_clusts,
      lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      col = list(MC = col_clust_col),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
  } else {
    top_column_col_ha <- HeatmapAnnotation(
      which = "column",
      lab = anno_text(which = "column", col_ord, just = "left", location = unit(0, 'npc'), gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
    bottom_column_col_ha <- HeatmapAnnotation(
      which = "column",
      lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
  }
  
  # intersect
  cl_int <- overlap(mat1,mat2)
  ovrl_mat <- cl_int[row_ord,col_ord]
  sf <- max(ovrl_mat,na.rm=TRUE)
  cor_mat_sc <- ovrl_mat/sf
  
  # heatmap
  if (is.null(name1)) {
    row_title=""
  } else {
    row_title=name1
  }
  
  if (is.null(name2)) {
    column_title=""
  } else {
    column_title=name2
  }
  
  if (plot_type=="dotplot") {
    if(is.null(cor_color))
      cor_color <- circlize::colorRamp2(c(0, cor_max), c("white", "red"))
    hm <- Heatmap(
      cor_mat, col = cor_color, name = cor_name, border = TRUE, rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (grid == TRUE)
          grid.rect(
            x = x, y = y, width = width, height = height, 
            gp = gpar(col = "gray50", fill = NA, lty = 1, lwd = 0.1)
          )
        grid.circle(
          x = x, y = y, r = abs(cor_mat_sc[i, j])/2 * max(unit.c(width, height)) * cex_dot, 
          gp = gpar(fill = cor_color(cor_mat[i, j]), col = "white")
        )
      },
      cluster_rows = FALSE, cluster_columns = FALSE, 
      show_row_names = FALSE, show_column_names = FALSE,
      row_title = row_title, column_title = column_title,
      right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,  
      top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
      heatmap_legend_param = list(
        #col_fun = cor_color, at = cols_range,
        title = cor_name, border = TRUE,
        legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
        title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
        labels_gp = gpar(fontsize = label_font_size)
      )
    )
    
  } else if (plot_type=="heatmap") {
    if(is.null(cor_color))
      cor_color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
    hm <- Heatmap(
      pmin(cor_mat,cor_max), col = cor_color, name = cor_name, border = TRUE,
      rect_gp = gpar(col = "gray50", lwd = 0.2), 
      cluster_rows = FALSE, cluster_columns = FALSE, 
      show_row_names = FALSE, show_column_names = FALSE,
      row_title = row_title, column_title = column_title,
      right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,  
      top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
      heatmap_legend_param = list(
        title = cor_name, border = TRUE,
        legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
        title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
        labels_gp = gpar(fontsize = label_font_size)
      )
    )
  }
  
  if (!is.null(output_file)) {
    
    png(output_file,h=height,w=width,res=res)
    ht_opt(
      COLUMN_ANNO_PADDING=unit(15,"mm"), 
      ROW_ANNO_PADDING=unit(15,"mm"), 
      DIMNAME_PADDING=unit(5,"mm")
      #HEATMAP_LEGEND_PADDING=unit(10,"mm")
    )
    draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
    
    if(!is.null(annotation_file_1) & annotation_grid_1==TRUE){
      mat2 <- rbind(row_clust_col)
      change_clust_row <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
      decorate_heatmap_body(cor_name, {
        for (i in change_clust_row) {
          grid.lines(x = c(0,1), y = 1-i/ncol(mat2),  gp = gpar(col = "gray50", lwd = 0.05))
        }
      })
    }
    if(!is.null(annotation_file_2) & annotation_grid_2==TRUE){
      mat2 <- rbind(col_clust_col)
      change_clust_col <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
      decorate_heatmap_body(cor_name, {
        for (i in change_clust_col) {
          grid.lines(x = i/ncol(mat2), y = c(0,1), gp = gpar(col = "gray50", lwd = 0.05))
        }
      })
    }
    
    dev.off()
  }
  
  ht_opt(RESET=TRUE)
  
  return(list(
    heatmap=hm, comp=cor_mat, overlap_mat=ovrl_mat, overlapping_genes=out_list
  ))
  
}

#' Correlation between gene expression of two clusterig solutions
#'
#' @param matrix1,matrix2 matrices, with genes in rows and clusters in columns; 
#'    shold have both row names and column names
#' @param output_file output heatmap file name
#' @param markers character, optional list of marker genes to use for 
#'    correlation calculation
#' @param master_gene
#' 
sca_clustering_correlation <- function(
  matrix1, matrix2, output_file, markers=NULL, master_gene=NULL, 
  pmin=1, pmax=0, cor_method="pearson",
  original_ordering_1st=FALSE, original_ordering_2nd=FALSE, 
  annotation_file_1=NULL, annotation_file_2=NULL,
  width=3000, height=3000, 
  annotation_size = 10,  label_font_size = 12
  # cex_row=1, cex_column=1, x_lab_rot=FALSE
){
  
  if(!is.null(markers)){
    markers=intersect(markers,rownames(matrix1))
    markers=intersect(markers,rownames(matrix2))
    mat1=matrix1[markers,]
    mat2=matrix2[markers,]
  }
  if(is.null(markers)){
    mat1=matrix1
    mat2=matrix2
  }
  
  #colnames(scr_markers_km$centers)
  cl_cor=cor(mat1,mat2,method=cor_method)
  rownames(cl_cor)=as.character(colnames(mat1))
  colnames(cl_cor)=as.character(colnames(mat2))
  
  
  
  if(is.null(master_gene) & !original_ordering_1st){
    hc=hclust(dist(cor(t(cl_cor))),method="ward.D2")
    tmp_cor=cl_cor[hc$order,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
    order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
  } else if(original_ordering_1st & !original_ordering_2nd){
    hc=list()
    hc$order=colnames(mat1)
    tmp_cor=cl_cor[,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
    order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
  }else if(original_ordering_1st & original_ordering_2nd){
    hc=list()
    hc$order=colnames(mat1)
    order_cols=colnames(matrix2)
  } else if (!is.null(master_gene)){
    order_cols=names(sort(matrix2[master_gene,],decreasing=T))
    hc=hclust(dist(cor(t(cl_cor[,order_cols]))),method="ward.D2")
  }
  
  cor_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  ht_opt(
    COLUMN_ANNO_PADDING=unit(15,"mm"), 
    ROW_ANNO_PADDING=unit(15,"mm"), 
    DIMNAME_PADDING=unit(5,"mm"),
    HEATMAP_LEGEND_PADDING=unit(15,"mm")
  )
  
  cor_mat <- (pmin(pmax(cl_cor[hc$order,order_cols],pmax),pmin))
  cor_name <- sprintf("%s \n",cor_method)
  
  col_ord <- colnames(cor_mat)
  row_ord <- rownames(cor_mat)
  
  if(!is.null(annotation_file_1)){
    clust_anno_size  <- unit(annotation_size,"mm")
    annr <- read.table(annotation_file_1,header=TRUE)
    rid <- unlist(lapply(row_ord,match,table=annr$metacell))
    row_clusts <- annr[rid,]$cell_type
    row_clust_col <- annr[rid,]$color
    names(row_clust_col) <- row_clusts
    row_col_ha <- HeatmapAnnotation(
      which = "row",
      MC = row_clusts, col = list(MC = row_clust_col),
      lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
  } else {
    row_col_ha <- HeatmapAnnotation(
      which = "row",
      lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
  }
  
  
  if(!is.null(annotation_file_2)){
    clust_anno_size  <- unit(annotation_size,"mm")
    annc <- read.table(annotation_file_2, header=TRUE)
    cid <- unlist(lapply(col_ord,match,table=annc$metacell))
    col_clusts <- annc[cid,]$cell_type
    col_clust_col <- annc[cid,]$color
    names(col_clust_col) <- col_clusts
    column_col_ha <- HeatmapAnnotation(
      which = "column",
      MC = col_clusts, col = list(MC = col_clust_col),
      lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
  } else {
    column_col_ha <- HeatmapAnnotation(
      which = "column",
      lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
  }
  
  hm <- Heatmap(
    cor_mat, col = cor_color, name = cor_name, border = TRUE,
    rect_gp = gpar(col = "gray50", lwd = 0.2), 
    cluster_rows = FALSE, cluster_columns = FALSE, 
    show_row_names = FALSE, show_column_names = FALSE,
    right_annotation = row_col_ha, left_annotation = row_col_ha,  
    top_annotation = column_col_ha, bottom_annotation = column_col_ha,
    heatmap_legend_param = list(
      title = cor_name, border = TRUE,
      legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
      title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
      labels_gp = gpar(fontsize = label_font_size)
    )
  )
  
  png(output_file,h=height,w=width)
  draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
  if(!is.null(annotation_file_1)){
    mat2 <- rbind(row_clust_col)
    change_clust_row <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
    decorate_heatmap_body(cor_name, {
      for (i in change_clust_row) {
        grid.lines(x = c(0,1), y = 1-i/ncol(mat2), gp = gpar(col = "gray50", lwd = 0.2))
      }
    })
  }
  if(!is.null(annotation_file_2)){
    mat2 <- rbind(col_clust_col)
    change_clust_col <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
    decorate_heatmap_body(cor_name, {
      for (i in change_clust_col) {
        grid.lines(x = i/ncol(mat2), y = c(0,1), gp = gpar(col = "gray50", lwd = 0.2))
      }
    })
  }
  
  dev.off()
  ht_opt(RESET=TRUE)
  
  return(cl_cor[hc$order,order_cols])
  
}

