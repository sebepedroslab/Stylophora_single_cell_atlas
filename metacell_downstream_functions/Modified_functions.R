require(tgconfig)
require(metacell)

#' Plot color bar
#'
plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=NULL)
{
  if (!is.null(fig_fn)) {
    .plot_start(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

  if (is.null(show_vals_ind)) {
    show_vals_ind = rep(T, length(cols))
  }
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
}


#' Plot 2d projection with more options
#' 
mcell_mc2d_plot_mas = function(
  mc2d_id, legend_pos="topleft", fn_suf="mc", 
  show_mc=TRUE, show_mcid=TRUE, colors=NULL,
  plot_edges=TRUE, min_edge_l=0, edge_w = 1, short_edge_w=0, 
  cell_outline=F, sc_colors=NULL, sc_cex=1, sc_alpha=1, mcp_2d_id_cex=NULL,
  filt_mc=NULL, plot_format="png"
)
{
  mcp_2d_height = get_param("mcell_mc2d_height",package="metacell")
  mcp_2d_width = get_param("mcell_mc2d_width",package="metacell")
  mcp_2d_plot_key = get_param("mcell_mc2d_plot_key",package="metacell")
  mcp_2d_cex = get_param("mcell_mc2d_cex",package="metacell")
  if (is.null(mcp_2d_id_cex)) mcp_2d_id_cex = mcp_2d_cex
  mcp_2d_legend_cex = get_param("mcell_mc2d_legend_cex",package="metacell")
  
  mc2d = scdb_mc2d(mc2d_id)
  if(is.null(mc2d)) {
    stop("missing mc2d when trying to plot, id ", mc2d_id)
  }
  mc = scdb_mc(mc2d@mc_id)
  if(is.null(mc)) {
    stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
  }
  if(!is.null(filt_mc)) {
    f_sc = filt_mc[mc@mc[names(mc2d@sc_x)]]
    mc2d@sc_x[!f_sc] = NA
    mc2d@sc_y[!f_sc] = NA
    mc2d@mc_x[!filt_mc] = NA
    mc2d@mc_y[!filt_mc] = NA
  }
  fig_nm = scfigs_fn(fn_suf, ifelse(plot_edges, "2d_graph_proj", "2d_proj"), ext=plot_format)
  #.plot_start(fig_nm, w=mcp_2d_width, h = mcp_2d_height)
  if (plot_format=="png") {
    png(fig_nm, width = mcp_2d_width, height = mcp_2d_height)
  } else if (plot_format=="pdf") {
    pdf(fig_nm, width = mcp_2d_width/72, height = mcp_2d_height/72)
  } else {
    stop("plot_format should be either 'png' or 'pdf'")
  }
  par(bg=NA)
  if(is.null(colors)) {
    cols = mc@colors
  } else {
    cols = colors
  }
  cols[is.na(cols)] = "gray"
  if (is.null(sc_colors)) {
    sc_colors=cols[mc@mc[names(mc2d@sc_x)]]
  } else {
    if (!is.null(names(sc_colors)))
      sc_colors <- sc_colors[names(mc2d@sc_x)]
  }
  if(cell_outline) {
    raster::plot(mc2d@sc_x, mc2d@sc_y, pch=21, bg=alpha(sc_colors,sc_alpha), cex=sc_cex, lwd=0.5,
         axes=FALSE, frame.plot=TRUE, xlab="", ylab="")
  } else {
    raster::plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=alpha(sc_colors,sc_alpha), cex=sc_cex,
         axes=FALSE, frame.plot=TRUE, xlab="", ylab="")
  }
  if(show_mc) {
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    if (plot_edges) {
      dx = mc2d@mc_x[fr]-mc2d@mc_x[to]
      dy = mc2d@mc_y[fr]-mc2d@mc_y[to]
      f = sqrt(dx*dx+dy*dy) > min_edge_l
      segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
               lwd=ifelse(f, edge_w, short_edge_w))
    }  
    points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=cols)
    if(show_mcid) {
      text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_id_cex)
    }  
  }
  if(nrow(mc@color_key)!=0 & mcp_2d_plot_key) {
    key = mc@color_key[ mc@color_key$color %in% mc@colors, ]
    #		if(nrow(key!=0)) {
    if(!is.null(key) & is.vector(key) & nrow(key) != 0) {
      #group	gene	color	priority	T_fold
      gmark = tapply(key$gene, key$group, paste, collapse=", ")
      gcol = unique(data.frame(col=key$color, group=key$group))
      rownames(gcol) = gcol$group
      if(is.vector(gmark)) {
        gmark = gmark[order(names(gmark))]
      }
      if(legend_pos == "panel") {
        dev.off()
        fig_nm = scfigs_fn(mc2d_id, "2d_proj_legend")
        pdf(fig_nm, width = 600/72, height=(length(gmark)*40+400)/72)
        plot.new()
        legend_pos = "topleft"
      }
      legend(legend_pos,
             legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)),
             pch=19, cex=mcp_2d_legend_cex,
             col=as.character(gcol[names(gmark), 'col']), bty='n')
    }
  }
  
  dev.off()
}

#' Plot 2d projection coloured by values
#'
mcell_mc2d_plot_values <- function (
    mc2d_id, mc_values=NULL, sc_values=NULL, name="", show_mc_ids = F, show_legend = T, neto_points = F, 
    max_v = NA, min_v = NA, color_cells = F,
    zero_sc_v = 0, one_sc_v = 1, two_sc_v = 2, filt_mc = NULL, 
    plot_format="png", raster=TRUE, verbose=FALSE) 
{
    height = get_param("mcell_mc2d_gene_height",package="metacell")
    width = get_param("mcell_mc2d_gene_width",package="metacell")
    mc_cex = get_param("mcell_mc2d_gene_mc_cex",package="metacell")
    sc_cex = get_param("mcell_mc2d_gene_cell_cex",package="metacell")
    colspec = get_param("mcell_mc2d_gene_shades",package="metacell")
    if (is.na(max_v)) {
        max_v = get_param("mcell_mc2d_gene_max_v",package="metacell")
        min_v = -max_v
    }
    mc2d = scdb_mc2d(mc2d_id)
    if (is.null(mc2d)) {
        stop("missing mc2d when trying to plot, id ", mc2d_id)
    }
    mc = scdb_mc(mc2d@mc_id)
    if (is.null(mc)) {
        stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
    }
    if (!is.null(filt_mc)) {
        f_sc = filt_mc[mc@mc[names(mc2d@sc_x)]]
        mc2d@sc_x[!f_sc] = NA
        mc2d@sc_y[!f_sc] = NA
        mc2d@mc_x[!filt_mc] = NA
        mc2d@mc_y[!filt_mc] = NA
    }
    if (is.null(mc_values) & !is.null(sc_values)) {
       mc_values <- unlist(tapply(names(mc@mc),mc@mc, FUN=function(x) median(sc_values[x]),simplify=FALSE),use.names=FALSE) 
       if (verbose) message("mc_values: ", paste(mc_values,collaps=" "))
    } else if (is.null(mc_values) & is.null(sc_values)) {
      stop("You need to specify either mc_values or sc_values, or both!")
    }
    x = pmin(pmax(log2(mc_values), min_v), max_v) - 
        min_v
    shades = colorRampPalette(colspec)(100 * (max_v - min_v) + 1)
    mc_cols = shades[round(100 * x) + 1]
    if (verbose) message("mc_cols: ", paste(mc_cols,collaps=" "))
    fig_nm = scfigs_fn(mc2d_id, sub("\\/", "", name), sprintf("%s/%s.values", 
        .scfigs_base, mc2d_id))
    if (plot_format=="png") {
      png(paste0(fig_nm,".png"), w = width * ifelse(show_legend & !neto_points, 1.25, 1), h = height, res=200)
    } else if (plot_format=="pdf") {
      if(raster==TRUE){
        rasterpdf::raster_pdf(paste0(fig_nm,".pdf"), w = width * ifelse(show_legend & !neto_points, 1.25, 1) / 72, h = height / 72)
      } else {
        pdf(paste0(fig_nm,".pdf"), w = width * ifelse(show_legend & !neto_points, 1.25, 1) / 72, h = height / 72)
      }
    } else {
      stop("plot_format should be either 'png' or 'pdf'")
    }
    if (show_legend & !neto_points) {
        layout(matrix(c(1, 1:3), nrow = 2, ncol = 2), widths = c(4, 
            1))
    }
    if (neto_points) {
        par(mar = c(1, 1, 1, 1))
    }
    else {
        par(mar = c(4, 4, 4, 1))
    }
    sc_cols = "gray80"
    if (color_cells & !is.null(sc_values)) {
        cnms = intersect(names(mc2d@sc_x), names(sc_values))
        sc_umi = rep(NA, length(mc2d@sc_x))
        names(sc_umi) = names(mc2d@sc_x)
        sc_umi[cnms] = sc_values[cnms]
        sc_umi[is.na(sc_umi)] = 0
        base_shade = 1 + floor(length(shades) * max_v/(max_v - 
            min_v))
        l_shade = length(shades) - base_shade - 1
        collow = shades[base_shade + floor(l_shade/4)]
        colmid = shades[base_shade + floor(l_shade/2)]
        colhigh = shades[base_shade + floor(3 * l_shade/4)]
        sc_cols = ifelse(sc_umi <= zero_sc_v, "gray80", ifelse(sc_umi <= 
            one_sc_v, collow, ifelse(sc_umi <= two_sc_v, colmid, 
            colhigh)))
       if (verbose) message("sc_cols: ", paste(unique(sc_cols),collaps=" "))
    }
    plot(mc2d@sc_x, mc2d@sc_y, pch = 19, cex = sc_cex, col = sc_cols, 
        xlab = "", ylab = "", main = ifelse(neto_points, "", 
            name), cex.main = mc_cex, bty = ifelse(neto_points, 
            "n", "o"), xaxt = ifelse(neto_points, "n", "s"), 
        yaxt = ifelse(neto_points, "n", "s"))
    points(mc2d@mc_x, mc2d@mc_y, pch = 21, bg = mc_cols, cex = mc_cex)
    if (show_mc_ids) {
        text(mc2d@mc_x, mc2d@mc_y, seq_along(mc2d@mc_y), cex = mc_cex * 
            0.5)
    }
    if (show_legend & !neto_points) {
        par(mar = c(4, 1, 4, 1))
        plot_color_bar(seq(min_v, max_v, l = length(shades)), 
            shades, show_vals_ind = c(1, 100 * max_v + 1, 200 * 
                max_v + 1))
    }
    dev.off()
}

#' Calculate confusion matrix - needed for plotting function
#' 
mc_compute_norm_confu_matrix <- function(
  mc_id, graph_id, max_deg=NULL
){
  cgraph <- scdb_cgraph(graph_id)
  if(is.null(max_deg))	max_deg <- nrow(cgraph@edges)
  confu <- mcell_mc_confusion_mat(mc_id, graph_id, max_deg,ignore_mismatch=T)
  r_confu <- rowSums(confu)
  c_confu <- colSums(confu) 
  norm <- r_confu %*% t(c_confu)
  confu_n <- confu/norm
  confu_nodiag <- confu_n
  diag(confu_nodiag) <- 0
  confu_n <- pmin(confu_n, max(confu_nodiag))
  confu_n <- pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
  return(confu_n)
}

#' Cluster confusion matrix - needed for plotting function
#' 
mc_confusion_clustering <- function(
  confu_n, clust_method="average"
){
  epsilon <- quantile(confu_n[confu_n != 0], 0.02)
  hc <- hclust(as.dist(-log10(epsilon + confu_n)), clust_method)
  return(hc)
}

#' Plot confusion matrix and save as raster pdf
#'
mcell_mc_plot_confusion_mas <- function (
  mc_id, graph_id, coc_id = NULL, use_orig_order = F,
  mc_order = NULL, fig_fn = NULL, 
  plot_format="png", raster=TRUE, w=7, h=7)
{
  mc = scdb_mc(mc_id)
  if (is.null(mc)) {
    stop("undefined meta cell object ", mc_id)
  }
  if (is.null(fig_fn)) {
    fig_fn = scfigs_fn(mc_id, sprintf("graph%s_confusion", graph_id), ext = plot_format)
  }
  if (!is.null(graph_id)) {
    if (!is.null(coc_id)) {
      stop("cannot specify both a graph and coclust graph when plotting confusion")
    }
    cgraph = scdb_cgraph(graph_id)
    if (is.null(cgraph)) {
      stop("undefined cgraph object when trying to plot confusion, id ",
           graph_id)
    }
    max_deg = nrow(cgraph@edges)
    confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg,
                                   ignore_mismatch = T)
  }
  else if (!is.null(coc_id)) {
    coc = scdb_coclust(coc_id)
    if (is.null(coc)) {
      stop("undefined coclust object when trying to plot confusion, id ",
           coc_id)
    }
    max_deg = median(table(mc@mc))/2
    confu = mcell_mc_coclust_confusion_mat(mc_id, coc_id = coc_id,
                                           K = max_deg, ignore_mismatch = T, alpha = 2)
  }
  r_confu = rowSums(confu)
  c_confu = colSums(confu)
  norm = r_confu %*% t(c_confu)
  confu_n = confu/norm
  colors = mc@colors
  if (!use_orig_order) {
    if (is.null(mc_order)) {
      hc = mcell_mc_hclust_confu(mc_id, graph_id = NULL,
                                 confu)
      mc_order = hc$order
    }
    confu = confu[mc_order, mc_order]
    confu_n = confu_n[mc_order, mc_order]
    colors = colors[mc_order]
    colnames(confu_n) = (1:ncol(confu_n))[mc_order]
    rownames(confu_n) = (1:ncol(confu_n))[mc_order]
  }
  colors[is.na(colors)] = "gray"
  shades = colorRampPalette(c("white", "pink", "red", "black",
                              "brown", "orange"))

  if (plot_format=="png") {
    png(fig_fn, width = w, height = h)
  } else if (plot_format=="pdf") {
    if(raster==TRUE){
      raster_pdf(fig_fn, width = w, height = h, units="in")
    } else {
      pdf(fig_fn, width = w, height = h, units="in")
    }
  } else {
    stop("plot_format should be either 'png' or 'pdf'")
  }

  layout(matrix(c(1, 4, 2, 3), nrow = 2), heights = c(800,
                                                      50), width = c(50, 800))
  tl_marg = c(0, 2, 5, 0)
  par(mar = tl_marg)
  n_mc = ncol(mc@mc_fp)
  image(t(as.matrix(1:n_mc, nrow = 1)), col = colors, yaxt = "n",
        xaxt = "n")
  top_marg = c(0, 0, 5, 5)
  par(mar = top_marg)
  log_scale = F
  if (log_scale) {
    image(log2(1 + confu), col = shades(1000), xaxt = "n",
          yaxt = "n")
  }
  else {
    confu_nodiag = confu_n
    diag(confu_nodiag) = 0
    confu_n = pmin(confu_n, max(confu_nodiag))
    confu_n = pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
    image(confu_n, col = shades(1000), xaxt = "n", yaxt = "n")
  }
  lower_marg = c(3, 0, 0, 5)
  par(mar = lower_marg)
  image(as.matrix(1:n_mc, nrow = 1), col = colors, yaxt = "n",
        xaxt = "n")
  dev.off()
}

