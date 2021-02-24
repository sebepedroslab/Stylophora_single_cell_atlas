require(data.table)
require(stringr)
require(ape)
require(tidytree)
require(ggtree)
require(treeio)
require(dendextend)
require(phangorn)
require(phytools)
require(ComplexHeatmap)
require(foreach)
require(doParallel)

#' Get tips labels in order in which they appear in the tree
#' 
#' @param tree object of class `phylo`
#' 
tipLabels <- function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}

#' Build trees by downsampling variable genes
#' 
treeDownsampling=function(x,vargenes,p=0.75,method="hclust") {
  if (method=="hclust") {
    as.phylo(hclust(as.dist(1-cor(x[sample(vargenes,p*length(vargenes)),],method="pearson")), "average")) 
  }else if (method=="nj") {
    as.phylo(nj(as.dist(1-cor(x[sample(vargenes,p*length(vargenes)),],method="pearson")))) 
  }
}

# Compare two trees
compareTrees <- function(tree1_fn, tree2_fn,outfile_prefix=NULL, labelcols=NULL) {
  # load trees
  ext1 <- tools::file_ext(tree1_fn)
  ext2 <- tools::file_ext(tree2_fn)
  if (ext1=="beast") {
    tree1=as.phylo(read.beast(tree1_fn))
  } else if (ext1=="RDS") {
    tree1=readRDS(tree1_fn)
  } else {
    stop("uNrecognized format of tree1_fn - should be either .beast or .RDS")
  }
  if (ext2=="beast") {
    tree2=as.phylo(read.beast(tree2_fn))
  } else if (ext2=="RDS") {
    tree2=readRDS(tree2_fn)
  } else {
    stop("uNrecognized format of tree2_fn - should be either .beast or .RDS")
  }
  # compare trees
  compare=phytools::cophylo(tree1, tree2)
  tree1_name=str_remove(basename(tree1_fn),"\\.((RDS)|(beast))")
  tree2_name=str_remove(basename(tree2_fn),"\\.((RDS)|(beast))")
  # save trees
  if (is.null(outfile_prefix))
    outfile_prefix="compare"
  out_fn <- sprintf("%s_%s_%s.pdf",outfile_prefix,tree1_name,tree2_name)
  pdf(out_fn,height=18,width=16,useDingbats=TRUE)
  par(mar=rep(2, 4))
  plot(
    compare,
    link.type="curved",link.lwd=4,link.lty="solid",link.col=make.transparent("red",0.25)
  )
  title(main=sprintf("%s vs %s",tree1_name,tree2_name))
  if (!is.null(labelcols)) {
    tiplabels.cophylo(pch=19,cex=3,col=labelcols[compare$trees[[1]]$tip.label],which="left")
    tiplabels.cophylo(pch=19,cex=3,col=labelcols[compare$trees[[2]]$tip.label],which="right")
  }
  dev.off()
}


#' Build tree by ensemble clustering using downsampling variable genes
#'
#' @param x matrix, genes in rows, cell types (or metacells) in columns
#' @param n numer of iterations of clutering that will be performed 
#'   to calculate coocurence
#' @param k numeric, number of clusters into which to cut trees 
#'   to get clusters for which the co-occurrence will be calculated
#' @param h numeric, height at which to cut trees to get clusters for
#'   which the co-occurrence will be calculated
#' @param vargenes character, optinal variable genes for sampling
#' @param p numeric, between 0 and 1, fraction of genes for sampling
#' @param clustering_method character, method to use for clustering, 
#'   for available options see method argument in `?hclust()`
#' @param cor_method character, see method argument in `?cor()`
#' 
treeFromEnsembleClustering=function(
  x, n=1000, k=NULL, h=NULL, vargenes=NULL, p=0.75, bootstrap=FALSE,
  clustering_algorithm="hclust", clustering_method="average",cor_method="pearson"
) {
  ensemble_hc=vector(mode="list",length=n)
  hs <- rep(h,n)[1:n]
  if (clustering_algorithm=="hclust") {
    for (i in 1:n) { 
      if (bootstrap==TRUE) p=1
      ids=sample(vargenes,p*length(vargenes),replace=bootstrap)
      hc=hclust(as.dist(1-cor(
        x[ids,],method=cor_method
      )), clustering_method)
      hc$height <- round(hc$height, 6)
      ensemble_hc[[i]]=cutree(hc,h=hs[i])
    }
  } else if (clustering_algorithm=="nj") {
    for (i in 1:n) {
      if (bootstrap==TRUE) p=1
      ids=sample(vargenes,p*length(vargenes),replace=bootstrap)
      hc=as.hclust(force.ultrametric(root(nj(as.dist(1-cor(
        x[ids,],method=cor_method
      ))),outgroup=1),method="extend"))
      ensemble_hc[[i]]=cutree(hc,h=hs[i])
    }
  }
  coocmat=matrix(0,nrow=ncol(x),ncol=ncol(x))
  colnames(coocmat)=colnames(x)
  rownames(coocmat)=colnames(x)
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(x)) {
      if (!i<j) {
        a=colnames(x)[i]; b=colnames(x)[j]
        cooc=sum(unlist(lapply(ensemble_hc, function(el) el[a]==el[b])))
        coocmat[i,j]=cooc
        coocmat[j,i]=cooc
      }
    }
  }
  if (clustering_algorithm=="hclust") {
    tree=as.phylo(hclust(dist(coocmat), method=clustering_method))
  } else if (clustering_algorithm=="nj") {
    tree=as.phylo(nj(dist(coocmat)))
  }
  return(list(
    tree=tree,
    cooccurrence=coocmat
  ))
}

#' Convert list of trees to multiPhylo object
#' 
listTreesAsMultiPhylo <- function(x) {
  mp = x[[1]]
  for (i in 2:length(x))
    mp = c(mp, x[[i]])
  return(mp)
}

#' Turn low supported nodes to politomies
#' 
#' @param tree object of class 'phylo'
#' @param support numeric, vector of support values
#' @param thrs numeric, support threshold value, nodes with support 
#'   lower than this value will be turned into politomies
#'   
#' @return object of class 'phylo'
#' 
di2multi_support <- function(tree, support, thrs) {
  
  # add support to tree
  support[is.na(support)] <- 0
  tree_tb <- as_tibble(tree)
  tree_tb$support <- c(rep(NA, length(tree$tip.label)), support)
  
  # set low support branches to 0
  topoly <- !is.na(tree_tb$support) & tree_tb$support<thrs
  tree_tb$branch.length[topoly] <- 0
  tree_new <- as.phylo(tree_tb)
  
  # turn short branches into politomies
  polytree <- di2multi(tree_new, tol=0.0)
  #plot(polytree)
  return(polytree)
  
}

#' Get gap genes for internal nodes in a tree
#' 
#' @param tree
#' @param feature_matrix, rows should be features (genes or motifs) and columns should be 
#'   tips of the tree; data should be log2 scaled
# @param feature_inclusion_ths numeric, threshold to consider a feature (gene, motif) (default: 1)
#' @param branch_length_ths numeric, threshold for selecting long branches (default: 0.1)
#' @param feature_in_thrs min threshold median value for feature in 
#'   the columns of feature_matrix in the selected node split
#' @param feature_out_thrs max threshold median value for features in 
#'   the columns of feature_matrix outside of selected node split
#' @param method,methodbg character, "absolute" or "median"
#' @param ncores integer (default: `detectCores()-1`)
#' @param verbose logical (default: FALSE)
#' 
tree_gap_genes <- function(
  tree, feature_matrix, branch_length_thrs=0.1,
  feature_in_thrs = 1.5, feature_out_thrs = 1,
  method="absolute", methodbg="absolute", abs_leakyness=0, abs_leakynessbg=0.05,
  ncores=detectCores()-1, verbose=FALSE
) {
  
  # get tree data
  tree_tb=as_tibble(tree)
  treet=tidytree::as.treedata(tree)
  
  # get the order of the tips
  tip_labels <- treet@phylo$tip.label
  is_tip <- treet@phylo$edge[,2] <= length(tip_labels)
  ordered_tips <- treet@phylo$edge[is_tip, 2]
  
  # get the ordered tip labels
  tip_labels <- treet@phylo$tip.label[ordered_tips]
  
  # select iternal nodes with large branch lengths. CRITICAL!!
  # (start at n_tips + 2 (because the first internal node is the whole tree!)
  nodes <- 1:nrow(tree_tb)
  nodes_above_ths <- tree_tb$branch.length[(length(tip_labels)+2):nrow(tree_tb)] > branch_length_thrs
  long_nodes <- which(nodes_above_ths)+length(tip_labels)+1
  
  calc_nodes=c(ordered_tips,long_nodes)
  
  # matrix to store the binary split info  1/0
  m_splits=matrix(0,ncol=length(tip_labels),nrow=length(calc_nodes))
  colnames(m_splits)=tip_labels
  rownames(m_splits)=calc_nodes #1:nrow(m_splits)
  # list to store the enriched feature in the tips under the node vs all other tips
  top_features=vector("list",length=length(nodes))
  # list to store the enriched features in the tips under the node vs sister clade (branching from mrca)
  top_feature_sister=vector("list",length=length(nodes))
  # list to store the anti-enriched features (or enriched oustide) the node
  top_feature_anti=vector("list",length=length(nodes))
  
  # require(foreach)
  # require(doParallel)
  # maxcores <- detectCores()
  # if (ncores > maxcores) {
  #   warning(sprintf(
  #     "Specified numer of cores (%s) exceeds the number of available cores (%s)!\nContinuing using %s cores.",
  #     ncores, maxcores, maxcores
  #   ))
  #   ncores=maxcores
  # }
  #registerDoParallel(ncores)
  #i=0
  #foreach (node=nodes) %dopar% {
  for(node in nodes){
    if (node %in% calc_nodes) {
      if (verbose==TRUE) print(sprintf("Starting analysis for node %s", node))
      #i=i+1
      i=match(node,calc_nodes)
      node_parent <- treeio::parent(tree,node)
      node_children <- treeio::child(tree,node_parent)
      node_sibling <- setdiff(node_children,node)
      if (node %in% ordered_tips) {
        tips_in <- tip_labels[node]
      } else {
        tips_in <- treeio::tree_subset(treet,node,levels_back=0)@phylo$tip.label
      }
      if (verbose==TRUE) print(sprintf("Tips in: %s", paste(tips_in,collapse=", ")))
      sister_label <- tree_tb[tree_tb$node%in%node_sibling,]$label
      tips_sister <- tryCatch({
        if (any(!is.na(sister_label))) {
          sister_label[!is.na(sister_label)]
        } else {
          unlist(lapply(node_sibling, function(node_sibling) 
            treeio::tree_subset(treet,node_sibling,levels_back=0)@phylo$tip.label
          ))
        }
      })#, error=function(e) NULL)
      if (verbose==TRUE) print(sprintf("Sister tips: %s", paste(tips_sister,collapse=", ")))
      tips_out <- setdiff(treet@phylo$tip.label,tips_in)
      if (verbose==TRUE) print(sprintf("Out tips: %s", paste(tips_out,collapse=", ")))
      
      # fill the reference split matrix
      m_splits[i,tips_in] <- 2
      m_splits[i,tips_out] <- 1
      .calc_leakyness <- function(abs_leakyness,n) {
        if (abs_leakyness<1) {
          pmin(round(c(1-abs_leakyness)*n), n)
        } else {
          pmin(round(n-abs_leakyness), n)
        }
      }
      if (method=="median") {
        fin=apply(feature_matrix[,tips_in,drop=F],1,median) > feature_in_thrs
        fin_inv=apply(feature_matrix[,tips_out,drop=F],1,median) > feature_in_thrs
      } else if (method=="absolute") {
        fin=apply(feature_matrix[,tips_in,drop=F], 1, function(x) 
          #all(x > feature_in_thrs)
          !(sum(x > feature_in_thrs) < .calc_leakyness(abs_leakyness,n=length(tips_in)))
        )
        fin_inv=apply(feature_matrix[,tips_out,drop=F], 1, function(x) 
          #all(x > feature_in_thrs)
          !(sum(x > feature_in_thrs) <  .calc_leakyness(abs_leakyness,n=length(tips_out)))
        )
      }
      if (methodbg=="median") {
        fout=apply(feature_matrix[,tips_out,drop=F],1,median) < feature_out_thrs
        fout_inv=apply(feature_matrix[,tips_in,drop=F],1,median) < feature_out_thrs
        fout_sister=apply(feature_matrix[,tips_sister,drop=F],1,median) < feature_out_thrs
      } else if (methodbg=="absolute") {
        fout=apply(feature_matrix[,tips_out,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_out)))
        )
        fout_inv=apply(feature_matrix[,tips_in,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_in)))
        )
        fout_sister=apply(feature_matrix[,tips_sister,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_sister)))
        )
      }
      
      f_in=which(fin & fout)
      f_sister=which(fin & fout_sister)
      f_out=which(fin_inv & fout_inv)
      
      # now add each feature to the corresponding list (or "none")
      if(length(f_in)>0){
        top_feat=rownames(feature_matrix)[f_in]
        if(length(f_in)>1)
          top_feat=top_feat[order(apply(feature_matrix[top_feat,tips_in,drop=F],1,median),decreasing=T)] 
        top_features[[node]]=top_feat
      } else { top_features[[node]]="NONE"}	
      
      if(length(f_sister)>0){
        top_feat=rownames(feature_matrix)[f_sister]
        if(length(f_sister)>1)
          top_feat=top_feat[order(apply(feature_matrix[top_feat,tips_in,drop=F],1,median),decreasing=T)] 
        top_feature_sister[[node]]=top_feat
      } else { top_feature_sister[[node]]="NONE"}	
      
      if(length(f_out)>0){
        top_anti=rownames(feature_matrix)[f_out]
        if(length(f_out)>1)
          top_anti=top_anti[order(apply(feature_matrix[top_anti,tips_out,drop=F],1,median),decreasing=T)]
        top_feature_anti[[node]]=top_anti
      } else { top_feature_anti[[node]]="NONE"}
      
      if (verbose==TRUE) print(sprintf("Finished analysis for node %s", node))
      
    } else {
      if (verbose==TRUE) 
        print(sprintf(
          "Node %s is shorther than threshold branch length (%s); skipping.", 
          node, branch_length_thrs
        ))
      top_features[[node]] <- "NOT CALCULATED"
      top_feature_sister[[node]] <- "NOT CALCULATED"
      top_feature_anti[[node]] <- "NOT CALCULATED"
    }
  }
  #stopImplicitCluster()
  hm <- Heatmap(
    m_splits, col=c("0"="gray88","1"="lightgrey","2"="black"), show_heatmap_legend=FALSE,
    border = TRUE, rect_gp = gpar(col="gray80"),
    cluster_columns = FALSE, cluster_rows = FALSE
  )
  
  # top_features_v <- vector("list",nrow(tree_tb))
  # top_feature_anti_v <- vector("list",nrow(tree_tb))
  # for (i in 1:length(calc_nodes)) {
  #   ln <- calc_nodes[i]
  #   top_features_v[[ln]] <- top_features[[i]]
  #   top_feature_anti_v[[ln]] <- top_feature_anti[[i]]
  # }
  
  return(list(
    top_features=top_features, 
    top_feature_anti=top_feature_anti, 
    top_feature_sister=top_feature_sister,
    splits_mat=m_splits, heatmap=hm
  ))
}


#' Calculate support on nodes and get tree with politomies based on support
#' 
saveTreeFiles <- function(
  tree_rds, ctcol, out_dir, save_nexus=FALSE, scale_branches=FALSE, 
  width_tree=10, height_tree=14, width_polytree=8, height_polytree=12
) {
  method <- str_extract(tree_rds,"nj|hclust(?=\\.RDS)"); if(is.na(method)) method="hclust"
  tree <- readRDS(tree_rds)
  treename <- stringr::str_remove(basename(tree_rds),"\\.RDS")
  iter <- 1000
  support_trees <- vector(mode="list",length=iter)
  p <- 0.5
  for(i in 1:iter){ support_trees[[i]]=treeDownsampling(mc_fp,vargenes=var_genes,p=p,method=method) }
  support=prop.clades(tree, support_trees)
  tip_color=ctcol[tree$tip.label]
  if (scale_branches==TRUE)
    tree <- compute.brlen(tree, 1)
  pdf(sprintf("%s/%s.pdf",outdir,treename),height=height_tree,width=width_tree,useDingbats=TRUE)
  plot(tree,tip.color=tip_color,cex=1.5)
  drawSupportOnEdges(support,col="black",bg="white",frame="circle")
  dev.off()
  if (save_nexus) {
    col_tib <- tibble(color=gplots::col2hex(tip_color),label=names(tip_color))
    tree_tib <- as_tibble(tree)
    tree_tib$support <- c(rep(NA, length(tree$tip.label)), support)
    tree_full_tib <- treeio::full_join(tree_tib,col_tib, by="label")
    col_labs <- unlist(lapply(1:nrow(tree_full_tib), function(i) {
      lb=tree_full_tib$label[i]
      cl=tree_full_tib$color[i]
      if (!is.na(cl)) {
        sprintf("%s[&!color=%s]",lb,cl)
      } else {
        lb
      }
    }))
    tree_full_tib$label <- col_labs
    tree_full_dat <- as.treedata(tree_full_tib)
    write.beast(tree_full_dat, file=file.path(sprintf("%s/%s.beast",outdir,treename)))
  }
  polit <- 0.1
  supp <- polit*iter
  polytree <- di2multi_support(tree, support, supp)
  tree_tib <- as_tibble(polytree)
  missing_terminal_nodes <- tree_tib$branch.length<0.01 & !is.na(tree_tib$label)
  tree_tib$branch.length[missing_terminal_nodes] <- tree_tib$branch.length[missing_terminal_nodes]+1
  polytree <- as.phylo(tree_tib)
  if (scale_branches==TRUE)
    polytree <- compute.brlen(polytree, 1)
  pdf(
    sprintf("%s/%s_polytomies_support%.2f.pdf",outdir,treename,polit),
    height=height_polytree,width=width_polytree,useDingbats=TRUE
  )
  par(mar=rep(2, 4))
  plot(polytree, main=sprintf("Polytomies: support < %.0f%%",polit*100), cex=1.5, tip.color=tip_color)
  dev.off()
  if (save_nexus) {
    col_tib <- tibble(color=gplots::col2hex(tip_color),label=names(tip_color))
    tree_col_tib <- treeio::full_join(polytree,col_tib,by="label")
    tree_tib <- as_tibble(polytree)
    tree_tib_support <- c(rep(NA, length(polytree$tip.label)), support)
    if (length(tree_tib_support)-nrow(tree_tib)==1) {
      tree_tib_support <- tree_tib_support[2:length(tree_tib_support)]
    }
    tree_tib$support <- tree_tib_support
    tree_full_tib <- treeio::full_join(tree_tib,col_tib, by="label")
    col_labs <- unlist(lapply(1:nrow(tree_full_tib), function(i) {
      lb=tree_full_tib$label[i]
      cl=tree_full_tib$color[i]
      if (!is.na(cl)) {
        sprintf("%s[&!color=%s]",lb,cl)
      } else {
        lb
      }
    }))
    tree_full_tib$label <- col_labs
    tree_full_dat <- as.treedata(tree_full_tib)
    write.beast(tree_full_dat,file=file.path(sprintf("%s/%s_politomies%.2f.beast",outdir,treename,polit)))
  }
}

#' Plot correlation/distance/co-ocurrence matrix based on which the tree is built
#' 
#' @param build_tree logical, should the tree be built using `ggtree`, and 
#'   order infered from this object? If FALSE, order of tips in `tree$tip.label` 
#'   is used (default: TRUE)
#'   
plotTreeMatrix <- function(
  mat, tree, build_tree=TRUE, rev=TRUE,
  ctcol, hmapcol=NULL, name, 
  out_name, width=7, height=5,
  show_row_names = TRUE, row_names_side = "left", 
  show_column_names = TRUE, column_names_side="bottom",
  annotation_side = c("right", "top")
) {
  if (is.null(hmapcol))
    hmapcol=colorRamp2(
      c(0,100,200,400,600,700,800,900,1000),
      colors=c(c("white",'#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#990000'))
    )
  if (build_tree) {
    gt <- ggplot_build(ggtree(tree)+theme_tree2()+geom_tiplab()); gtd <- gt$data[[3]]; setDT(gtd); setorder(gtd,"y")
    labels_ord <- gtd$label
  } else {
    labels_ord <- tree$tip.label
  }
  if (rev) labels_ord <- rev(labels_ord)
  row_ha = rowAnnotation(
    ct = labels_ord, col= list(ct = ctcol[labels_ord]),
    border = TRUE, show_legend = FALSE, show_annotation_name = FALSE
  )
  col_ha = columnAnnotation(
    ct = labels_ord, col= list(ct = ctcol[labels_ord]), 
    border = TRUE, show_legend = FALSE, show_annotation_name = FALSE
  )
  bottom_annotation <- NULL
  top_annotation <- NULL
  right_annotation <- NULL
  left_annotation <- NULL
  if(length(annotation_side)>0) {
    if (any(grepl("top",annotation_side))) top_annotation <- col_ha
    if (any(grepl("bottom",annotation_side))) bottom_annotation <- col_ha
    if (any(grepl("right",annotation_side))) right_annotation <- row_ha
    if (any(grepl("left",annotation_side))) left_annotation <- row_ha
  }
  hm <- ComplexHeatmap::Heatmap(
    mat[labels_ord,labels_ord],
    name = name,
    col = hmapcol,
    border = TRUE, rect_gp =  gpar(col = "gray66", lwd = 0.2),
    cluster_rows = FALSE, cluster_columns = FALSE, 
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_gp = gpar(col = ctcol[labels_ord]),
    column_names_gp = gpar(col = ctcol[labels_ord]),
    bottom_annotation = bottom_annotation, top_annotation = top_annotation,
    right_annotation = right_annotation, left_annotation = left_annotation,
    show_row_names=show_row_names, row_names_side = row_names_side,
    show_column_names = show_column_names, column_names_side=column_names_side
  )
  pdf(out_name,width=width,height=height)
  print(hm)
  dev.off()
}

#' Add colors to tips of the tree for display in FigTree
#' 
#' @param tree_file tree file in nexus format
#' @param ctcol character, vector of colors for tips, names should be `tip.label` of the `tree`
#' @value `phylo` object which is also saved to a file same as input nexus file, but with beast extension
#' 
treeColorTips <- function(tree=NULL, tree_file=NULL, out_file=NULL, ctcol) {
  if (is.null(tree)) {
    if (is.null(tree_file)) 
      stop("If no tree object is passed to function, a tree file must be specified!")
    message("Loading tree from ", tree_file)
    tree=as.phylo(read.beast(tree_file))
    if (is.null(out_file))
      out_file <- str_replace(tree_file,"\\.nexus","\\.beast")
    if (all.equal(tree_file,out_file)==TRUE)
      out_file <- str_replace(tree_file,"\\.beast","\\.2.beast")
  }
  tip_color=ctcol[tree$tip.label]
  col_tib <- tibble(color=gplots::col2hex(tip_color),label=names(tip_color))
  tree_col_tib <- full_join(tree,col_tib,by="label")
  tree_tib <- as_tibble(tree)
  #tree_tib_support <- c(rep(NA, length(tree$tip.label)), support)
  tree_full_tib <- full_join(tree_tib,col_tib, by="label")
  col_labs <- unlist(lapply(1:nrow(tree_full_tib), function(i) {
    lb=tree_full_tib$label[i]
    cl=tree_full_tib$color[i]
    if (!is.na(cl)) {
      sprintf("%s[&!color=%s]",lb,cl)
    } else {
      lb
    }
  }))
  tree_full_tib$label <- col_labs
  tree_full_dat <- as.treedata(tree_full_tib)
  
  if (is.null(out_file))
    out_file <- "tree.beast"
  write.beast(tree_full_dat, file=out_file)
  message("Saved the tree to ", out_file)
  return(as.phylo(tree_full_dat))
}
