require(zoo)
require(scales)
require(RColorBrewer)
require(viridis)
require(data.table)
require(stringr)
require(circlize)
require(ComplexHeatmap)
require(LaplacesDemon)
require(preprocessCore)
library(ggpubr)


# # # # # # # # # # # # # # # # #
#                               #
#       UTILITY FUNCTIONS       #
#                               #
# # # # # # # # # # # # # # # # #

#' Overalap between all cells in two matrices
overlap <- function(M1, M2) {
  # make matrices conformable
  missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
  missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
  add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
  add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
  rownames(add_M2) <- missing_in_M2
  rownames(add_M1) <- missing_in_M1
  M1 <- data.matrix(rbind(M1, add_M1))
  M2 <- data.matrix(rbind(M2, add_M2))
  # intersect: matrix multiplication t(M1) x M2
  outM <- t(M1) %*% M2
  rownames(outM) <- colnames(M1)
  colnames(outM) <- colnames(M2)
  outM
}

#' Intersect normalized by smaller group size
intersect_by_min <- function(M1, M2) {
  # make matrices conformable
  missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
  missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
  add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
  add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
  rownames(add_M2) <- missing_in_M2
  rownames(add_M1) <- missing_in_M1
  M1 <- data.matrix(rbind(M1, add_M1))
  M2 <- data.matrix(rbind(M2, add_M2))
  # intersect: matrix multiplication t(M1) x M2
  intM <- t(M1) %*% M2
  rownames(intM) <- colnames(M1)
  colnames(intM) <- colnames(M2)
  
  modlen1 <- colSums(M1)
  modlen2 <- colSums(M2)
  minM <- intM
  for(i in 1:nrow(minM)) {
    for(j in 1:ncol(minM)) {
      minM[i,j] <- min(modlen1[i], modlen2[j])
    }
  }
  
  outM <- intM/minM
  return(outM)
}

#' Jaccard distance between columns of two binary matrices
#' @param M1,M2 binary matrix
jaccard <- function(M1, M2) {
  # make matrices conformable
  missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
  missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
  add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
  add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
  rownames(add_M2) <- missing_in_M2
  rownames(add_M1) <- missing_in_M1
  M1 <- data.matrix(rbind(M1, add_M1))
  M2 <- data.matrix(rbind(M2, add_M2))
  # intersect: matrix multiplication t(M1) x M2
  intersectM <- t(M1) %*% M2
  # union: sum by rows and columns
  unionM <- matrix(nrow = ncol(M1), ncol = ncol(M2))
  for(i in 1:ncol(M1)) {
    unionM[i,] <- colSums((M1[,i] + M2) > 0)
  }
  # jaccard
  outM <- intersectM/unionM
  outM[is.na(outM)] <- 0
  rownames(outM) <- colnames(M1)
  colnames(outM) <- colnames(M2)
  outM
}

#' Jensen-Shannon Divergence
#library(philentropy)
#library(gtools)
#mc_umifrac <- apply(mc_umifrac, 2, function(x) x/1000)
#jsd <- JSD(as.matrix(t(mc_umifrac[var_genes,])),unit="log2")
#colnames(jsd) <- rownames(jsd) <- colnames(mc_umifrac)
#jsd_dist <- as.matrix(1-sqrt(jsd))

#' KLD between columns (metacells, cell types) of matrix (footprint, genes in rows)
calcKLD=function(mc_fp) {
  nc=ncol(mc_fp)
  kld_mat=matrix(NA,nrow=nc,ncol=nc)
  kld_list=vector("list",nc)
  for (i in 1:nc) {
    message(sprintf("%s of %s",i,nc))
    for (j in 1:nc) {
      if (!i>j) {
        kld=LaplacesDemon::KLD(px=mc_fp[,i],py=mc_fp[,j])
        kld_mat[i,j]=kld$sum.KLD.px.py
        kld_mat[j,i]=kld$sum.KLD.py.px
        kld_list[[i]][[j]]=kld$KLD.px.py
        kld_list[[j]][[i]]=kld$KLD.py.px
      }
    }
  }
  colnames(kld_mat)=colnames(mc_fp);rownames(kld_mat)=colnames(mc_fp)
  list(mat=kld_mat,probs=kld_list)
}

#' Quantile normalization
#' @param x matrix or data.frame
quantile_normalisation <- function(x){
  df_rank <- data.frame(apply(x,2,rank,ties.method="min"))
  df_sorted <- data.frame(apply(x, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(x)
  return(df_final)
}

# # # # # # # # # # # # # # # # #
#                               #
#     CSPS INIT FUNCTIONS       #
# create csps object as output  #
#                               #
# # # # # # # # # # # # # # # # #

#' Create cross-species comparison object.
#' 
#' @param sp1_fp_fn,sp2_fp_fn path to expression matrices for the analyzed species, 
#'   rows are genes and columns are cells, metacells or cell types (usually these 
#'   are metacell footprint matrices, `mc@mc_fp`); this can be either raw fold changes 
#'   UMI counts or UMI fraction matrices in either RDS or tab-separated text format
#' @param OG_pairs_fn path to file with broccoli-style pairwise orthologies in
#'   tab-separated text format; make sure they are in the right order
#' @param sp_names character of length 2, optional species names/abbreviations;
#'   if not specified, will use `c("sp1","sp2")`
#' @param make_sp_colnames logical, append the species abbreviations at the 
#'   begining of column names (default: TRUE)
#' @param quant_norm logical (default: TRUE)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` 
#'   we require in each species (default: 1, set to NULL to skip filtering by fc)
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL 
#'   to skip filtering by fc)
#' @param one2one logical (default: TRUE)
#'
#' @return a list with the following elements:
#'   * merged, a matrix with selected orthologous genes in rows and metacells 
#'     (or cell types) from both species in columns
#'   * og_pairs, data.frame with orthologous pairs in columns
#'   * sp1, a matrix with selected orthologous genes in rows and metacells 
#'     of first species in columns
#'   * sp2, a matrix with selected orthologous genes in rows and metacells 
#'     of the second species in columns
#'   * top_cross_sp1, top orthologous genes in the first species
#'   * top_cross_sp2, top orthologous genes in the second species
#' This will be a standard object for cross-species analyses
#' 
csps_create_crossspecies_object <- function(
  sp1_fp_fn, sp2_fp_fn, OG_pairs_fn, sp_names=c("sp1","sp2"), 
  make_sp_colnames=TRUE, quant_norm=TRUE, cross_n=1, cross_fc_thrs=2, one2one=TRUE
){ 
  
  # Read and parse data  
  nc=nchar(sp1_fp_fn)
  if (substr(sp1_fp_fn,nc-3,nc)==".RDS") {
    sp1_fp=readRDS(sp1_fp_fn)
  } else {
    sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  nc=nchar(sp2_fp_fn)
  if (substr(sp2_fp_fn,nc-3,nc)==".RDS") {
    sp2_fp=readRDS(sp2_fp_fn)
  } else {
    sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  OG_pairs=fread(OG_pairs_fn,header=FALSE,sep="\t",stringsAsFactors=FALSE)
  
  setnames(OG_pairs, c(sp_names))
  sp1=sp_names[1]
  sp2=sp_names[2]
  sp1_fp=sp1_fp[intersect(rownames(sp1_fp),OG_pairs[[sp1]]),]
  sp2_fp=sp2_fp[intersect(rownames(sp2_fp),OG_pairs[[sp2]]),]
  if (make_sp_colnames) {
    colnames(sp1_fp) <- paste0(sp1,colnames(sp1_fp))
    colnames(sp2_fp) <- paste0(sp2,colnames(sp2_fp))
  }
  
  # Reduce orthology table to only expressed genes
  OG_pairs=OG_pairs[which(OG_pairs[[sp1]] %in% rownames(sp1_fp) & OG_pairs[[sp2]] %in% rownames(sp2_fp)),]
  
  # Eliminate non-one2one orthologs (relaxed non-expression criterion)  
  OG_pairs_one2one=as.data.frame(
    OG_pairs[!duplicated(OG_pairs[[sp1]]) & !duplicated(OG_pairs[[sp2]]),]
  )
  
  if (one2one==TRUE) {
    
    OG_pairs = OG_pairs_one2one
    cross1 <- OG_pairs[[sp1]]
    
  } else { 
    # Allow for one2many (2-3 max) relationships, DUPLICATING entries for the paralogs.
    # one2one
    one2one <- OG_pairs_one2one[[sp1]]
    # one2many (2 or three)
    OG_pairs[[sp1]] <- make.names(OG_pairs[[sp1]],unique=TRUE)
    one2many <- grep("\\.[12]$",OG_pairs[[sp1]],value=TRUE)
    one2many <- unique(sort(c(one2many, str_remove(one2many,"\\.[12]$"))))
    
    OG_pairs = OG_pairs[OG_pairs[[sp1]] %in% c(one2one,one2many),]
    cross1 <- str_remove(OG_pairs[[sp1]],"\\.[12]$")
  }
  
  # Reorder matrices, quantile_norm? and merge matrices
  if(class(sp1_fp)=="dgCMatrix"|class(sp2_fp)=="dgCMatrix") {
    merged=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[OG_pairs[[sp2]],]),sparse=TRUE)
  } else {
    merged=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[OG_pairs[[sp2]],]))
  }
  rnm=OG_pairs[[sp1]]
  rownames(merged)=rnm
  cnm=colnames(merged)
  
  if(quant_norm){
    if(any(class(merged)=="dgCMatrix")) {
      sparse=TRUE
      merged=as.matrix(merged)
      copymat=FALSE
    } else {
      sparse=FALSE
      copymat=TRUE
    }
    merged=tryCatch({
      preprocessCore::normalize.quantiles(merged,copy=copymat)
    }, error = function(e){
      warning(e)
      quantile_normalisation(merged)
    })
    if (sparse==TRUE)
      merged=Matrix::Matrix(merged,sparse=TRUE)
    rownames(merged)=rnm
    colnames(merged)=cnm
  }
  
  # Select variable genes in BOTH species
  if (!is.null(cross_fc_thrs)&!is.null(cross_n)) {
    top_cross=names(which(
      apply(merged[,1:ncol(sp1_fp)],1,function(x) 
        sort(x,decreasing=T)[cross_n]) > cross_fc_thrs & 
        apply(merged[,(ncol(sp1_fp)+1):ncol(merged)],1,function(x) 
          sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    ))
  } else {
    top_cross=rownames(merged)
  }
  
  # Return a list with both matrices
  out_m1=merged[,1:ncol(sp1_fp)]
  rownames(out_m1)=OG_pairs[[sp1]]
  out_m2=merged[,(ncol(sp1_fp)+1):ncol(merged)]
  rownames(out_m2)=make.names(OG_pairs[[sp2]],unique=TRUE)
  
  ids <- unlist(lapply(top_cross,function(x) which(OG_pairs[[sp1]]==x)))
  top_cross_sp2=OG_pairs[ids,][[sp2]]
  #top_cross <- str_remove(top_cross,"\\.\\d+")
  #top_cross_sp2 <- str_remove(top_cross2,"\\.\\d+")
  
  csps <- list(
    merged=merged,
    og_pairs=OG_pairs,
    sp1=out_m1, sp2=out_m2,
    top_cross_sp1=top_cross, top_cross_sp2=top_cross_sp2
  )
  return(csps)
}

#' @inheritParams csps_create_crossspecies_object
#' @seealso [csps_create_crossspecies_object()]
csps_create_3way_crossspecies_object <- function(
  sp1_fp_fn, sp2_fp_fn, sp3_fp_fn, 
  OG_pairs1_2_fn, OG_pairs1_3_fn, OG_pairs2_3_fn, 
  sp_names=c("sp1","sp2","sp3"), make_sp_colnames=TRUE,
  one2one=TRUE, quant_norm=TRUE, cross_n=1, cross_fc_thrs=2
){
  # Read and parse data  
  if (length(sp_names)!=3)
    stop("Length of sp_names should be 3!")
  nc=nchar(sp1_fp_fn)
  if (substr(sp1_fp_fn,nc-3,nc)==".RDS") {
    sp1_fp=readRDS(sp1_fp_fn)
  } else {
    sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  nc=nchar(sp2_fp_fn)
  if (substr(sp2_fp_fn,nc-3,nc)==".RDS") {
    sp2_fp=readRDS(sp2_fp_fn)
  } else {
    sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  if (substr(sp3_fp_fn,nc-3,nc)==".RDS") {
    sp3_fp=readRDS(sp3_fp_fn)
  } else {
    sp3_fp=read.table(sp3_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  OG_pairs1_2=fread(OG_pairs1_2_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs1_3=fread(OG_pairs1_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs2_3=fread(OG_pairs2_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  sp1=sp_names[1]
  sp2=sp_names[2]
  sp3=sp_names[3]
  setnames(OG_pairs1_2, c(sp1,sp2))
  setnames(OG_pairs1_3, c(sp1,sp3))
  setnames(OG_pairs2_3, c(sp2,sp3))
  intersect1=intersect(intersect(rownames(sp1_fp),OG_pairs1_2[[sp1]]),OG_pairs1_3[[sp1]])
  intersect2=intersect(intersect(rownames(sp2_fp),OG_pairs1_2[[sp2]]),OG_pairs2_3[[sp2]])
  intersect3=intersect(intersect(rownames(sp3_fp),OG_pairs1_3[[sp3]]),OG_pairs2_3[[sp3]])
  sp1_fp=sp1_fp[intersect1,]
  sp2_fp=sp2_fp[intersect2,]
  sp3_fp=sp3_fp[intersect3,]
  if (make_sp_colnames) {
    colnames(sp1_fp) <- paste0(sp1,colnames(sp1_fp))
    colnames(sp2_fp) <- paste0(sp2,colnames(sp2_fp))
    colnames(sp3_fp) <- paste0(sp3,colnames(sp3_fp))
  }
  
  # Reduce orthology table to only expressed genes
  OG_pairs1_2=OG_pairs1_2[which(OG_pairs1_2[[sp1]] %in% rownames(sp1_fp) & 
                                  OG_pairs1_2[[sp2]] %in% rownames(sp2_fp)),]
  OG_pairs1_3=OG_pairs1_3[which(OG_pairs1_3[[sp1]] %in% rownames(sp1_fp) & 
                                  OG_pairs1_3[[sp3]] %in% rownames(sp3_fp)),]
  OG_pairs2_3=OG_pairs2_3[which(OG_pairs2_3[[sp2]] %in% rownames(sp2_fp) & 
                                  OG_pairs2_3[[sp3]] %in% rownames(sp3_fp)),]
  
  # Eliminate non-one2one orthologs (relaxed non-expression criterion)  
  OG_pairs1_2_one2one=as.data.frame(
    OG_pairs1_2[!duplicated(OG_pairs1_2[[sp1]]) & !duplicated(OG_pairs1_2[[sp2]]),]
  )
  OG_pairs1_3_one2one=as.data.frame(
    OG_pairs1_3[!duplicated(OG_pairs1_3[[sp1]]) & !duplicated(OG_pairs1_3[[sp3]]),]
  )
  OG_pairs2_3_one2one=as.data.frame(
    OG_pairs2_3[!duplicated(OG_pairs2_3[[sp2]]) & !duplicated(OG_pairs2_3[[sp3]]),]
  )
  
  if (one2one==TRUE) {
    
    OG_pairs1_2 = OG_pairs1_2_one2one
    OG_pairs1_3 = OG_pairs1_3_one2one
    OG_pairs2_3 = OG_pairs2_3_one2one
    mdt1=merge.data.table(OG_pairs1_2,OG_pairs1_3,by=sp1)
    mdt=merge.data.table(mdt1,OG_pairs2_3,by=c(sp2,sp3))
    cross1=mdt[[sp1]]
    cross2=mdt[[sp2]]
    cross3=mdt[[sp3]]
    
  } else { 
    
    stop("one2many not implemented yet, set one2one=TRUE")
    
  }
  
  # Reorder matrices, quantile_norm and merge matrices
  if (any(c(class(sp1_fp),class(sp2_fp),class(sp3_fp))=="dgCMatrix")) {
    merged1=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]),sparse=TRUE)
    merged=Matrix::Matrix(cbind(merged1,sp3_fp[cross3,]),sparse=TRUE)
  } else {
    merged1=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]))
    merged=data.matrix(cbind(merged1,sp3_fp[cross3,]))
  }
  rnm=cross1
  rownames(merged)=rnm
  cnm=colnames(merged)
  
  if(quant_norm){
    if(any(class(merged)=="dgCMatrix")) {
      sparse=TRUE
      merged=as.matrix(merged)
      copymat=FALSE
    } else {
      sparse=FALSE
      copymat=TRUE
    }
    merged=tryCatch({
      preprocessCore::normalize.quantiles(merged,copy=copymat)
    }, error = function(e){
      quantile_normalisation(merged)
    })
    if (sparse==TRUE)
      merged=Matrix::Matrix(merged,sparse=TRUE)
    rownames(merged)=rnm
    colnames(merged)=cnm
  }
  
  # Select variable genes in ALL species
  nc1=ncol(sp1_fp)
  nc2=ncol(sp1_fp)+ncol(sp2_fp)
  nc3=ncol(merged)
  if (!is.null(cross_fc_thrs)&!is.null(cross_n)) {
    var1= apply(merged[,1:nc1],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    var2= apply(merged[,c(nc1+1):nc2],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    var3= apply(merged[,c(nc2+1):nc3],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    top_cross=names(which(var1 & var2 & var3))
  } else {
    top_cross=rownames(merged)
  }
  
  # Return a list with both matrices
  out_m1=merged[,1:nc1]
  rownames(out_m1)=cross1
  out_m2=merged[,(nc1+1):nc2]
  rownames(out_m2)=make.names(cross2,unique=TRUE)
  out_m3=merged[,(nc2+1):nc3]
  rownames(out_m3)=make.names(cross3,unique=TRUE)
  
  ids2 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_2[[sp1]]==x)))
  top_cross_sp2=OG_pairs1_2[ids2,][[sp2]]
  ids3 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_3[[sp1]]==x)))
  top_cross_sp3=OG_pairs1_3[ids3,][[sp3]]
  
  csps <- list(
    merged=merged,
    OG_pairs1_2=OG_pairs1_2, OG_pairs1_3=OG_pairs1_3, OG_pairs2_3=OG_pairs2_3,
    sp1=out_m1, sp2=out_m2, sp3=out_m3,
    top_cross_sp1=top_cross, top_cross_sp2=top_cross_sp2, top_cross_sp2=top_cross_sp3
  )
  return(csps)
  
}

#' @inheritParams csps_create_crossspecies_object
#' @param restrict_paralogs logical, restrict number of possible paralogs for a gene 
#'   (default: FALSE, if TRUE, genes with more than 4 paralogs will not be considered)
#' @seealso [csps_create_crossspecies_object()]
csps_create_4way_crossspecies_object <- function(
  sp1_fp_fn, sp2_fp_fn, sp3_fp_fn, sp4_fp_fn,  
  OG_pairs1_2_fn, OG_pairs1_3_fn,OG_pairs1_4_fn, OG_pairs2_3_fn, OG_pairs2_4_fn, OG_pairs3_4_fn,
  sp_names=c("sp1","sp2","sp3","sp4"), make_sp_colnames=TRUE,
  one2one=TRUE, restrict_paralogs=FALSE, quant_norm=TRUE, 
  cross_n=1, cross_fc_thrs=2
){
  # Read and parse data  
  if (length(sp_names)!=4)
    stop("Length of sp_names should be 4!")
  nc=nchar(sp1_fp_fn)
  if (substr(sp1_fp_fn,nc-3,nc)==".RDS") {
    sp1_fp=readRDS(sp1_fp_fn)
  } else {
    sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  nc=nchar(sp2_fp_fn)
  if (substr(sp2_fp_fn,nc-3,nc)==".RDS") {
    sp2_fp=readRDS(sp2_fp_fn)
  } else {
    sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  if (substr(sp3_fp_fn,nc-3,nc)==".RDS") {
    sp3_fp=readRDS(sp3_fp_fn)
  } else {
    sp3_fp=read.table(sp3_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  if (substr(sp4_fp_fn,nc-3,nc)==".RDS") {
    sp4_fp=readRDS(sp4_fp_fn)
  } else {
    sp4_fp=read.table(sp4_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
  }
  OG_pairs1_2=fread(OG_pairs1_2_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs1_3=fread(OG_pairs1_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs1_4=fread(OG_pairs1_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs2_3=fread(OG_pairs2_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs2_4=fread(OG_pairs2_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  OG_pairs3_4=fread(OG_pairs3_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  sp1=sp_names[1]
  sp2=sp_names[2]
  sp3=sp_names[3]
  sp4=sp_names[4]
  setnames(OG_pairs1_2, c(sp1,sp2))
  setnames(OG_pairs1_3, c(sp1,sp3))
  setnames(OG_pairs1_4, c(sp1,sp4))
  setnames(OG_pairs2_3, c(sp2,sp3))
  setnames(OG_pairs2_4, c(sp2,sp4))
  setnames(OG_pairs3_4, c(sp3,sp4))
  intersect1=intersect(intersect(intersect(rownames(sp1_fp),OG_pairs1_2[[sp1]]),OG_pairs1_3[[sp1]]),OG_pairs1_4[[sp1]])
  intersect2=intersect(intersect(intersect(rownames(sp2_fp),OG_pairs1_2[[sp2]]),OG_pairs2_3[[sp2]]),OG_pairs2_4[[sp2]])
  intersect3=intersect(intersect(intersect(rownames(sp3_fp),OG_pairs1_3[[sp3]]),OG_pairs2_3[[sp3]]),OG_pairs3_4[[sp3]])
  intersect4=intersect(intersect(intersect(rownames(sp4_fp),OG_pairs1_4[[sp4]]),OG_pairs2_4[[sp4]]),OG_pairs3_4[[sp4]])
  sp1_fp=sp1_fp[intersect1,]
  sp2_fp=sp2_fp[intersect2,]
  sp3_fp=sp3_fp[intersect3,]
  sp4_fp=sp4_fp[intersect4,]
  if (make_sp_colnames) {
    colnames(sp1_fp) <- paste0(sp1,colnames(sp1_fp))
    colnames(sp2_fp) <- paste0(sp2,colnames(sp2_fp))
    colnames(sp3_fp) <- paste0(sp3,colnames(sp3_fp))
    colnames(sp4_fp) <- paste0(sp4,colnames(sp4_fp))
  }
  
  # Reduce orthology table to only expressed genes
  OG_pairs1_2=OG_pairs1_2[which(OG_pairs1_2[[sp1]] %in% rownames(sp1_fp) & 
                                  OG_pairs1_2[[sp2]] %in% rownames(sp2_fp)),]
  OG_pairs1_3=OG_pairs1_3[which(OG_pairs1_3[[sp1]] %in% rownames(sp1_fp) & 
                                  OG_pairs1_3[[sp3]] %in% rownames(sp3_fp)),]
  OG_pairs1_4=OG_pairs1_4[which(OG_pairs1_4[[sp1]] %in% rownames(sp1_fp) & 
                                  OG_pairs1_4[[sp4]] %in% rownames(sp4_fp)),]
  OG_pairs2_3=OG_pairs2_3[which(OG_pairs2_3[[sp2]] %in% rownames(sp2_fp) & 
                                  OG_pairs2_3[[sp3]] %in% rownames(sp3_fp)),]
  OG_pairs2_4=OG_pairs2_4[which(OG_pairs2_4[[sp2]] %in% rownames(sp2_fp) & 
                                  OG_pairs2_4[[sp4]] %in% rownames(sp4_fp)),]
  OG_pairs3_4=OG_pairs3_4[which(OG_pairs3_4[[sp3]] %in% rownames(sp3_fp) & 
                                  OG_pairs3_4[[sp4]] %in% rownames(sp4_fp)),]
  
  # Eliminate non-one2one orthologs (relaxed non-expression criterion)  
  OG_pairs1_2_one2one=as.data.frame(
    OG_pairs1_2[!duplicated(OG_pairs1_2[[sp1]]) & !duplicated(OG_pairs1_2[[sp2]]),]
  )
  OG_pairs1_3_one2one=as.data.frame(
    OG_pairs1_3[!duplicated(OG_pairs1_3[[sp1]]) & !duplicated(OG_pairs1_3[[sp3]]),]
  )
  OG_pairs1_4_one2one=as.data.frame(
    OG_pairs1_4[!duplicated(OG_pairs1_4[[sp1]]) & !duplicated(OG_pairs1_4[[sp4]]),]
  )
  OG_pairs2_3_one2one=as.data.frame(
    OG_pairs2_3[!duplicated(OG_pairs2_3[[sp2]]) & !duplicated(OG_pairs2_3[[sp3]]),]
  )
  OG_pairs2_4_one2one=as.data.frame(
    OG_pairs2_4[!duplicated(OG_pairs2_4[[sp2]]) & !duplicated(OG_pairs2_4[[sp4]]),]
  )
  OG_pairs3_4_one2one=as.data.frame(
    OG_pairs3_4[!duplicated(OG_pairs3_4[[sp3]]) & !duplicated(OG_pairs3_4[[sp4]]),]
  )
  
  if (one2one==TRUE) {
    
    OG_pairs1_2 = OG_pairs1_2_one2one
    OG_pairs1_3 = OG_pairs1_3_one2one
    OG_pairs1_4 = OG_pairs1_4_one2one
    OG_pairs2_3 = OG_pairs2_3_one2one
    OG_pairs2_4 = OG_pairs2_4_one2one
    OG_pairs3_4 = OG_pairs3_4_one2one
    mdt1=merge.data.table(OG_pairs1_2,OG_pairs1_3,by=sp1)
    mdt2=merge.data.table(mdt1,OG_pairs1_4,by=sp1)
    mdt3=merge.data.table(mdt2,OG_pairs2_3,by=c(sp2,sp3))
    mdt4=merge.data.table(mdt3,OG_pairs2_4,by=c(sp2,sp4))
    mdt=merge.data.table(mdt4,OG_pairs3_4,by=c(sp3,sp4))
    cross1=mdt[[sp1]]
    cross2=mdt[[sp2]]
    cross3=mdt[[sp3]]
    cross4=mdt[[sp4]]
    
  } else { 
    
    #stop("one2many not implemented yet, set one2one=TRUE")
    
    # Allow for one2many (3 max) relationships, DUPLICATING entries for the paralogs.
    multiOGs <- function(OG_pairs) {
      # one2one
      one2one <- OG_pairs[[1]]
      # one2many (2 or three)
      OG_pairs[[1]] <- make.names(OG_pairs[[1]],unique=TRUE)
      one2many <- grep("\\.[12]$",OG_pairs[[1]],value=TRUE)
      one2many <- unique(sort(c(one2many, str_remove(one2many,"\\.[12]$"))))
      
      OG_pairs = OG_pairs[OG_pairs[[1]] %in% c(one2one,one2many),]
      cross <- str_remove(OG_pairs[[1]],"\\.[12]$")
      #list(OG_pairs, cross)
      OG_pairs[[1]] <- cross
      OG_pairs
    }
    
    ogs_dts <- list(
      OG_pairs1_2_one2one, OG_pairs1_3_one2one, OG_pairs1_4_one2one,
      OG_pairs2_3_one2one, OG_pairs2_4_one2one, OG_pairs3_4_one2one
    )
    og_dts <- lapply(list(
      OG_pairs1_2, OG_pairs1_3, OG_pairs1_4, 
      OG_pairs2_3, OG_pairs2_4, OG_pairs3_4
    ), multiOGs)
    
    og_list <- lapply(sp_names, function(sp) {
      ogids <- grep(sp,lapply(og_dts, colnames))
      multis <- og_dts[ogids]
      ogsids <- setdiff(1:length(og_dts),ogids)
      singles <- ogs_dts[ogsids]
      
      all_multis <- na.omit(Reduce(function(...) merge(..., all=TRUE, by=sp, allow.cartesian=TRUE), multis))
      mdt1 <- merge.data.table(singles[[1]],singles[[2]])
      mdt2 <- merge.data.table(mdt1,singles[[3]], by=intersect(colnames(mdt1),colnames(singles[[3]])))
      dt <- merge.data.table(all_multis,mdt2,by=setdiff(sp_names,sp))
      setcolorder(dt,sp_names)
      dt
    })
    
    mdt <- unique(rbindlist(og_list))
    
    cross1=mdt[[sp1]]
    cross2=mdt[[sp2]]
    cross3=mdt[[sp3]]
    cross4=mdt[[sp4]]
    
    mdt[[sp1]] <- make.names(mdt[[sp1]], unique=TRUE)
    mdt[[sp2]] <- make.names(mdt[[sp2]], unique=TRUE)
    mdt[[sp3]] <- make.names(mdt[[sp3]], unique=TRUE)
    mdt[[sp4]] <- make.names(mdt[[sp4]], unique=TRUE)
    
    if (restrict_paralogs) {
      ll <- lapply(sp_names, function(sp) grepl("*\\.[3-9]",mdt[[sp]]))
      ol <- ll[[1]] | ll[[2]] | ll[[3]] | ll[[4]]
      remove_genes <- unique(str_remove(mdt[[sp1]][ol],"\\.\\d+"))
      remove_genes_length <- length(remove_genes)
      remove_genes_message <- paste(head(remove_genes,pmin(remove_genes_length,6)),collapse=", ")
      if (remove_genes_length>6) remove_genes_message <- paste0(remove_genes_message, ", ...")
      message("Removed ", remove_genes_length, " genes with more than three paralogs.")
      message(remove_genes_message)
      mdt <- mdt[grep(paste(remove_genes,collapse="|"),mdt[[sp1]],invert=TRUE),]
      
      cid <- !(cross1 %in% remove_genes)
      cross1=cross1[cid]
      cross2=cross2[cid]
      cross3=cross3[cid]
      cross4=cross4[cid]
    }
    
  }
  
  # Reorder matrices, quantile_norm and merge matrices
  if (any(c(class(sp1_fp),class(sp2_fp),class(sp3_fp))=="dgCMatrix")) {
    merged1=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]),sparse=TRUE)
    merged2=Matrix::Matrix(cbind(merged1,sp3_fp[cross3,]),sparse=TRUE)
    merged=Matrix::Matrix(cbind(merged2,sp4_fp[cross4,]),sparse=TRUE)
  } else {
    merged1=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]))
    merged2=data.matrix(cbind(merged1,sp3_fp[cross3,]))
    merged=data.matrix(cbind(merged2,sp4_fp[cross4,]))
  }
  rnm=mdt[[sp1]]
  rownames(merged)=rnm
  cnm=colnames(merged)
  
  if(quant_norm){
    if(any(class(merged)=="dgCMatrix")) {
      sparse=TRUE
      merged=as.matrix(merged)
      copymat=FALSE
    } else {
      sparse=FALSE
      copymat=TRUE
    }
    merged=tryCatch({
      preprocessCore::normalize.quantiles(merged,copy=copymat)
    }, error = function(e){
      quantile_normalisation(merged)
    })
    if (sparse==TRUE)
      merged=Matrix::Matrix(merged,sparse=TRUE)
    rownames(merged)=rnm
    colnames(merged)=cnm
  }
  
  # Select variable genes in ALL species
  nc1=ncol(sp1_fp)
  nc2=ncol(sp1_fp)+ncol(sp2_fp)
  nc3=nc2+ncol(sp3_fp)
  nc4=ncol(merged)
  if (!is.null(cross_fc_thrs)&!is.null(cross_n)) {
    var1= apply(merged[,1:nc1],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    var2= apply(merged[,c(nc1+1):nc2],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    var3= apply(merged[,c(nc2+1):nc3],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    var4= apply(merged[,c(nc3+1):nc4],1,function(x) 
      sort(x,decreasing=T)[cross_n]) > cross_fc_thrs
    top_cross=names(which(var1 & var2 & var3 & var4))
  } else {
    top_cross=rownames(merged)
  }
  
  # Return a list with both matrices
  out_m1=merged[,1:nc1]
  rownames(out_m1)=cross1
  out_m2=merged[,(nc1+1):nc2]
  rownames(out_m2)=make.names(cross2,unique=TRUE)
  out_m3=merged[,(nc2+1):nc3]
  rownames(out_m3)=make.names(cross3,unique=TRUE)
  out_m4=merged[,(nc3+1):nc4]
  rownames(out_m4)=make.names(cross4,unique=TRUE)
  
  ids2 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_2[[sp1]]==x)))
  top_cross_sp2=OG_pairs1_2[ids2,][[sp2]]
  ids3 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_3[[sp1]]==x)))
  top_cross_sp3=OG_pairs1_3[ids3,][[sp3]]
  ids4 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_4[[sp1]]==x)))
  top_cross_sp4=OG_pairs1_4[ids4,][[sp4]]
  
  csps <- list(
    merged=merged,
    OG_pairs1_2=OG_pairs1_2, OG_pairs1_3=OG_pairs1_3, OG_pairs1_4=OG_pairs1_4, 
    OG_pairs2_3=OG_pairs2_3, OG_pairs2_4=OG_pairs2_4, OG_pairs3_4=OG_pairs3_4,
    sp1=out_m1, sp2=out_m2, sp3=out_m3, sp4=out_m4,
    top_cross_sp1=top_cross, top_cross_sp2=top_cross_sp2, top_cross_sp2=top_cross_sp3, top_cross_sp4=top_cross_sp4
  )
  return(csps)
  
}


# # # # # # # # # # # # # # # # #
#                               #
#    CSPS ANALYSIS FUNCTIONS    #
#   take csps object as input   #
#                               #
# # # # # # # # # # # # # # # # #

#' Plots cross-species correlation matrix
#' 
#' Takes a csps object, and optionally a path to metacell annotation table files.
#' Saves heatmap of correlation between two species and return correlation matrix
#' used for plotting, as well as overlapping genes supporting each pairwise relation 
#' (defined by binarizing gene expression with `fc_thrs`).
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param output_file output png image filename
#' @param cluster_together logical, whether to cluster together metacells from 
#'   both species (i.e. rows and columns)
#' @param var_genes either character, vector of genes for comparison, or logical, 
#'   if TURE, using the csps-defined top-cross genes, if FALE, using all genes
#' @param annotation_file_1,annotation_file_2 metacell annotation tsv files with 
#'   three columns: metacell, cell_type, color 
#' @param cor_method character, corelation method to use, one of the following: 
#'   `c("pearson","spearman","kendall","jaccard")` (default: "jaccard")
#' @param cor_max numeric, color scaling max value (default: 1)
#' @param cex_dot numeric, dot size scaling factor (default: 1)
#' @param pow numeric, for plotting, raise the correlation to the power of 
#'   `pow` (default: 1)
#' @param fc_thrs numeric, fold change threshold for binarizing gene expression 
#'   (default: 1.2)
#'   when calculating Jaccard distance (default: 2); only used when `cor_method=='jaccard'`
#' @param reorder_sp2 logical, whether to plot reordered metacells in second species
#'   (default: FALSE)
#' @param reorder_by_ann1,reorder_by_ann2 logical, whether to plot metacells in order 
#'   in which they appear in annotation files (default: FALSE)
#' @param width numeric, figure height (in px)
#' @param height numeric, figure width (in px)
#' @param res numeric (default: NA)
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size umeric, size of annotation labels
#' 
#' @return a list with following elements: 
#'   1) `heatmap` complex heatmap object
#'   2) `cor_matrix` correlation matrix used for plotting
#'   3) `overlap_matrix` matrix with number of overlapping genes
#'   4) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
csps_plot_correlation_matrix=function(
  csps, cluster_together=FALSE, var_genes=TRUE, 
  annotation_file_1=NULL, annotation_file_2=NULL,
  cor_method="jaccard", fc_thrs=1.2, 
  output_file=NULL, width=3000, height=3000, res=NA,
  reorder_sps2=FALSE, reorder_by_ann1=FALSE, reorder_by_ann2=FALSE,
  pow = 1, cor_color = NULL, plot_type="heatmap", 
  annotation_size = 10, label_font_size = 12, cor_max=1, cex_dot=1,
  grid = TRUE,  annotation_grid_1 = TRUE, annotation_grid_2 = TRUE
){ 
  
  merged=cbind.data.frame(csps$sp1,csps$sp2)
  annotation=as.data.frame(substr(colnames(merged),1,4))
  colnames(annotation)="species"
  rownames(annotation)=colnames(merged)
  #annotation_colors=c("black","lightgrey")
  #names(annotation_colors)=c(substr(csps$top_cross_sp1[1],1,4),substr(csps$top_cross_sp2[1],1,4))
  
  if(length(var_genes)>1) {
    var_genes=glist
  } else if (var_genes==TRUE) {
    var_genes=csps$top_cross_sp1
  } else {
    var_genes=rownames(csps$merged)
  }
  #if(!is.null(glist)) var_genes=glist
  
  # heatmap color
  # htp_color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(1000)
  
  # ordering
  if(cluster_together){
    cor_joint=cor(merged[var_genes,],method=cor_method) ^ pow
    diag(cor_joint)=NA
    hc=hclust(dist(cor_joint),method="ward.D2")
    cor_mat <- cor_joint[hc$order,hc$order]
  } else {  
    m1 <- merged[var_genes,1:ncol(csps$sp1)]
    m2 <- merged[var_genes,(ncol(csps$sp1)+1):ncol(merged)]
    # binarize for overlap
    m1j <- m1; m1j[] <- 0 
    m2j <- m2; m2j[] <- 0
    m1j[!(m1<fc_thrs)] <- 1
    m2j[!(m2<fc_thrs)] <- 1
    m1j <- data.matrix(m1j)
    m2j <- data.matrix(m2j)
    # genes in common
    out_list <- vector("list",length = ncol(m1j))
    for (i in 1:ncol(m1j)) {
      intglist <- lapply(1:ncol(m2j), function(j) {
        ints <- which(m1j[,i] > 0 & m2j[,j] > 0)
        if (length(ints)==0) { NA } else { 
          pairs_ids <- match(rownames(m1j)[ints], rownames(csps$sp1))
          list(rownames(csps$sp1)[pairs_ids],rownames(csps$sp2)[pairs_ids])
        }
      })
      names(intglist) <- colnames(m2j)
      out_list[[i]] <- intglist
    }
    names(out_list) <- colnames(m1j)
    # correlation
    if (cor_method=="jaccard") {
      cor_m=jaccard(m1j,m2j) ^ pow
      # } else if (cor_method=="mahalanobis") {
      #   mt <- t(merged[var_genes,])
      #   mhdist <- biotools::D2.dist(data=mt,cov=cov(mt))
      #   hcl <- hclust(mhdist)
    } else if (cor_method=="kld") {
      mm=cbind(m1,m2)
      message("Calculating KLD, this might take a while")
      cor_m=calcKLD(mm)$mat ^ pow
    } else {
      cor_m=cor(m1,m2,method=cor_method) ^ pow
    }
    if(reorder_sps2){ 
      order_cols=names(sort(apply(cor_m,2,function(x) which.max(x)))) 
    } else { 
      order_cols=colnames(cor_m) 
    }
    cor_mat <- cor_m[,order_cols]
  }
  # plotting matrix
  col_ord <- colnames(cor_mat)
  row_ord <- rownames(cor_mat)
  cor_mat_name <- sprintf("%s\n", cor_method)
  
  
  # METACELL ANNOTATIONS
  if (!is.null(annotation_file_1)) {
    clust_anno_size  <- unit(annotation_size,"mm")
    annr <- fread(annotation_file_1,header=TRUE)
    setnames(annr,c("metacell","cell_type","color"))
    if (!any(annr$metacell %in% colnames(csps$merged))) {
      preffix1 <- str_extract(colnames(csps$merged)[1],"^[A-z]{4}_*")
      annr$metacell <- paste0(preffix1,annr$metacell)
    }
    if (reorder_by_ann1==TRUE) {
      tryCatch({
        row_ord <- annr$metacell
        cor_mat <- cor_mat[row_ord,]
        
      }, error=function(e) stop(
        "Different annotation and matrix names: ", 
        setdiff(row_ord,rownames(cor_mat)), " vs ",
        setdiff(rownames(cor_mat),row_ord),
      )
      )
    } 
    rid <- unlist(lapply(row_ord,match,table=annr$metacell))
    row_clusts <- annr[rid,]$cell_type
    row_clust_col <- annr[rid,]$color
    names(row_clust_col) <- row_clusts
    right_row_col_ha <- HeatmapAnnotation(
      which = "row", MC = row_clusts, col = list(MC = row_clust_col),
      lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
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
    left_row_col_ha <- right_row_col_ha
    
  }
  
  if(!is.null(annotation_file_2)){
    clust_anno_size  <- unit(annotation_size,"mm")
    annc <- fread(annotation_file_2,header=TRUE)
    setnames(annr,c("metacell","cell_type","color"))
    if (!any(annc$metacell %in% colnames(csps$merged))) {
      preffix2 <- str_extract(colnames(csps$merged)[ncol(csps$merged)],"^[A-z]{4}_*")
      annc$metacell <- paste0(preffix2,annc$metacell)
    }
    if (reorder_by_ann2==TRUE) {
      col_ord <- annc$metacell
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
      which = "column", col = list(MC = col_clust_col), MC = col_clusts,
      lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      border = TRUE, simple_anno_size = clust_anno_size,
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(annotation_size/3,"mm")
    )
  } else {
    top_column_col_ha <- HeatmapAnnotation(
      which = "column", lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
      border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
    )
    bottom_column_col_ha <- top_column_col_ha
  }
  
  # intersect
  cl_int <- overlap(m1j,m2j)
  ovrl_mat <- cl_int[row_ord,col_ord]
  sf <- max(ovrl_mat,na.rm=TRUE)
  cor_mat_sc <- ovrl_mat/sf
  
  # heatmap
  if (plot_type=="dotplot") {
    if(is.null(cor_color)) {
      cols<-c("white", "red")
      if (cor_method=="kld") col<-rev(cols)
      cor_color <- circlize::colorRamp2(c(0, cor_max), cols)
    }
    hm <- Heatmap(
      pmin(cor_mat,cor_max), col = cor_color, name = cor_mat_name, border = TRUE, 
      rect_gp = gpar(type = "none"),
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
      right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,  
      top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
      heatmap_legend_param = list(
        #col_fun = cor_color, at = cols_range,
        title = cor_mat_name, border = TRUE,
        legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
        title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
        labels_gp = gpar(fontsize = label_font_size)
      )
    )
    
  } else if (plot_type=="heatmap") {
    
    if(is.null(cor_color)) {
      #cor_color <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
      col_fun = circlize::colorRamp2(c(-2,0,2), c("blue","white","red"))
      cor_color <- col_fun(seq(-2,2,0.1))
      if (cor_method=="kld") cor_color <- rev(cor_color)
    }
    cv <- as.numeric(stringr::str_extract(as.character(packageVersion("ComplexHeatmap")),"\\d\\.\\d"))
    if (cv < 2) {
      hm <- Heatmap(
        pmin(cor_mat,cor_max), name=cor_mat_name,# col=cor_color, 
        rect_gp = gpar(col = "gray50", lwd = 0.2), 
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
        heatmap_legend_param = list(
          title = cor_mat_name, border = TRUE,
          legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
          title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
          labels_gp = gpar(fontsize = label_font_size)
        )
      )
    } else {
      hm <- Heatmap(
        pmin(cor_mat,cor_max), name=cor_mat_name,# col=cor_color, 
        rect_gp = gpar(col = "gray50", lwd = 0.2), border=TRUE, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,  
        top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
        heatmap_legend_param = list(
          title = cor_mat_name, border = TRUE,
          legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
          title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
          labels_gp = gpar(fontsize = label_font_size)
        )
      )
    }

  }
  if (!is.null(output_file)) {
    
    png(output_file,h=height,w=width,res=res)
    ht_opt(
      COLUMN_ANNO_PADDING=unit(10,"mm"), 
      ROW_ANNO_PADDING=unit(10,"mm"), 
      DIMNAME_PADDING=unit(5,"mm"),
      HEATMAP_LEGEND_PADDING=unit(15,"mm")
    )
    draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
    if(!is.null(annotation_file_1) & annotation_grid_1==TRUE){
      mat2 <- rbind(row_clust_col)
      change_clust_row <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
      decorate_heatmap_body(cor_mat_name, {
        for (i in change_clust_row) {
          grid.lines(x = c(0,1), y = 1-i/ncol(mat2), gp = gpar(lty = 1, lwd = 0.5))
        }
      })
    }
    if(!is.null(annotation_file_2) & annotation_grid_2==TRUE){
      mat2 <- rbind(col_clust_col)
      change_clust_col <- which(sapply(2:ncol(mat2), function(i) mat2[,i]!=mat2[,i-1]))
      decorate_heatmap_body(cor_mat_name, {
        for (i in change_clust_col) {
          grid.lines(x = i/ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
        }
      })
    }
    dev.off()
  }
  
  ht_opt(RESET = TRUE)
  
  return(list(heatmap=hm, cor_matrix=cor_mat, overlap_matrix=ovrl_mat, overlapping_genes=out_list))
}

#' Identify co-expressed genes in a defined set of metacells/cell types in two compared species
#' Optionally can be subseted by gene list (e.g. showing only TFs)
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param sp1_focus,sp2_focus character, one or more metacell (or cell type) names 
#'   (i.e. columns in `csps$merged`)
#' @param glist character, optional gene list for subsetting the genes (default: NULL)
#' @param fc_thrs numeric, gene fold change threshold (default: 1.5)
#' @param gene_annot_sp1_fn,gene_annot_sp1_fn character, path to gene annotation file for 
#'   the first and the second species in `csps`
#' @param out_fn character, path to png file where expression heatmap will be saved, same filename
#'   will be used to generate other output files (barplot figure and txt file)
#'   (default: NULL, saves files to current working directory)
#' 
csps_identify_coexpressed_genes=function(
  csps, sp1_focus, sp2_focus,
  glist=NULL, fc_thrs=1.5,
  gene_annot_sp1_fn, gene_annot_sp2_fn,
  out_fn=NULL, width_cex=20, height_cex=20, res=NA,
  fc_max=3, annotation_size = 10, label_font_size=12
){ 
  
  if (is.null(out_fn)) {
    out_bfn=paste0("Coexpressed_genes_",paste(c(sp1_focus,sp2_focus),collapse="_"))
  } else {
    out_bfn=str_remove(out_fn,".png")
  }
  message(out_bfn)
  
  annot_sp1=read.table(gene_annot_sp1_fn,h=T,row.names=1,sep="\t",quote="",stringsAsFactors=F,fill=T)
  annot_sp2=read.table(gene_annot_sp2_fn,h=T,row.names=1,sep="\t",quote="",stringsAsFactors=F,fill=T)
  bckg_sp1=setdiff(colnames(csps$sp1),sp1_focus)
  bckg_sp2=setdiff(colnames(csps$sp2),sp2_focus)
  
  f_sp1=apply(csps$merged[,sp1_focus,drop=F],1,function(x) sort(x,decreasing=T)[1]) >= fc_thrs & apply(csps$merged[,bckg_sp1],1,median) < fc_thrs
  f_sp2=apply(csps$merged[,sp2_focus,drop=F],1, function(x) sort(x,decreasing=T)[1]) >= fc_thrs & apply(csps$merged[,bckg_sp2],1,median) < fc_thrs
  
  sps1_selected_ids=names(which(f_sp1 & f_sp2))
  hc=hclust(dist(cor(t(csps$merged[sps1_selected_ids,]))),method="ward.D2")
  gene_order=sps1_selected_ids[hc$order]
  sps2_selected_ids=csps$og_pairs[match(gene_order,csps$og_pairs[[1]])][[2]]
  
  output_table=cbind.data.frame(annot_sp1[str_remove(gene_order,"\\.\\d"),],annot_sp2[sps2_selected_ids,])
  if (!is.null(out_fn))
    write.table(output_table,file=paste0(out_bfn,"_table.txt"),quote=F,sep="\t",col.names=F,row.names=T)
  
  shades2=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
  hmat <- pmin(csps$merged[gene_order,], fc_max)
  hmtitle=sprintf("Coexspressed genes in %s and %s", sp1_focus, sp2_focus)
  col_ord <- colnames(hmat)
  row_ord <- str_remove(rownames(hmat),"\\.\\d") #rownames(hmat)
  row_ord_1 <- annot_sp1[row_ord,1]
  row_ord_1[nchar(row_ord_1)>30] <- paste0(substr(row_ord_1[nchar(row_ord_1)>30],1,27),"...")
  row_ord_1 <- str_replace(row_ord_1,"\"\"","")
  row_ord_2 <- annot_sp1[row_ord,2]
  row_ord_2[nchar(row_ord_2)>30] <- paste0(substr(row_ord_2[nchar(row_ord_2)>30],1,27),"...")
  row_ord_2 <- str_replace(row_ord_2,"\"\"","")
  top_column_col_ha <- HeatmapAnnotation(
    which = "column",
    lab = anno_text(which = "column", just = "left", location = unit(0, 'npc'), col_ord, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  bottom_column_col_ha <- HeatmapAnnotation(
    which = "column",
    lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  left_row_col_ha <- HeatmapAnnotation(
    which = "row",
    lab = anno_text(which = "row", row_ord_2, just = "right", location = unit(1, 'npc'), gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  right_row_col_ha <- HeatmapAnnotation(
    which = "row",
    lab = anno_text(which = "row", row_ord_1, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  hm <- Heatmap(
    hmat, col = shades2, name = "expression", border = TRUE,
    rect_gp = gpar(col = "gray50", lwd = 0.2), 
    show_row_names = FALSE, show_column_names = FALSE, 
    row_title = hmtitle, row_title_gp = gpar(fontsize = 2*label_font_size),
    cluster_rows = FALSE, cluster_columns = FALSE, 
    left_annotation = left_row_col_ha, right_annotation = right_row_col_ha, 
    bottom_annotation = bottom_column_col_ha, top_annotation = top_column_col_ha, 
    heatmap_legend_param = list(
      title = "expression\n", border = TRUE,
      legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
      title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
      labels_gp = gpar(fontsize = label_font_size)
    )
  )
  
  if(!is.null(out_fn)) {
    
    png(paste0(out_bfn,"_boxplot.png"),h=500,w=2000)
    par(mar=c(15,7,2,2))
    boxplot(log2(csps$merged[sps1_selected_ids,]),las=2,pch=20,outcol=alpha("black",0.2),ylab="log2FC")
    abline(h=0)
    dev.off()
    
    height=max(length(gene_order)*height_cex,2000)
    width=ncol(csps$merged)*width_cex
    png(paste0(out_bfn,"_heatmap.png"),height=height,width=width)
    ht_opt(
      COLUMN_ANNO_PADDING=unit(5,"mm"), 
      ROW_ANNO_PADDING=unit(5,"mm"), 
      DIMNAME_PADDING=unit(5,"mm"),
      HEATMAP_LEGEND_PADDING=unit(10,"mm"),
      TITLE_PADDING=unit(10,"mm")
    )
    draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
    dev.off()
    
  }
  return(list(heatmap=hm, output_table=output_table))
}


#' Select one or more MCs in the first species and project expression 
#' of top genes across both species in comparison. 
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param focus_mcs character, one or more metacell (or cell type) names (i.e. columns in `csps$merged`)
#' @param gene_annot_sp1_fn character, path to gene annotation file for the first species 
#'   in `csps`
#' @param sps1_tfs_fn character, path to the matrix of gene fold change
#'   in metacells (or cell types) for the first species in `csps` (usually metacell 
#'   footprint matrix, `mc@mc_fp`)
#' @param max_n_genes integer, number of top genes to project (default: 100)
#' @param fc_thrs numeric, fold change threshold (default: 2)
#' @param out_fn character, path to png file where expression heatmap will be saved
#'   (default: NULL)
#' @param width_cex numeric, figure height scalling
#' @param height_cex numeric, figure width scalling
#' @param res numeric (default: NA)
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size numeric, size of annotation labels
#' @param fc_max numeric, fc color scaling max value (default: 3)
#' 
#' @return list with following elements
#'   1) `heatmap` ComplexHeatmap object
#'   2) `expr_mat` matrix with the expression of top genes
#'   
csps_focused_coexpression=function(
  csps, focus_mcs, gene_annot_sp1_fn, sps1_tfs_fn, 
  max_n_genes=100, fc_thrs=2, 
  out_fn=NULL, width_cex=20, height_cex=20, res=NA,
  fc_max=3, annotation_size = 10, label_font_size=12
) {
  
  tfs_sp1=read.table(sps1_tfs_fn,h=T,row.names=1,sep="\t",quote="",stringsAsFactors=F,fill=T)
  annot_sp1=read.table(gene_annot_sp1_fn,row.names=1,sep="\t",quote="",stringsAsFactors=F,fill=T,header=F)
  
  ffall=sort(apply(csps$merged[,focus_mcs,drop=F],1,median))
  #ff=tail(ffall[!grepl("\\.\\d",names(ffall))],max_n_genes)
  ff=tail(ffall,max_n_genes)
  focus_genes=names(which(ff>fc_thrs))
  
  shades2=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
  #label_col=ifelse(focus_genes %in% rownames(tfs_sp1),"red","black")
  hmtitle=sprintf("Top %s genes in %s",max_n_genes, focus_mcs)
  hmat <- pmin(csps$merged[rev(focus_genes),],fc_max)
  col_ord <- colnames(hmat)
  row_ord <- str_remove(rownames(hmat),"\\.\\d") #rownames(hmat)
  row_ord_1 <- annot_sp1[row_ord,1]
  row_ord_1[nchar(row_ord_1)>30] <- paste0(substr(row_ord_1[nchar(row_ord_1)>30],1,27),"...")
  row_ord_1 <- str_replace(row_ord_1,"\"\"","")
  row_ord_2 <- annot_sp1[row_ord,2]
  row_ord_2[nchar(row_ord_2)>30] <- paste0(substr(row_ord_2[nchar(row_ord_2)>30],1,27),"...")
  row_ord_2 <- str_replace(row_ord_2,"\"\"","")
  top_column_col_ha <- HeatmapAnnotation(
    which = "column",
    lab = anno_text(which = "column", just = "left", location = unit(0, 'npc'), col_ord, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  bottom_column_col_ha <- HeatmapAnnotation(
    which = "column",
    lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  left_row_col_ha <- HeatmapAnnotation(
    which = "row",
    lab = anno_text(which = "row", row_ord_2, just = "right", location = unit(1, 'npc'), gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  right_row_col_ha <- HeatmapAnnotation(
    which = "row",
    lab = anno_text(which = "row", row_ord_1, gp = gpar(fontsize = label_font_size)),
    border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
  )
  hm <- Heatmap(
    hmat, col = shades2, name = "expression", border = TRUE,
    rect_gp = gpar(col = "gray50", lwd = 0.2), 
    show_row_names = FALSE, show_column_names = FALSE, 
    row_title = hmtitle, row_title_gp = gpar(fontsize = 2*label_font_size),
    cluster_rows = FALSE, cluster_columns = FALSE, 
    left_annotation = left_row_col_ha, right_annotation = right_row_col_ha, 
    bottom_annotation = bottom_column_col_ha, top_annotation = top_column_col_ha, 
    heatmap_legend_param = list(
      title = "expression\n", border = TRUE,
      legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
      title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
      labels_gp = gpar(fontsize = label_font_size)
    )
  )
  
  if (!is.null(out_fn)) {
    out_fn_prefix <- str_remove(out_fn,"\\.png")
    out_bfn=paste0(out_fn_prefix,"_",strtrim(paste(focus_mcs,collapse="_"),50),".png")
    height=max(length(focus_genes)*height_cex,2000)
    width=ncol(csps$merged)*width_cex
    
    png(out_bfn,h=height,w=width,res=NA)
    
    ht_opt(
      COLUMN_ANNO_PADDING=unit(5,"mm"), 
      ROW_ANNO_PADDING=unit(5,"mm"), 
      DIMNAME_PADDING=unit(5,"mm"),
      HEATMAP_LEGEND_PADDING=unit(10,"mm"),
      TITLE_PADDING=unit(10,"mm")
    )
    draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
    dev.off()
    
  }
  
  return(list(heatmap=hm, expr_mat=hmat))
}


# # # # # # # # # # # # # # # # #
#                               #
#    DOWNSTREAM FUNCTIONS       #
#    take matrix as input       #
#                               #
# # # # # # # # # # # # # # # # #

#' Identify genes with conserved expression in broad cell types of two species 
#' by calulating Jaccard index for expression in cell types
#' 
#' @param mc_fp_1,mc_fp_2 cell type gene expression matrix, with genes in rows 
#'   and cell types in columns. At least some of the column names should either 
#'   contain the pattern specified with `cts` argument, or start with common broad 
#'   cell type particle (four letters) - this will be used for overlap calculation.
#' @param cts pattern to use for matching cell typs; if NULL (default), broad cell type
#'   particle is used (i.e. the first four letters of expression matrix column names)
#' @param fc_thrs numeric, fold change threshold for genes selection
#' @param jacc_thrs numeric, Jaccard index threshold for genes selection
#' @return data.table with conserved cell types for each gene
#' 
ctConservedGenes <- function(mc_fp_1, mc_fp_2, cts=NULL, fc_thrs=1.7, jacc_thrs=0.5) {
  
  # jaccard distance for overlap in cell types
  .jacc_ct <- function(cell_types_1, cell_types_2) {
    cell_types_int <- intersect(cell_types_1,cell_types_2)
    cell_types_uni <- c(cell_types_1,cell_types_2)
    length_int <- length(cell_types_uni[cell_types_uni %in% cell_types_int])
    length_uni <- length(cell_types_uni)
    length_int/length_uni
  }
  
  # cts for all genes
  gen_ct_1 <- apply(mc_fp_1, 1, function(x) {
    id <- x>fc_thrs
    if (is.null(cts)) {
      cts <- substr(colnames(mc_fp_1)[id],1,4)
    } else {
      cts <- str_extract(colnames(mc_fp_1)[id],pattern=paste(cts,collapse="|"))
    }
    if (length(cts>0)) {cts} else {NULL}
  })
  gen_ct_1 <- gen_ct_1[-which(lapply(gen_ct_1,is.null)==T)]
  gen_ct_1 <- sapply(gen_ct_1, function(g) g[!is.na(g)], USE.NAMES=TRUE, simplify=FALSE)
  gen_ct_2 <- apply(mc_fp_2, 1, function(x) {
    id <- x>fc_thrs
    if (is.null(cts)) {
      cts <- substr(colnames(mc_fp_2)[id],1,4)
    } else {
      cts <- str_extract(colnames(mc_fp_2)[id],pattern=paste(cts,collapse="|"))
    }
    
    if (length(cts>0)) {cts} else {NULL}
  })
  gen_ct_2 <- gen_ct_2[-which(lapply(gen_ct_2,is.null)==T)]
  gen_ct_2 <- sapply(gen_ct_2, function(g) g[!is.na(g)], USE.NAMES=TRUE, simplify=FALSE)
  
  # common genes
  cg <- intersect(names(gen_ct_1),names(gen_ct_2))
  message(length(cg), " genes in intersect")
  
  # found common cts for genes
  ctj <- sapply(cg, function(g) {
    jacc <- .jacc_ct(gen_ct_1[[g]], gen_ct_2[[g]])
    if (length(jacc)>0 & !is.na(jacc)) {
      if (jacc>jacc_thrs) {
        table(c(gen_ct_1[[g]],gen_ct_2[[g]]))
      } else {NULL}
    } else {NULL}
  }, USE.NAMES=TRUE, simplify=FALSE)
  nulls <- lapply(ctj,is.null)==T
  if (any(nulls))
    ctj <- ctj[-which(nulls)]
  message(length(ctj))
  ctj_tbl <- table(sapply(ctj,length))
  message(sprintf(
    "Found %s conserved genes (Jaccard > 0.5); 
    \n%s genes with unique conserved ct;
    \n%s genes with non-unique conserved cts",
    length(ctj), ctj_tbl["1"], sum(ctj_tbl)-ctj_tbl["1"] 
  ))
  dt <- rbindlist(lapply(ctj, function(x) {
    x <- as.matrix(unclass(x))
    data.table(ct=rownames(x),occurence=x[,1])
  }),idcol="gene")
  dt
}

#' Identify genes with conserved expression in selected columns of expression matrix 
#' (metacells, cell types) by expression fold change thresholding.
#' 
#' @param mc_fp cell type gene expression matrix, with genes in rows and cell types in columns
#' @param cols integer or character specifiying poisitions or names of columns to use
#' @param feature_in_thrs (default: 1.5)
#' @param feature_out_thrs (default: 1)
#' @param method,methodbg character specifying how to summarize gene expression in selected columns, 
#'   one of "absolute" (default) or "median"
#' @param abs_leakyness,abs_leakynessbg numeric, percent of selected columns or background columns 
#'   in which the expression value can be below or above the specified threshold, respectively;
#'   this is only used if `method`` or `methodbg`` is "absolute"; (default 0 and 0.05, respectively)
#' 
commonGenes <- function(
  mc_fp, cols, feature_in_thrs = 1.5, feature_out_thrs = 1,
  method="absolute", methodbg="absolute", abs_leakyness=0, abs_leakynessbg=0.05
) {
  
  # check that given cols are in the matrix and get in and out cols
  if (all(class(cols)=="integer")) {
    cols=colnames(mc_fp)[cols]
  } else if (all(class(cols)=="character")) {
    cols=intersect(cols,colnames(mc_fp))
  }
  cols_out <- setdiff(colnames(mc_fp),cols)
  
  # function to calculate leakynes
  .calc_leakyness <- function(abs_leakyness,n) {
    if (abs_leakyness<1) {
      pmin(round(c(1-abs_leakyness)*n), n)
    } else {
      pmin(round(n-abs_leakyness), n)
    }
  }
  # featues in and out
  if (method=="median") {
    fin=apply(mc_fp[,cols,drop=F],1,median) > feature_in_thrs
    fin_inv=apply(mc_fp[,cols_out,drop=F],1,median) > feature_in_thrs
  } else if (method=="absolute") {
    fin=apply(mc_fp[,cols,drop=F], 1, function(x) 
      !(sum(x > feature_in_thrs) < .calc_leakyness(abs_leakyness,n=length(cols)))
    )
    fin_inv=apply(mc_fp[,cols_out,drop=F], 1, function(x) 
      !(sum(x > feature_in_thrs) <  .calc_leakyness(abs_leakyness,n=length(cols_out)))
    )
  } else {
    stop("method should be either 'median' or 'absolute'")
  }
  if (methodbg=="median") {
    fout=apply(mc_fp[,cols_out,drop=F],1,median) < feature_out_thrs
    fout_inv=apply(mc_fp[,cols,drop=F],1,median) < feature_out_thrs
  } else if (methodbg=="absolute") {
    fout=apply(mc_fp[,cols_out,drop=F], 1, function(x) 
      !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(cols_out)))
    )
    fout_inv=apply(mc_fp[,cols,drop=F], 1, function(x) 
      !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(cols)))
    )
  } else {
    stop("methodbg should be either 'median' or 'absolute'")
  }
  
  f_in=which(fin & fout)
  f_out=which(fin_inv & fout_inv)
  
  list(genes_in=names(f_in), genes_out=names(f_out))
  
}


#' Plot violin plots of expression of selected genes in selected cell type, versus all other cell types
#' 
#' @param mc_fp expression matrix, rows are genes, columns are cell types with species prefix
#' @param ct one of colnames of mc_fp
#' @param species vector of species four-letter abbreviations
#' @param logical, whether cell type names are preceeded by species four-letter abbreviations
#' @param genes cell type specific genes 
#' @param scale.fc.max numeric, scale gene fc to this max value for plotting (default NULL)
#' @param sign.test function to calculate significance test, either `t.test` or `wilcox.test`
#' @param sign.label "pval","padj","p.adj.signif"
#' @param title logical, show cell type names as title
#' @param x.names logical, hide cell type names on x axis
#' @param x.names.replacement named character, optional replacement pattern for cell type labes on x axis
#' 
csps_plot_groupped_gene_expression <- function(
  mc_fp, ct, species=c("Spis","Nvec","Xesp","Hvul"), species_cell_type_names=TRUE, 
  genes, min.genes=3, scale.fc.max=NULL, sign=TRUE, sign.test=t.test, sign.label="p.adj.signif",
  col, title=TRUE, x.names=TRUE, x.names.angle=30, x.names.replacement=NULL,
  text.size=12, label.size=12, title.size=12, tip.length=0, step.increase=0.1
) {
  
  spreg=paste(species,collapse="|")
  ctl <- paste(ct,collapse=" ")
  message ("Plotting expression for ", ctl)
  
  # transform matrix to data.table
  mc_fp_dt=melt.data.table(as.data.table(
    mc_fp,keep.rownames="gene"
  ),id.vars="gene",variable.name="cell_type",value.name="fc")
  mc_fp_dt[,species:=str_extract(cell_type,spreg)]
  mc_fp_dt[,species:=factor(species,levels=unique(species))]
  if (!species_cell_type_names) 
    mc_fp_dt[,cell_type:=str_remove(cell_type,sprintf("(%s)_",spreg))]
  
  # genes to plot
  if (length(genes)>0) {
    features_in_mc_fp=mc_fp_dt[gene %in% genes]
    features_in_mc_fp[,selected_cell_type:="other"]
    features_in_mc_fp[cell_type%in%ct, selected_cell_type:=cell_type]
    features_in_mc_fp[,selected_cell_type:=factor(selected_cell_type,levels=c(unique(ct),"other"))]
    uniqgenes <- min(unique(features_in_mc_fp[,.N,.(cell_type,species)]$N))
    if (is.infinite(uniqgenes)) uniqgenes <- 0
    if (uniqgenes > min.genes) {
      
      stattest=copy(features_in_mc_fp)
      
      # species must have at least two levels
      ss=stattest[,unique(.SD[,.(species,selected_cell_type)]),.(species)][,.N,species][N>1]$species
      stattest=stattest[species %in% ss]
      
      # cell types to compare to "other"
      levels_ct <- levels(stattest$selected_cell_type)
      levels_compare <- levels_ct[levels_ct != "other"]
      
      dt <- stattest[,.(fc=mean(fc)),.(gene,selected_cell_type,species)]
      dt <- dt[,cell_type:=as.character(selected_cell_type)]
      st <- dt[,rstatix::pairwise_sign_test(data=.SD, formula=fc~cell_type, ref.group="other"),species]
      st[,species:=factor(species,levels=unique(species))]
      
      # plot
      if (!"other" %in% names(col)) {
        col <- append(col,"gray60")
        names(col)[length(col)] <- "other"
      }
      features_plot <- dt[species %in% ss]
      if (is.null(scale.fc.max)) {
        dt[,y:=max(fc)+max(fc)*0.1,species]
      }else {
        features_plot[,fc:=pmin(fc,scale.fc.max)]
        st[,y:=scale.fc.max+scale.fc.max*0.1]
      }
      features_plot[,selected_cell_type_labels:=str_remove_all(selected_cell_type,sprintf("(%s)_*",spreg))]
      cell_types_labels <- features_plot$selected_cell_type_labels
      if (!is.null(x.names.replacement) & length(x.names.replacement)>0)
        cell_types_labels <- str_replace_all(cell_types_labels,x.names.replacement)
      names(cell_types_labels) <- features_plot$selected_cell_type
      features_plot[,selected_cell_type:=droplevels(selected_cell_type)]
      selected_cell_type <- unique(features_plot$selected_cell_type)
      message("groups: ",paste(selected_cell_type,collapse=", "))
      message("colors: ",paste(col[selected_cell_type],collapse=", "))
      gp=ggplot(features_plot,aes(selected_cell_type,fc,color=selected_cell_type)) + 
        facet_grid(.~species,scales="free_x",space="free_x",switch="x",drop=TRUE) + #strip.position="bottom"
        geom_jitter(size=1, alpha=0.5, width=0.4) +
        geom_violin(alpha=1, lwd=0.2) +
        scale_color_manual(values=col) + scale_fill_manual(values=col) + 
        scale_x_discrete(labels=cell_types_labels) +
        geom_text(aes(y=fc*1.2, label="")) + 
        theme(
          legend.position="none", text=element_text(size=text.size),
          axis.line.x=element_blank(), #axis.ticks.x=element_blank(), 
          axis.text.x=element_text(angle=x.names.angle,hjust=1), # vjust=0.5 if angle=90
          plot.title=element_text(size=title.size),
          strip.background=element_blank(), strip.placement="outside",
          panel.border=element_blank()
        ) + 
        labs(x=NULL,y="gene fold change",title=ctl)
      if (sign==TRUE) {
        gp <- gp + stat_pvalue_manual(
          st, label=sign.label, y.position="y",
          step.increase=step.increase, step.group.by="species",
          label.size = label.size/3.88, tip.length=tip.length
        )
      }
      if (title!=TRUE) {
        gp <- gp + theme(plot.title=element_blank())
      }
      if (x.names!=TRUE)
        gp <- gp + theme(axis.text.x=element_blank())
    } else {
      gp <- NULL
      message(sprintf("Not enough genes (%s provided, more than %s required)", uniqgenes, min.genes))
    }
  } else {
    gp <- NULL
    message(sprintf("Not enough genes in expression data (%s provided, more than %s required)", length(genes), min.genes))
  }
  return(gp)
}


#' Plot cross-species chord diagram
#' 
#' @param mat matrix with similarities to plot on the circos; colnames and 
#'   rownames should start with four-letter species abbreviations
#' @param sp character, four-letter abbreviations of the species for which to 
#'   plot the circos, (at least some of) `colnames(mat)` and `rownames(mat)` 
#'   should start with these
#' @param revert logical, wheter to use `1/mat` for plotting links, set this to 
#'   TRUE if the values in matrix are a measure of divergence, and to FALSE if 
#'   they measure similarity
#' @param threshold numeric, between 0 and 1, a threshold to use for selecting 
#'   links to be shown on circos plot - the values in the matrix ae scaled to 
#'   [0,1] range and only those with scaled value above specified threshold are 
#'   plotted
#' @param name.suffix character, optional suffix to be appended to the plot
#'   filename which by default will include species names
#' @param outdir character, directory in which the plot will be saved
#' @param width numeric, width of the plot in inches
#' @param height numeric, height of the plot in inches
#' @param sectors.order charcter vector with names of sectors to be plotted 
#'   on the circos, should be in `c(colnames(mat)`, `rownames(mat))`
#' @param sectors.labels named character, optional labels for sectors, the names
#'   should be `sectors.order`
#' @param sectors.colors named character, optional colors for sectors, the names 
#'   should be `sectors.order`
#' @param sectors.groups named factor or character, indicates grouping of sectors, 
#'   levels of factor determine ordering on the circos, and the names should be 
#'   `sectors.order`
#' @param sectors.groups.labels named character, optional labels for groups of 
#'   sectors, the names should be `unique(sectors.groups)`
#' @param start.degree numeric
#' @param annotation.track logical
#' @param annotation.names logical
#' 
csps_plot_chord_diagram <- function(
  mat, sp, revert=FALSE, threshold=NULL, threshold.quantile=NULL,
  scale.values=TRUE, regularize.values=TRUE, reg.param=1,
  save.plot=TRUE, save.data=TRUE, name.suffix="", outdir=".", width=18, height=18, mar=c(10,5,10,5),
  sectors.order=NULL, sectors.labels=NULL, sectors.colors=NULL, grid.border=NULL,
  sectors.groups=NULL, sectors.groups.labels=NULL,
  sector.scale=TRUE, sector.width.link=FALSE, remove.empty.sectors=FALSE, self.links=FALSE, 
  sector.groups.padding=c(0.1,0,1.8,0), sector.groups.col="gray98", sector.groups.border="gray88", 
  start.degree=360, small.gap=0, big.gap=5, annotation.track=TRUE, annotation.names=TRUE,
  sector.text.cex=1, sector.groups.text.cex=2, sector.groups.text.vjust=-6,
  title.main=NULL, title.sub=NULL, title.main.cex=2, title.sub.cex=2
) {
  # functions
  col_fun = function(x,min,max) {
    f=circlize::colorRamp2(c(min,max), c("white","red"))
    f(x)
  }
  transparency_fun = function(x,n=4,scale=TRUE,pseudocount=0.0001) {
    if(scale==TRUE) {
      x=(x-min(x))/(max(x)-min(x))
    }
    x=-n*log(x+pseudocount)
    (x-min(x))/(max(x)-min(x))
  }
  highlightSector = function(
    group.name,group,track.index=1,
    col="gray98",border="black",padding=c(0,0,0,0),
    text=group.name,text.col="black",text.cex=1,text.vjust=0.5,niceFacing=TRUE
  ) {
    highlight.sector(
      sector.index=names(group[group==group.name]), track.index=1, 
      col=col, border=border, padding=padding,
      text=text, text.vjust=text.vjust, text.col=text.col, cex=text.cex,
      niceFacing=niceFacing, facing="bending.inside"
    )
  }
  # data
  if (is.null(sectors.order))
    sectors.order=colnames(mat)
  if (is.null(sectors.labels)) {
    sectors.labels=sectors.order
    names(sectors.labels)=sectors.order
  }
  if (is.null(sectors.colors)) {
    sectors.colors=rand_color(length(sectors.labels))
    names(sectors.colors)=sectors.order
  }
  dt=as.data.table(mat,keep.rownames="ct1")
  dtm=melt.data.table(dt,id.vars="ct1",variable.name="ct2",value.name="value")
  dtm[,ct1:=factor(ct1,levels=sectors.order[sectors.order %in% ct1])]
  dtm[,ct2:=factor(ct2,levels=sectors.order[sectors.order %in% ct2])]
  dtm[,sp1:=stringr::str_extract(ct1,"[A-z]{4}")]
  dtm[,sp2:=stringr::str_extract(ct2,"[A-z]{4}")]
  if (revert==TRUE) {
    dtm[,value:=1/value]
    dtm[is.infinite(value),value:=0]
  }
  #setcolorder(dtm,c("ct1","ct2","value","sp1","sp2"))
  dtm[,linkvalue:=value]
  if (regularize.values==TRUE) {
    regdt=dtm[sp1==sp2][ct1!=ct2][,.(rg=max(value)),ct1]#[,.(rg=max(value)),.(sp1,sp2)]
    dtm[regdt,on="ct1",rg:=i.rg]
    dtm[,linkvalue:=linkvalue/(reg.param*rg)]
  }
  if (scale.values==TRUE) 
    dtm[,linkvalue:=linkvalue/max(linkvalue)]
  if (self.links==TRUE) {
    dtms=dtm[sp1==sp2]
  } else {
    dtms=NULL
  }
  if (length(sp)>4) {
    # use the first group as a reference and plot only links from it to others
    dtms=rbindlist(list(
      dtms,
      dtm[grepl(sp[1],dtm[["ct1"]])][grepl(paste(sp[2:length(sp)],collapse="|"),dtm[["ct2"]])]
    ))
  } else {
    if (!is.na(sp[2])) { # plot links: 1-2
      dtms=rbindlist(list(
        dtms,dtm[grepl(sp[1],dtm[["ct1"]]) & grepl(sp[2],dtm[["ct2"]])]
      ))
    } 
    if (!is.na(sp[3])) { # plot links: 1-2, 1-3, 2-3
      dtms=rbindlist(list(
        dtms,
        dtm[grepl(sp[1],dtm[["ct1"]]) & grepl(sp[3],dtm[["ct2"]])],
        dtm[grepl(sp[2],dtm[["ct1"]]) & grepl(sp[3],dtm[["ct2"]])]
      ))
    }
    if (!is.na(sp[4])) { # plot links: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
      dtms=rbindlist(list(
        dtms,
        dtm[grepl(sp[2],dtm[["ct1"]]) & grepl(sp[4],dtm[["ct2"]])],
        dtm[grepl(sp[3],dtm[["ct1"]]) & grepl(sp[4],dtm[["ct2"]])]
      ))
    }
  }
  dtcd=dtms[,.(ct1,ct2,linkvalue)][order(ct1,ct2)]
  # select links above threshold
  if (!any(!is.null(threshold),!is.null(threshold.quantile)) |
      all(is.null(threshold),is.null(threshold.quantile))) {
    stop("You must specify either threshold or threshold.quantile!")
  }
  if (!is.null(threshold.quantile)) {
    threshold=quantile(dtcd$linkvalue,threshold.quantile)
  }
  dtcds=dtcd[linkvalue>threshold]
  # add empty space at the link origin for links below threshold
  dtcds=rbindlist(list(
    dtcds,
    dtcd[!ct1%in%dtcds$ct1][,.SD[order(linkvalue,decreasing=TRUE)][1],ct1][,.(ct1,ct2,linkvalue)]
  ))
  # add empty space at the link destination for links below threshold
  if (remove.empty.sectors==FALSE)
    dtcall=rbindlist(list(
      dtcds,
      dtcd[!ct2%in%dtcds$ct2][,.SD[order(linkvalue,decreasing=TRUE)][1],ct2][,.(ct1,ct2,linkvalue)]
    ))
  dtcall=dtcall[order(ct1,ct2)]
  dtcall[,ct1:=as.character(ct1)][,ct2:=as.character(ct2)]
  # link color by link value
  # dtcall[,col:=col_fun(linkvalue,min(linkvalue),max(linkvalue))]
  # link color by originating sector, hide links below threshold
  dtcall[,col:=sectors.colors[as.character(ct1)]]
  dtcall[!(linkvalue>threshold),col:="white"]
  # remove below threshold links (i.e. make above threshold links the same width as sector)
  if (sector.width.link==TRUE) {
    single_link_cts=dtcall[,.N,ct1][N==1]$ct1
    empty_cts=dtcall[ct1 %in% single_link_cts][col=="white"]$ct1[1]
    dtcall[col=="white",ct1:=empty_cts]
  }
  # remove sectors that only have below threshold links
  if (remove.empty.sectors==TRUE) {
    multiple_link_cts=dtcall[,.N,ct1][N>1]$ct1
    dtcall=dtcall[!(ct1 %in% multiple_link_cts & col=="white")]
  }
  # set transparency
  dtcall[,transparency:=transparency_fun(linkvalue)]
  dtcall[!(linkvalue>threshold),transparency:=1]
  # set order
  dtcall[,zindex:=rank(linkvalue)]
  setcolorder(dtcall,c("ct1","ct2","linkvalue","col","transparency","zindex"))
  if (save.data==TRUE)
    saveRDS(
      dtcall, 
      file.path(outdir,sprintf("circos_%s_%s.RDS",paste0(sp,collapse="_"),name.suffix))
    )
  
  # circos
  if (save.plot==TRUE) {
    pdf(
      file.path(outdir,sprintf("circos_%s_%s.pdf",paste0(sp,collapse="_"),name.suffix)),
      w=width,h=height,useDingbats=TRUE
    )
    par(mar=mar)
  }
  circos.par(
    "canvas.xlim"=c(-2, 2),"canvas.ylim"=c(-2, 2),
    "track.height"=0.2, "track.margin"=c(0,0), cell.padding=c(0,1,0,1),
    "start.degree"=start.degree
  )
  annotationTrack=ifelse(annotation.track,"grid",NULL)
  chordDiagram(
    dtcall, grid.col=sectors.colors, grid.border=grid.border,
    annotationTrack=annotationTrack, small.gap=small.gap, big.gap=big.gap, 
    directional=2, #direction.type=c("arrows"), link.arr.type = "big.arrow", 
    #annotationTrackHeight = c(0.05,-0.05), 
    scale=sector.scale, col=dtcall$col, transparency=dtcall$transparency+0.1,
    link.zindex=dtcall$zindex, group=sectors.groups, preAllocateTracks=1,
  )
  groupchord=sectors.groups[names(sectors.groups)%in%c(as.character(dtcall$ct1),as.character(dtcall$ct2))]
  for (group.name in unique(c(dtm$sp1,dtm$sp2))) {
    if (!is.null(sectors.groups.labels)) {
      sectorlabs=sectors.groups.labels[group.name]
    } else {
      sectorlabs=group.name
    }
    highlightSector(
      group.name=group.name, group=groupchord,
      padding=sector.groups.padding, col=sector.groups.col, border=sector.groups.border,
      text=sectorlabs, text.vjust=sector.groups.text.vjust, text.cex=sector.groups.text.cex
    )
  }
  circos.trackPlotRegion(track.index=1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    if (annotation.names==TRUE) {
      sector.text=sectors.labels[sector.name]
    } else {
      sector.text=""
    }
    circos.text(
      mean(xlim), ylim[1]+.1, sector.text, 
      facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), 
      col=sectors.colors[sector.name], 
      cex = sector.text.cex
    )
  }, bg.border = NA)
  circos.clear()
  if(any(!is.null(title.main),!is.null(title.sub))) {
    title(
      main=title.main, sub=title.sub,
      cex.main=title.main.cex, cex.sub=title.sub.cex
    )
  }
  if (save.plot==TRUE) 
    dev.off()
  return(dtcall)
}
