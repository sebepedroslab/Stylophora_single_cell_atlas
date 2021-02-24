library(topGO)
#library(clusterProfiler)

read_eggnog <- function (file){
  # col_names <- c(
  #   "query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", 
  #   "seed_ortholog_score", "predicted_gene_name", "GO_terms", "KEGG_KOs", 
  #   "BiGG_reactions", "Annotation_tax_scope", "OGs", "bestOG_evalue_score", 
  #   "COG_cat", "eggNOG_annot"
  # )
  col_names <- c(
    "query_name","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score",
    "best_tax_level","Preferred_name","GOs","EC",
    "KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
    "BRITE","KEGG_TC","CAZy","BiGG_Reaction"
  )
  fread(cmd=sprintf("grep -v '^#' %s",file), col.names=col_names, sep="\t", select=1:length(col_names))
}

#' @param list_interest character, gene names
#' @param gomap named list of GO annotations for genes, names of list should be gene names
#' @param output_name character, prefix for output file names
#' @param name_geneset character, used to construct file name for output files
#' @param ontology_set character(s) indicating GO ontology to use, `c("BP","CC","MF")` 
#' @param tg_test character, which test to use, one of `c("fisher","t")`, see `statistic` in `?topGO::runTest`
#' @param tg_algorithm character, which algorithm to use, see `algorithm` in `?topGO::runTest`
#' @param printfile logical, whether to save plot and table
#' @param p_adj_method character, multiple correction method to use
topgofun  <- function(
  list_interest, gomap, output_name, name_geneset, ontology_set, tg_test="fisher", tg_algorithm="classic", 
  topnum=20, nodesize=10, printfile=TRUE, p_adj_method="BH", firstSigNodes=10
) {
  
  library(topGO)
  
  # Input 
  list_interest = unique(list_interest)
  genom = names(gomap)
  gesel = factor(as.integer(genom %in% list_interest))
  names(gesel) = genom
  
  # shortened go mappings without empty transcripts
  gomap_nonempty = gomap[lapply(gomap,length)>0]
  
  namepref <- paste0(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm)
  # if(printfile){
  #   pdf(file=paste0(namepref,".pdf"),height=4.5,width=4)
  # }
  par(mar=c(5,12,5,2))
  
  topgo_tau_tot = data.frame()
  
  if (length(list_interest[list_interest %in% names(gomap_nonempty)])>1) {
    
    for (ontology_seti in ontology_set) {
      # topGO setup 
      
      GOdata = new(
        "topGOdata", ontology=ontology_seti, allGenes=gesel,
        annot=annFUN.gene2GO, gene2GO=gomap
      )
      
      num_interest_feasible = sum(GOdata@feasible & genom %in% list_interest)
      
      # topGO analysis
      topgo_res = runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test)
      topgo_tau = GenTable(
        GOdata, pval_test = topgo_res, orderBy = "pval_test", 
        topNodes = length(usedGO(object = GOdata))
      )
      topGO::printGraph(
        GOdata, result=topgo_res, firstSigNodes=firstSigNodes, # all the nodes in the graph: length(usedGO(object = GOdata)) -- a mess
        useInfo="all", fn.prefix=paste(namepref,ontology_seti,sep="."), pdfSW=TRUE
      )
      topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
      topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
      topgo_tau$ontology = ontology_seti
      topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
      
      # Output 
      # ploti=barplot(height = rev(head(log(topgo_tau$pval_test,10),topnum)),
      #               names.arg = rev(head(paste(topgo_tau$Term,topgo_tau$GO.ID),topnum)),
      #               xlim=c(0,-5),horiz=T,las=1,col="slategray3",border=NA,
      #               cex.names=0.35,cex.axis=0.6,cex.lab=0.6,cex.sub=0.6,cex.main=0.6,
      #               main=paste(name_geneset,"top GO:",ontology_seti,tg_test,tg_algorithm),
      #               sub =paste("n=",num_interest_feasible,"/",length(list_interest), sep=""),
      #               xlab="log(p)")
      # abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
      # abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
      # text(x=0,ploti,labels = paste("p =",signif(rev(head(topgo_tau$pval_test,topnum)),3)),
      #      col="red",pos=4,cex=0.35)
    }
    
  }else {
    print("skip, no annotations in interest list!")
  }
  
  if(printfile){
    write.table(
      topgo_tau_tot,
      file=paste(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm,".txt",sep=""),
      sep="\t", quote=F, col.names=T, row.names=F, append = F)
    dev.off()
  }
  
  return(topgo_tau_tot)
}
