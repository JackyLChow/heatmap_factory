library(org.Hs.eg.db)
# input is results of ranked list GSEA
# function to replace ENTREZID with SYMBOL in GSEA outputs
entrezid_to_symbol <- function(gsea_results){
  gsea_out_ <- gsea_results
  gsea_out_[, "core_enrichment_symbol"] <- NULL
  if(nrow(gsea_out_) == 0){
    cat("\nno enriched pathways\n")
    return(gsea_out_)
  } else {
    for(i in rownames(gsea_out_)){
      entrezids_ <- unlist(str_split(gsea_out_[i, "core_enrichment"], "/"))
      symbols_ <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = entrezids_,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID")[, "SYMBOL"]
      gsea_out_[i, "core_enrichment_symbol"] <- paste(symbols_, collapse = "/")
    }
    return(gsea_out_)
  }
}