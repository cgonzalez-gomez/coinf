#' Enrichment score using the signed one-sample KS test
#' @param Sx vector of strings representing the subset of
#' query (ex. disease) genes (S) with the most extreme gene scores
#' (either from the top or bottom) defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the drug perturbation ranks from the highest
#' to the lowest gene scores.
#' @return enrichment score
# one_sample_KS <- function(Sx, R) {
#   Ns <- length(Sx)
#   Nr <- length(R)
#   
#   m <- match(Sx, names(R))
#   m <- m[!is.na(m)]
#   if(length(m) == 0) return(0)
#   subR <- R[m]
#   subR <- subR[order(subR)]
#   
#   a <- max((1:Ns) / Ns - subR / Nr)
#   b <- max(subR / Nr - ((1:Ns) - 1) / Ns)
#   if (a >= b) {
#     return(a)
#   } else{
#     return(-b)
#   }
# }

#' CMAP1 raw connectivity score
#'
#' @param S_up vector of strings representing the subset of
#' up regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S_down vector of strings representing the subset of
#' down regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the drug perturbation ranks from the highest
#' to the lowest gene scores.
#' @return raw connectivity score
# fcmap1<- function(S_up, S_down, R){
#   ES_up <- one_sample_KS(S_up, R)
#   ES_down <- one_sample_KS(S_down, R)
#   if (sign(ES_up) != sign(ES_down)) {
#     return(ES_up - ES_down)
#   } else {
#     return(0)
#   }
# }

source(file = "~/Documents/tools_bioinfo/coinf/R/CMap1.R")
source(file = "~/Documents/tools_bioinfo/coinf/R/CMap2.R")
source(file = "~/Documents/tools_bioinfo/coinf/R/CSS.R")
source(file = "~/Documents/tools_bioinfo/coinf/R/eXtremeScores.R")

#' @param S query signature (differential expression signature) named vertor.
#' @param Ssig significance values associated to the query expression signature. Default `NULL`
#' @param enrichment_method name of the clusterProfiler function to be used (enrichGO, enrichKEGG, 
#' enricher, GSEA, gseGO, gseKEGG, ...)
#' @param queryCutoff [1] LFC cutoff [2]significance cutoff
#' @param similarity_score default "cmap1" the one used in the original version of GEM2
#' @param strategy default "ORA"
#' @param ... further parameters pass to the enrichment functions in the clusterProfiler package. Exemple 
#' OrgDb="org.Hs.eg.db", ont = "ALL". Keep readable = FALSE for all the ORA (enrichX methods)
#' @import clusterProfiler

GEMS2 <- function(S,list_R,Ssig= NULL, queryCutoff= c(2,0.05), strategy = c("ORA","GSEA"),
                enrichment_method="enrichGO", readable = FALSE,
                pvalueCutoff = 0.01, similarity_score="cmap1",ncores=1,...){
  if(length(strategy)==2) strategy = "ORA"
  S_up <- names(S)[S > queryCutoff[1]]
  S_down <- names(S)[S < (-queryCutoff[1])]
  enrich_meth <- get(enrichment_method)
  if(strategy == "ORA"){
    ego <- enrich_meth(gene = c(S_up, S_down),
                       universe = names(S),
                       readable = readable,
                       ...)
    enriched <- which(ego@result$pvalue < pvalueCutoff)
    enrichedGS <- ego@result$ID[enriched]
  }else if(strategy == "GSEA"){
    ego <- enrich_meth(geneList = sort(S, decreasing = TRUE),
                       pvalueCutoff = pvalueCutoff,...)
    # ego <- enrich_meth(geneList = sort(S,decreasing = TRUE),OrgDb = OrgDb,
    #                    pvalueCutoff = pvalueCutoff)
    enriched <- which(ego@result$pvalue < pvalueCutoff)
    enrichedGS <- ego@result$ID[enriched]
  }
  
  similarity_fn <- switch(similarity_score,
           "cmap1"=cs,
           "cmap2"=wcs,
           "css"=css_ordered,
           "xsum"=XSum,
           "xcos"=XCos,
           "xpears"=XCor,
           "xspear"=function(S, R){
             XCor(S, R, method = "spearman")
           })
  
  
  if(similarity_score %in% c("cmap1","cmap2","xsum")){
    reduced_mat <- purrr::map_dfr(enrichedGS, function(GS){
      up_genes <- intersect(S_up,ego@geneSets[[GS]])
      down_genes <- intersect(S_down,ego@geneSets[[GS]])
      res <- matrix(unlist(parallel::mclapply(list_R,function(R){
        similarity_fn(S_up=up_genes, S_down=down_genes, R=R)
      },mc.cores=ncores)), nrow=1)
      res <- data.frame(res)
      if(is.null(names(list_R))){
        colnames(res) <- seq(1,ncol(res))
      }else{
        colnames(res) <- names(list_R)
      }
      return(res)
    })
  }else if(similarity_score %in% c("xcos","xpears","xspear")){
    reduced_mat <- purrr::map_dfr(enrichedGS, function(GS){
      res <- matrix(unlist(parallel::mclapply(list_R,function(R){
        similarity_fn(S= S[ego@geneSets[[GS]]] , R=R)
      },mc.cores=ncores)), nrow=1)
      res <- data.frame(res)
      if(is.null(names(list_R))){
        colnames(res) <- seq(1,ncol(res))
      }else{
        colnames(res) <- names(list_R)
      }
      return(res)
    })
  }else{ #css
    sr_S <- rank(abs(S)) * sign(S)
    reduced_mat <- purrr::map_dfr(enrichedGS, function(GS){
      res <- matrix(unlist(parallel::mclapply(list_R,function(R){
        similarity_fn(sr_S= sr_S[ego@geneSets[[GS]]], R=R)
      },mc.cores=ncores)), nrow=1)
      res <- data.frame(res)
      if(is.null(names(list_R))){
        colnames(res) <- seq(1,ncol(res))
      }else{
        colnames(res) <- names(list_R)
      }
      return(res)
    })
  }
  rownames(reduced_mat) <- enrichedGS
  colnames(reduced_mat) <- names(list_R)
  return(reduced_mat)
}
