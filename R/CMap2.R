#' Enrichment score based on GSEA
#'
#' @description
#' The disease-drug ES in CMap 2.0 is based directly on GSEA’s
#' weighted signed two-sample KS statistic that compares the positions
#' of Sx genes to those of R − SX genes with the weight wES set to 1.
#' @param Sx vector of strings representing the subset of
#' query (ex. disease) genes (S) with the most extreme gene scores
#' (either from the top or bottom) defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the drug perturbation expression data
#' rank-ordered from the highest
#' to the lowest gene scores.
#' @param wES int weight assigned to the positions in R when
#' calculating the ES (default 1).
#' @return enrichment score
ES_gsea <- function(Sx, R, wES) {
  genes_intersect <- sign(match(names(R), Sx, nomatch = 0))
  genes_diff <- 1 - genes_intersect
  Nsx <- length(Sx)
  Nr <- length(R)
  Nm <- Nr - Nsx
  # if (all(genes_intersect == 0)) {
  #   return(-Nr / Nm)
  # }else{
  abs_R <- abs(R)^wES
  Ns_et <- sum(abs_R[genes_intersect == 1])
  es_vect <- cumsum((genes_intersect * abs_R / Ns_et) - (genes_diff / Nm))
  min_ES <- min(es_vect)
  min_ES <- ifelse(is.na(min_ES),0,min_ES)
  max_ES <- max(es_vect)
  max_ES <- ifelse(is.na(max_ES),0,max_ES)
  if(max_ES > - min_ES){
    return(max_ES)
  }else{
    return(min_ES)
  }
  # }
}

#' CMAP2 weighted connectivity score
#'
#' @description
#' WCS ranges from -1 to +1, it is calculated using the enrichment scores for
#' the up and down regulated query lists. 
#' @param S_up vector of strings representing the subset of
#' up regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S_down vector of strings representing the subset of
#' down regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the drug perturbation expression data
#' rank-ordered from the highest to the lowest gene scores.
#' @param wES int weight assigned to the positions in R when
#' calculating the ES (default 1).
#' @return weighted connectivity score
wcs <- function(S_up, S_down, R, wES=1) {
  ES_up <- ES_gsea(S_up, R, wES)
  ES_down <- ES_gsea(S_down, R, wES)
  if (sign(ES_up) != sign(ES_down)) {
    return((ES_up - ES_down) / 2)
  }
  return(0)
}

#' Normalized connectivity score for CMAP2
#' 
#' @description
#' The NCS was developed to enable the comparison of WCS across 
#' cell lines and drug type. Given the WCS for a disease in 
#' relation to a specific drug of a type dt, tested in cell line c.
#' 
#' @param wcs vector of weighted connectivity scores
#' @param touchstone vector of weighted connectivity scores associated
#' to the touchstone signatures.
#' @param cell_line character vector of the cell line in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @param drg_type character vector of the drug type in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @return normalized connectivity score
normalization_cmap2 <- function(wcs, touchstone,
  cell_line = NULL, drg_type = NULL){
  comb <- unique(cbind(cell_line,drg_type))
  if(is.null(comb)){
    wcs[wcs > 0]<- wcs[wcs > 0]/ abs(mean(touchstone[touchstone > 0]))
    wcs[wcs < 0]<- wcs[wcs < 0]/ abs(mean(touchstone[touchstone < 0]))
  }else{
    if(ncol(comb) == 2){
      for(i in 1:nrow(comb)){
        selection<-which(cell_line == comb[i,1] & drg_type==comb[i,2])
        wcs[selection][wcs[selection]>0]<-wcs[selection][wcs[selection]>0]/ abs(mean(wcs[selection][wcs[selection]>0]))
        wcs[selection][wcs[selection]<0]<-wcs[selection][wcs[selection]<0]/ abs(mean(wcs[selection][wcs[selection]<0]))
      }
    }else{
      var <- get(colnames(comb))
      for(i in 1:nrow(comb)){
        selection<-which(var == comb[i,1])
        wcs[selection][wcs[selection]>0]<-wcs[selection][wcs[selection]>0]/ abs(mean(wcs[selection][wcs[selection]>0]))
        wcs[selection][wcs[selection]<0]<-wcs[selection][wcs[selection]<0]/ abs(mean(wcs[selection][wcs[selection]<0]))
      }
    }
  }
  return(wcs)
}
#' CMAP2 connectivity score 
#' 
#' @description
#' Calculation of the weighted or normalized connectivity score for 
#' all the pertubation DB
#' 
#' @param S_up vector of strings representing the subset of
#' up regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S_down vector of strings representing the subset of
#' down regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param pert_names names of all the perturbagenes (colums of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`)
#' @param gene_names names of all the genes (rows of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`).
#' Default NULL.
#' @param cell_line character vector of the cell line in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @param drg_type character vector of the drug type in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @param mat_R expression matrix or path from which the expression 
#' matrix could be loaded. The matrix must have both colnames and
#' rownames. [test the loading from path]
#' @param list_ranked_R list of rank vectors or path from 
#' which the list of rank vectors could be loaded.
#' The vectors should be named (gene names) so as the list (perturbagene names)
#' [test the loading from path]
#' @param wES int weight assigned to the positions in R when
#' calculating the ES (default 1).
#' @param return_weighted boolean set to `TRUE` to return only the weighted
#' connectivity scores. (Default FALSE)
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param touchstone_signatures names of the DB signatures (ex. perturbagenes)
#' considered as touchstones.
#' @param save_wcs boolean set to `TRUE` to save the weighted connectivity score
#' vector as an .RDS in a specific `path`. (Default FALSE)
#' @param path path used to save the weighted connectivity scores if
#' save_wcs is set to TRUE. (Default NULL)
#' @param normalize boolean set to `TRUE` to calculate the 
#' normalized connectivity scores. (Default FALSE)
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#' @export
cmap2_score <- function(S_up, S_down, pert_names, gene_names =NULL,
  cell_line = NULL, drg_type =NULL, 
  mat_R = NULL, list_ranked_R = NULL, wES = 1,
  return_weighted = FALSE, ncpus=1,
  touchstone_signatures=NULL, save_wcs=FALSE, path=NULL, normalize = FALSE,
  permuted_pval = FALSE, permute_nb = 10000, padj_method = "BH") {
  if (is.null(mat_R) & is.null(list_ranked_R)){
    stop("mat_R or list_ranked_R should be defined")
  }
  if (!is.null(list_ranked_R)) {
    if (class(list_ranked_R) == "character") {
      list_R <- readRDS(list_ranked_R)
    }else{
      list_R <- list_ranked_R
    }
  }else{
    if (class(mat_R) == "character") {
      list_R <- matrixToRankedList(mat_R = readRDS(mat_R), expression_val = TRUE, ncpus = ncpus)
    }else{
      list_R <- matrixToRankedList(mat_R = mat_R, expression_val = TRUE, ncpus = ncpus)
    }
  }
  if (!is.null(mat_R) & 
    (is.null(colnames(mat_R)) || is.null(rownames(mat_R)))) {
    stop("mat_R should have both rownames and colnames!")
  }
  if (!is.character(S_up)) S_up <- as.character(S_up)
  if (!is.character(S_down)) S_down <- as.character(S_down)

  all_wcs <- unlist(parallel::mclapply(list_R, function(R) {
      #names(R) <- gene_names
      return(wcs(S_up = S_up, S_down = S_down, R = R, wES = wES))
  }, mc.cores = ncpus))

  names(all_wcs) <- pert_names
  if (return_weighted) {
    return(all_wcs)
  }

  if (save_wcs) {
    saveRDS(all_wcs,file = path)
  }
  if (normalize) {
    if (is.null(touchstone_signatures)) {
      touchstone <- all_wcs
    }else{
      touchstone <- all_wcs[touchstone_signatures]
    }
    norm_wcs <- normalization_cmap2(all_wcs, touchstone, cell_line, drg_type)
    names(norm_wcs) <- pert_names
  }

  if(permuted_pval){
    permuteScore <- purrr::map_dfc(1:permute_nb,function(i){
      boot_up <- sample(gene_names, size = length(S_up))
      boot_down <- sample(gene_names, size = length(S_down))
      boot_score <- unlist(parallel::mclapply(list_R, function(R){
        #names(R) <- gene_names
        return(wcs(S_up = boot_up, S_down = boot_down, R = R, wES = wES))
      },mc.cores = ncpus))
      # names(boot_score) <-pert_names
      if (normalize) {
        if (is.null(touchstone_signatures)) {
          touchstone <- boot_score
        }else{
          touchstone <- boot_score[touchstone_signatures]
        }
        norm_wcs <- normalization_cmap2(boot_score, touchstone, cell_line, drg_type)
        # names(norm_wcs) <- pert_names
        return(norm_wcs)
      }
      return(boot_score)
    })
    permuteScore[is.na(permuteScore)] <- 0
    if(normalize){
      pval <- rowSums(sweep(abs(permuteScore),1,abs(norm_wcs),">="))/permute_nb
      score_df <- data.frame(score = norm_wcs, pval = pval)
    }else{
      pval <- rowSums(sweep(abs(permuteScore),1,abs(all_wcs),">="))/permute_nb
      score_df <- data.frame(score = all_wcs, pval = pval)
    }
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(norm_wcs)
}