
#' XSum extreme sum score
#'
#' @param S_up vector of strings representing the subset of
#' up regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S_down vector of strings representing the subset of
#' down regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the extreme drug perturbation
#' expression data rank-ordered from the highest to the lowest gene scores.
#' @return Xsum
XSum <- function(S_up, S_down, R){
  score_up <- sum(R[match(S_up, names(R))], na.rm = TRUE)
  score_down <- sum(R[match(S_down, names(R))], na.rm = TRUE)
  return(score_up - score_down)
}

#' XSum score 
#' 
#' @description
#' Calculation of the XSum score for all the pertubation DB
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
#' @param mat_R expression matrix or path from which the expression 
#' matrix could be loaded. The matrix must have both colnames and
#' rownames. [test the loading from path]
#' @param list_ranked_R list of rank vectors or path from 
#' which the extreme list of rank vectors could be loaded.
#' The vectors should be named (gene names) so as the list (perturbagene names)
#' [test the loading from path]
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param topN integer to select the N extreme drug
#' perturbation expression data
#' rank-ordered from the highest to the lowest gene scores (eXtreme scores)
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#' @export
XSum_score <- function(S_up, S_down,pert_names, gene_names =NULL, mat_R = NULL, list_ranked_R = NULL,
  topN = NULL, ncpus=1,permuted_pval = FALSE, permute_nb = 10000, padj_method = "BH") {
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
      list_R <- matrixToRankedList(mat_R = readRDS(mat_R), ncpus = ncpus)
    }else{
      list_R <- matrixToRankedList(mat_R = mat_R, ncpus = ncpus)
    }
  }
  if (!is.null(mat_R) & 
    (is.null(colnames(mat_R)) || is.null(rownames(mat_R)))) {
    stop("mat_R should have both rownames and colnames!")
  }
  if (!is.character(S_up)) S_up <- as.character(S_up)
  if (!is.character(S_down)) S_down <- as.character(S_down)

  all_xsum <- unlist(parallel::mclapply(list_R, function(R) {
    return(XSum(S_up = S_up, S_down = S_down, R = R))
  }, mc.cores = ncpus))
  names(all_xsum) <- pert_names

  if(permuted_pval){
    permuteScore <- purrr::map_dfc(1:permute_nb, function(i){
      boot_up <- sample(gene_names, size = length(S_up))
      boot_down <- sample(gene_names, size = length(S_down))
      boot_score <- unlist(parallel::mclapply(list_R, function(R){
        return(XSum(S_up = boot_up, S_down = boot_down, R = R))
      }, mc.cores = ncpus))
      return(boot_score)
    })
    permuteScore[is.na(permuteScore)] <- 0

    pval <- rowSums(sweep(abs(permuteScore), 1,
      abs(all_xsum), ">=")) / permute_nb
    score_df <- data.frame(score = all_xsum, pval = pval)
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(all_xsum)
}

#' Cosine and eXtreme Cosine correlation
#'
#' @param S named vector with disease expression data
#' @param R named vector with drug perturbation expression data
#' @return XCos score (cosine similarity)
XCos <- function(S, R) {
    genes_intersect <- intersect(names(S), names(R))
    if(length(genes_intersect) == 0) return(0)

    S_inter <- S[genes_intersect]
    R_inter <- R[genes_intersect]
    cos_sim <- crossprod(R_inter, S_inter) /
             sqrt(crossprod(R_inter) * crossprod(S_inter))
    return(cos_sim[1, 1])
}

#' XCos score
#' 
#' @description
#' Calculation of the XCos score for all the pertubation DB
#' 
#' @param S named vector with disease expression data
#' @param pert_names names of all the perturbagenes (colums of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`)
#' @param gene_names names of all the genes (rows of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`).
#' Default NULL.
#' @param mat_R expression matrix or path from which the expression 
#' matrix could be loaded. The matrix must have both colnames and
#' rownames. [test the loading from path]
#' @param list_ranked_R list of rank vectors or path from 
#' which the extreme list of rank vectors could be loaded.
#' The vectors should be named (gene names) so as the list (perturbagene names)
#' [test the loading from path]
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param topN integer to select the N extreme drug
#' perturbation expression data
#' rank-ordered from the highest to the lowest gene scores (eXtreme scores)
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#' @export
XCos_score <- function(S, pert_names, gene_names =NULL, mat_R = NULL, list_ranked_R = NULL,
  topN = NULL, ncpus=1,permuted_pval = FALSE, permute_nb = 10000, padj_method = "BH") {
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
      list_R <- matrixToRankedList(mat_R = readRDS(mat_R), ncpus = ncpus)
    }else{
      list_R <- matrixToRankedList(mat_R = mat_R, ncpus = ncpus)
    }
  }
  if (!is.null(mat_R) & 
    (is.null(colnames(mat_R)) || is.null(rownames(mat_R)))) {
    stop("mat_R should have both rownames and colnames!")
  }
  if (is.null(names(S))) stop("S should be a named vector")
  all_xcos <- unlist(parallel::mclapply(list_R, function(R){
    return(XCos(S = S, R = R))
  },mc.cores = ncpus))
  names(all_xcos) <- pert_names

  if(permuted_pval){
    permuteScore <- purrr::map_dfc(1:permute_nb, function(i) {
      names(S) <- sample(gene_names, size = length(S))
      boot_score <- unlist(parallel::mclapply(list_R, function(R) {
        return(XCos(S = S, R = R))
      }, mc.cores = ncpus))
      # names(boot_score) <- pert_names
      return(boot_score)
    })
    permuteScore[is.na(permuteScore)] <- 0

    pval <- rowSums(sweep(abs(permuteScore), 1,
      abs(all_xcos), ">=")) / permute_nb
    score_df <- data.frame(score = all_xcos, pval = pval)
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(all_xcos)
}

#' Pearson and eXtreme Pearson correlation
#'
#' @param S named vector with disease expression data
#' @param R named vector with drug perturbation expression data
#' @param method character to specify the use of `pearson` or `spearman``
#' correlation. (Default `pearson`).
#' @return XCor score (Pearson or Spearman correlation)
XCor <- function(S, R, method = "pearson"){
    genes_intersect <- intersect(names(S), names(R))
    if(length(genes_intersect) == 0) return(0)
    S_inter <- S[genes_intersect]
    R_inter <- R[genes_intersect]
    return(stats::cor(R_inter, S_inter, method = method, use = "na.or.complete"))
}

#' XCor score
#' 
#' @description
#' Calculation of the XCor (Pearson or Spearman) 
#' score for all the pertubation DB
#' 
#' @param S named vector with disease expression data
#' @param pert_names names of all the perturbagenes (colums of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`)
#' @param gene_names names of all the genes (rows of
#' the matrix `mat_R` or names of the elements in the list `list_ranked_R`).
#' Default NULL.
#' @param method character to specify the use of `pearson` or `spearman``
#' correlation. (Default `pearson`).
#' @param mat_R expression matrix or path from which the expression 
#' matrix could be loaded. The matrix must have both colnames and
#' rownames. [test the loading from path]
#' @param list_ranked_R list of rank vectors or path from 
#' which the extreme list of rank vectors could be loaded.
#' The vectors should be named (gene names) so as the list (perturbagene names)
#' [test the loading from path]
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param topN integer to select the N extreme drug
#' perturbation expression data
#' rank-ordered from the highest to the lowest gene scores (eXtreme scores)
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#' @export
XCor_score <- function(S, pert_names, method = "pearson",
  gene_names =NULL, mat_R = NULL, list_ranked_R = NULL,
  topN = NULL, ncpus=1,permuted_pval = FALSE, permute_nb = 10000,
  padj_method = "BH") {
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
      list_R <- matrixToRankedList(mat_R = readRDS(mat_R), ncpus = ncpus)
    }else{
      list_R <- matrixToRankedList(mat_R = mat_R, ncpus = ncpus)
    }
  }
  if (!is.null(mat_R) & 
    (is.null(colnames(mat_R)) || is.null(rownames(mat_R)))) {
    stop("mat_R should have both rownames and colnames!")
  }
  if (is.null(names(S))) stop("S should be a named vector")

  all_xcor <- unlist(parallel::mclapply(list_R, function(R) {
    return(XCor(S = S, R = R,method = method))
  }, mc.cores = ncpus))
  names(all_xcor) <- pert_names

  if(permuted_pval){
    permuteScore <- purrr::map_dfc(1:permute_nb, function(i) {
      names(S) <- sample(gene_names, size = length(S))
      boot_score <- unlist(parallel::mclapply(list_R, function(R) {
        return(XCor(S = S, R = R,method = method))
      }, mc.cores = ncpus))
      return(boot_score)
    })
    permuteScore[is.na(permuteScore)] <- 0

    pval <- rowSums(sweep(abs(permuteScore), 1,
      abs(all_xcor), ">=")) / permute_nb
    score_df <- data.frame(score = all_xcor, pval = pval)
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(all_xcor)
}