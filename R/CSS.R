#' CSS score for ordered disease query
#' @param S named vector with disease expression data 
#' @param R named vector with drug perturbation data representing the drug gene
#' list of interest with expression information.
#'
#' @return CSS score
#' (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-258)
#'
#' @importFrom purrr map_dbl
css_ordered <- function(R, sr_S) {
  Ns <- length(sr_S)
  Nr <- length(R)
  genes_intersect <- intersect(names(R), names(sr_S))
  if(length(genes_intersect) == 0) return(0)
  raw <- sum(R[genes_intersect] * sr_S[genes_intersect], na.rm = T)
  css_max <- sum(purrr::map_dbl(seq(1, Ns), function(i) {
    (Nr - i + 1) * (Ns - i + 1)
  }))

  return(raw / css_max)
}
#' CSS score for unordered disease query
#'
#' @param S unordered named vector with the sign of the disease expression
#' of each gene (-1 or +1).
#' @param R named vector with drug perturbation data representing the drug gene
#' list of interest with expression information.
#'
#' @return CSS score
#' (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-258)
#'
#' @importFrom purrr map_dbl
css_unordered <- function(R, S) {
  Ns <- length(S)
  Nr <- length(R)
  genes_intersect <- intersect(names(R), names(S))
  raw <- sum(R[genes_intersect]*S[genes_intersect], na.rm = T)
  css_max <- sum(purrr::map_dbl(seq(1, Ns), function(i) {
    Nr - i + 1
  }))
  return(raw / css_max)
}

#' CSS score 
#' 
#' @description
#' Calculation of the CSS for all the pertubation DB
#' 
#' @param S_up vector of strings representing the subset of
#' up regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S_down vector of strings representing the subset of
#' down regulated query (ex. disease) genes defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param S named vector with disease expression data.
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
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#' @param unordered_S boolean set to TRUE if S is unordered and just contain
#' gene names
#' @export
css_score <- function(pert_names, S = NULL, S_up = NULL, S_down = NULL, 
  gene_names =NULL, mat_R = NULL,
  list_ranked_R = NULL, ncpus=1,permuted_pval = FALSE,
  permute_nb = 10000, padj_method = "BH", unordered_S = FALSE){
  if (is.null(mat_R) & is.null(list_ranked_R)) {
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
      list_R <- matrixToRankedList(mat_R = readRDS(mat_R), signed_rank = TRUE, ncpus = ncpus)
    }else{
      list_R <- matrixToRankedList(mat_R = mat_R, signed_rank = TRUE, ncpus = ncpus)
    }
  }
  if (!is.null(mat_R) & 
    (is.null(colnames(mat_R)) || is.null(rownames(mat_R)))) {
    stop("mat_R should have both rownames and colnames!")
  }
  if(is.null(S) & is.null(S_up) & is.null(S_down)){
    stop("An input query must be defined either ordered S, or at least
      one of the character vectors S_up or S_down")
  }
  if (!unordered_S & is.null(names(S))) {
    stop("if S is ordered it should be a named vector")
  }

  if(unordered_S) {
    S <- c(rep(c(1, -1), c(length(S_up), length(S_down))))
    names(S) <- c(S_up, S_down)
    all_css <- unlist(parallel::mclapply(list_R, function(R) {
      return(css_unordered(S = S, R = R))
    }, mc.cores = ncpus))
  }else{
    sr_S <- rank(abs(S)) * sign(S)
    all_css <-unlist(parallel::mclapply(list_R, function(R) {
      return(css_ordered(sr_S = sr_S, R = R))
    }, mc.cores = ncpus))
  }
  names(all_css) <- pert_names
  if (permuted_pval) {
    if (unordered_S) {
      permuteScore <- purrr::map_dfc(1:permute_nb, function(i) {
        boot_query <- sample(c(-1,1), replace = TRUE, size = length(S))
        names(boot_query) <- sample(gene_names, size = length(S))
        boot_score <- unlist(parallel::mclapply(list_R, function(R) {
          return(css_unordered(S = boot_query, R = R))
        }, mc.cores = ncpus))
        return(boot_score)
      })
    }else{
      permuteScore <- purrr::map_dfc(1:permute_nb, function(i) {
        names(sr_S) <- sample(gene_names, size = length(sr_S))
        boot_score <- unlist(parallel::mclapply(list_R, function(R) {
          return(css_unordered(S = sr_S, R = R))
        }, mc.cores = ncpus))
        return(boot_score)
      })
    }
    permuteScore[is.na(permuteScore)] <- 0

    pval <- rowSums(sweep(abs(permuteScore), 1,
      abs(all_css), ">=")) / permute_nb
    score_df <- data.frame(score = all_css, pval = pval)
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(all_css)
}
