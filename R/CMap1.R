#' Enrichment score using the signed one-sample KS test
#'
#' @description
#' CMap 1.0 uses a signed one-sample KS test to compare the
#' empirical distribution of the positions of SX genes in R
#' compared to a reference uniform distribution (of disease
#' genes in the drug gene list)
#' @param Sx vector of strings representing the subset of
#' query (ex. disease) genes (S) with the most extreme gene scores
#' (either from the top or bottom) defined based on a
#' user-specified threshold of fold-change and/or significance.
#' @param R named vector with the drug perturbation ranks from the highest
#' to the lowest gene scores.
#' @return enrichment score
ES_cmap1 <- function(Sx, R) {
  Ns <- length(Sx)
  Nr <- length(R)
  
  m <- match(Sx, names(R))
  m <- m[!is.na(m)]
  if(length(m) == 0) return(0)
  subR <- R[m]
  subR <- subR[order(subR)]

  a <- max((1:Ns) / Ns - subR / Nr)
  b <- max(subR / Nr - ((1:Ns) - 1) / Ns)
  if (a >= b) {
    return(a)
  } else{
    return(-b)
  }
}

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
cs<- function(S_up, S_down, R){
  ES_up <- ES_cmap1(S_up, R)
  ES_down <- ES_cmap1(S_down, R)
  if (sign(ES_up) != sign(ES_down)) {
      return(ES_up - ES_down)
  } else {
      return(0)
  }
}

#' Normalized connectivity score for CMAP1. This normalization is
#' needed for CMAP1.
#'
#' @description
#' The final connectivity score is calculated by normalizing
#' the raw score by dividing by the maximum or minimum of raw
#' scores across treatment instances, depending on the sign of cs,
#' bringing it back to range between −1 and+ 1.
#' The use of `touchstone`, `cell_line` and `drg_type` parameters were 
#' inspired from CMAP2 score, not necessary used in the original CMAP1 method.
#'
#' @param cs vector of raw connectivity scores
#' @param touchstone vector of raw connectivity scores associated
#' to the touchstone signatures.
#' @param cell_line character vector of the cell line in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @param drg_type character vector of the drug type in which the
#' signature (ex. molecule) from the DB was generated. (Default NULL)
#' @return normalized connectivity score
normalization_cmap1 <- function(cs, touchstone,
  cell_line = NULL, drg_type = NULL) {
  comb <- unique(cbind(cell_line, drg_type))
  if(is.null(comb)) {
    max_cs <- max(touchstone)
    min_cs <- min(touchstone)
    cs[cs > 0] <- cs[cs > 0] / max_cs
    cs[cs < 0] <- -cs[cs < 0] / min_cs
  }else{
    if(ncol(comb) == 2) {
      for(i in 1:nrow(comb)) {
        selection <- which(cell_line == comb[i, 1] & drg_type == comb[i, 2])
        max_cs <- max(cs[selection])
        min_cs <- min(cs[selection])
        cs[selection][cs[selection] > 0] <-
          cs[selection][cs[selection] > 0] / max_cs
        cs[selection][cs[selection] < 0] <-
          -cs[selection][cs[selection] < 0] / min_cs
      }
    }else{
      var <- get(colnames(comb))
      for (i in 1:nrow(comb)) {
        selection <- which(var == comb[i, 1])
        max_cs <- max(cs[selection])
        min_cs <- min(cs[selection])
        cs[selection][cs[selection] > 0] <-
          cs[selection][cs[selection] > 0] / max_cs
        cs[selection][cs[selection] < 0] <-
          -cs[selection][cs[selection] < 0] / min_cs
      }
    }
  }
  return(cs)
}

#' CMAP1 connectivity score 
#' 
#' @description
#' Calculation of the raw or normalized connectivity score for 
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
#' @param return_raw boolean set to `TRUE` to return only the raw
#' connectivity scores. (Default FALSE)
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param touchstone_signatures names of the DB signatures (ex. perturbagenes)
#' considered as touchstones.
#' @param save_raw boolean set to `TRUE` to save the raw connectivity score
#' vector as an .RDS in a specific `path`. (Default FALSE)
#' @param path path used to save the raw connectivity scores if
#' save_raw is set to TRUE. (Default NULL)
#' @param normalize boolean set to `TRUE` to calculate the 
#' normalized connectivity scores. (Default FALSE)
#' @param permuted_pval boolean set to `TRUE` to calculate the permuted pvalue
#' of the scores. (Default FALSE)
#' @param permute_nb number of permutation to calcule the pvalue
#' (Default 10000).
#' @param padj_method Method of adjustment of the pvalue. See stats::p.adjust
#' `method` parameter. (Default 'BH')
#'@export 
cmap1_score <- function(S_up, S_down, gene_names, pert_names, mat_R = NULL,
  cell_line = NULL, drg_type = NULL,
  list_ranked_R = NULL, return_raw = FALSE, ncpus=1,
  touchstone_signatures = NULL, save_raw = FALSE, path=NULL, normalize = FALSE,
  permuted_pval = FALSE, permute_nb = 10000, padj_method = "BH"){
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

  all_cs <- unlist(parallel::mclapply(list_R, function(R) {
    tryCatch(
        expr = {
            return(cs(S_up = S_up, S_down = S_down, R = R))
        },
        error = function(e) {
            return(e$message)
        }
    )
}, mc.cores = ncpus))
  names(all_cs) <- pert_names

  if (return_raw) {
    return(all_cs)
  }

  if (save_raw) {
    saveRDS(all_cs,file = path)
  }

  if (normalize) {
    if (is.null(touchstone_signatures)) {
      touchstone <- all_cs
    }else{
      touchstone <- all_cs[touchstone_signatures]
    }
    norm_cs <- normalization_cmap1(all_cs, touchstone, cell_line, drg_type)
    names(norm_cs) <- pert_names
  }

  if (permuted_pval) {
    permuteScore <- purrr::map_dfc(1:permute_nb, function(i) {
      boot_up <- sample(gene_names, size = length(S_up))
      boot_down <- sample(gene_names, size = length(S_down))
      boot_score <- unlist(parallel::mclapply(list_R, function(R) {
        return(cs(S_up = boot_up, S_down = boot_down, R = R))
      },mc.cores = ncpus))
      if (normalize) {
        if (is.null(touchstone_signatures)) {
          touchstone <- boot_score
        }else{
          touchstone <- boot_score[touchstone_signatures]
        }
        norm_cs <- normalization_cmap1(boot_score,
          touchstone, cell_line, drg_type)
        return(norm_cs)
      }
      return(boot_score)
    })
    permuteScore[is.na(permuteScore)] <- 0
    if (normalize) {
      pval <- rowSums(sweep(abs(permuteScore),1,abs(norm_cs),">="))/permute_nb
      score_df <- data.frame(score = norm_cs, pval = pval)
    }else{
      pval <- rowSums(sweep(abs(permuteScore),1,abs(all_cs),">="))/permute_nb
      score_df <- data.frame(score = all_cs, pval = pval)
    }
    score_df$padj <- stats::p.adjust(pval, method = padj_method)
    rownames(score_df) <- pert_names
    return(score_df)
  }
  return(norm_cs)
}