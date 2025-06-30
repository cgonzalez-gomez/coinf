#' Transform an expression data matrix to a ranked list
#' 
#' @description 
#' drug perturbation ranks from the highest
#' to the lowest gene scores.
#' @param mat_R expression data matrix with colnames (perturbagene)
#' and colnames (genes).
#' @param ncpus number of cores used in the parallel calculations (Default 1).
#' @param expression_val boolean set to `TRUE` to return the drug
#' perturbation expression data rank-ordered from the highest to the lowest gene scores (CMAP2 and eXtremeScores).
#' @param signed_rank boolean set to `TRUE` to return the drug perturbation
#' signed rank from the highest to the lowest absolute gene scores. (CSS)
#' @export 
matrixToRankedList <- function(mat_R, ncpus = 1, expression_val = FALSE,
  signed_rank = FALSE) {
  ## Allocate memory for the refList
  if(is.null(rownames(mat_R))){
    stop("mat_R must have rownames")
  }
  rnames <- rownames(mat_R)
  if (signed_rank){
    list_R <- parallel::mclapply(mat_R, function(x){
      names(x) <- rnames
      x <- x[!is.na(x)]
      sorted <- x[order(x, decreasing = T)]
      return(rank(abs(sorted)) * sign(sorted))
    }, mc.cores = ncpus)
  }else if (expression_val){
    list_R <- parallel::mclapply(mat_R, function(x) {
      names(x) <- rnames
      x <- x[!is.na(x)]
      return(x[order(x, decreasing = T)])
    }, mc.cores = ncpus)
  }else{
    list_R <- parallel::mclapply(mat_R, function(x) {
      names(x) <- rnames
      x <- x[!is.na(x)]
      rank(-x)
    }, mc.cores = ncpus)
  }
  return(list_R)
}