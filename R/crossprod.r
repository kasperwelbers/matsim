
#' The tcrossprod function for sparse matrices, with output filters applied on the fly to reduce memory usage.
#'
#' @param m1 A dgCMatrix
#' @param m2 A dgCMatrix
#' @param min_value a numerical value, specifying the threshold for including a score in the output
#' @param only_upper if true, only the upper triangle of the matrix is returned. Only possible for symmetrical output (m1 and m2 have same number of columns)
#' @param diag if false, the diagonal of the matrix is not returned. Only possible for symmetrical output (m1 and m2 have same number of columns)
#' @param top_n an integer, specifying the top number of strongest scores for each column in m1
#' @param verbose if TRUE, report progress
#'
#' @return A dgCMatrix
#' @export
#'
#' @examples
#' set.seed(1)
#' m = Matrix::rsparsematrix(5,10,0.5)
#' tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = T)
#' tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = F)
#' tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F)
#' tcrossprod_sparse(m, min_value = 0.2, only_upper = T, diag = F)
#' tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F, top_n = 1)
tcrossprod_sparse <- function(m1, m2=NULL, min_value=0, only_upper=F, diag=T, top_n=NULL, verbose=F) {
  if (is.null(top_n)) top_n = 0
  if (is.null(m2)) m2 = m1
  tcrossprod_with_filters_cpp(m1, m2, min_value, only_upper, diag, top_n, verbose=verbose)
}

