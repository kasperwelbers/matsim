#' The tcrossprod function for sparse matrices, with output filters applied on the fly to reduce memory usage.
#'
#' @param m A dgCMatrix
#' @param m2 A dgCMatrix
#' @param min_value a numerical value, specifying the threshold for including a score in the output
#' @param only_upper if true, only the upper triangle of the matrix is returned. Only possible for symmetrical output (m and m2 have same number of columns)
#' @param diag if false, the diagonal of the matrix is not returned. Only possible for symmetrical output (m and m2 have same number of columns)
#' @param top_n an integer, specifying the top number of strongest scores for each column in m
#' @param rowsum_div if true, divide crossproduct by column sums of m. (this has to happen within the loop for min_value and top_n filtering)
#' @param crossfun The function used in the vector operations. Normally this is the "prod", for product (dot product). Here we also allow (currently only) the "min", for minimum value. We use this in our document overlap_pct score.
#' @param group Optionally, a character vector that specifies a group (e.g., source) for each row in m. If given, only pairs of rows with the same group are calculated. 
#' @param group2 If m2 and group are used, group2 has to be used to specify the groups for the rows in m2 (otherwise group will be ignored)
#' @param date Optionally, a character vector that specifies a date for each row in m. If given, only pairs of rows within a given date range (see lwindow, rwindow and date_unit) are calculated. 
#' @param date2 If m2 and date are used, date2 has to be used to specify the date for the rows in m2 (otherwise date will be ignored)
#' @param lwindow If date (and date2) are used, lwindow determines the left side of the date window. Also see date_unit.
#' @param rwindow Like lwindow, but for the right side. 
#' @param date_unit The date unit used in lwindow and rwindow. Supports "days", "hours", "minutes" and "seconds". Note that refers to the time distance between two rows ("days" doesn't refer to calendar days, but to a time of 24 hours)
#' @param batchsize experimental 
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
tcrossprod_sparse <- function(m, m2=NULL, min_value=0, only_upper=F, diag=T, top_n=NULL, rowsum_div=F, crossfun='prod', group=NULL, group2=NULL, date=NULL, date2=NULL, lwindow=1, rwindow=1, date_unit=c('days','hours','minutes','seconds'), batchsize=10000, verbose=F) {
  date_unit = match.arg(date_unit)
  if (is.null(top_n)) top_n = 0
  if (is.null(m2)) {
    m2 = m
    group2 = group
    date2 = date
  }
  if (!is.null(group) && !is.null(group2)) {
    if (!length(group) == nrow(m)) stop('group has to have the same length as the number of rows in m')
    if (!length(group2) == nrow(m2)) stop('group2 has to have the same length as the number of rows in m2')
    group = as.character(group)
    group2 = as.character(group2)
    unique_group = unique(c(group,group2))
    group = match(group, unique_group)
    group2 = match(group2, unique_group)
  } else {
    group = rep(1, nrow(m))
    group2 = rep(1, nrow(m2))
  }
  
  if (!is.null(date) && !is.null(date2)) {
    if (!length(date) == nrow(m)) stop('date has to have the same length as the number of rows in m')
    if (!length(date2) == nrow(m2)) stop('date2 has to have the same length as the number of rows in m2')
    
    date = as.POSIXct(date)
    date2 = as.POSIXct(date2)
    order1 = as.numeric(date)
    order2 = as.numeric(date2)
    startorder = min(c(order1, order2))
    order1 = order1 - startorder
    order2 = order2 - startorder
    if (date_unit == 'seconds') unit_multip = 1
    if (date_unit == 'minutes') unit_multip = 60
    if (date_unit == 'hours') unit_multip = 60 * 60    
    if (date_unit == 'days') unit_multip = 60 * 60 * 24
    lwindow = lwindow * unit_multip
    rwindow = rwindow * unit_multip

  } else {
    order1 = rep(1, nrow(m))
    order2 = rep(1, nrow(m2))
  }
  
  batched_tcrossprod_cpp(m, m2, group1=group, group2=group2, order1=order1, order2=order2, min_value=min_value, top_n=top_n, diag=diag, only_upper=only_upper, rowsum_div=rowsum_div, crossfun=crossfun,
                 lwindow=lwindow, rwindow=rwindow, verbose=verbose, batchsize=batchsize)
}

function(){
set.seed(1)
m = Matrix::rsparsematrix(10,10,0.5)
tcrossprod_sparse(m)
tcrossprod_sparse(m, group = c(1,1,1,2,2,2,3,3,3,3))
tcrossprod_sparse(m, date = seq.Date(as.Date('2010-01-01'), as.Date('2010-02-01'), by = 11)[c(1,1,1,2,2,2,3,3,3,3)])
}