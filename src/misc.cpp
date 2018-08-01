#include <Rcpp.h>
#include <RcppEigen.h>
#include "misc.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

Eigen::SparseMatrix<double> sm_sort(Eigen::SparseMatrix<double>& m, std::vector<std::tuple<double,double,int> > index, bool transpose) {
  if (m.rows() != index.size()) stop("number of rows is not equal to length of index");
  
  std::vector<int> tvec(index.size());
  for (int i = 0; i < index.size(); i++) tvec[std::get<2>(index[i])] = i;
  
  int row, col;
  std::vector<Eigen::Triplet<double>> tl(m.nonZeros());
  for (int k=0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      if (transpose) {
        row = it.col();
        col = tvec[it.row()];
      } else {
        row = tvec[it.row()];
        col = it.col();
      }
      tl.push_back(Eigen::Triplet<double>(row, col, it.value()));
    }
  }

  
  Eigen::SparseMatrix<double> out;
  if (transpose) 
    out = Eigen::SparseMatrix<double>(m.cols(), m.rows());
  else 
    out = Eigen::SparseMatrix<double>(m.rows(), m.cols());
    
  out.setFromTriplets(tl.begin(), tl.end());
  return(out);
}


std::vector<std::tuple<double,double,int> > create_index(Rcpp::IntegerVector group, 
                                                         Rcpp::NumericVector order) {
  std::vector<double> g;
  std::vector<double> o;
  g = as<std::vector<double> >(group);
  o = as<std::vector<double> >(order);
  return(index_and_sort<double,double>(g,o));
}


/*** R
*/
