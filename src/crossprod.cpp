// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


bool comp ( const std::pair<double,int>& a, const std::pair<double,int>& b) {return a.first > b.first;}

std::vector<std::pair<double,int>> index_and_sort(std::vector<double> x, int top_n) {
  // partial sort to get top_n highest scores, but after creating pair to remember the original indices.
  std::vector<std::pair<double,int>> xi(x.size());
  for (int i = 0; i < x.size(); i++) xi[i] = xi[i] = std::pair<double,int>(x[i], i);
  partial_sort(xi.begin(), xi.begin()+top_n, xi.end(), comp);
  return(xi);
}

void fill_triples(std::vector<Eigen::Triplet<double>>& tl, std::vector<double>& res, int i, double min_value, int top_n){
  if (top_n > 0 && top_n < res.size()) {
    std::vector<std::pair<double,int>> res_index = index_and_sort(res, top_n);
    for (int res_i = 0; res_i < top_n; res_i++) {
      if (res_index[res_i].first > min_value)
        tl.push_back(Eigen::Triplet<double>(i, res_index[res_i].second, res_index[res_i].first));
    }

  } else {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      if (res[res_i] > min_value)
        tl.push_back(Eigen::Triplet<double>(i, res_i, res[res_i]));
    }
  }
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> tcrossprod_with_filters_cpp(Eigen::SparseMatrix<double>& m1, Eigen::SparseMatrix<double>& m2,
                                double min_value=0, bool only_upper=false, bool diag=true, int top_n=0, bool verbose=false) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
  m1 = m1.transpose();

  int rows = m1.cols();
  int cols = m2.rows();
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if output is a symmetrical matrix");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if output is a symmetrical matrix");

  std::vector<Eigen::Triplet<double>> tl;
  if (top_n > 0) tl.reserve(rows*top_n); else tl.reserve(rows*cols*0.05);  // for top_n we know max required capacity, otherwise guess (and update in loop)

  Progress p(rows, verbose);
  for (int i = 0; i < rows; i++) {
    std::vector<double> res(cols);
    if (tl.capacity() < tl.size() + cols) {
      double pct_to_go = 1 - i/double(rows);
      tl.reserve(tl.capacity() * (1.2 + 0.8*pct_to_go));
    }
    for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
      for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {

        if (!diag) if (i == it2.row()) continue;
        if (only_upper) if (i > it2.row()) continue;

        res[it2.row()] += it1.value() * it2.value();
      }
    }

    fill_triples(tl, res, i, min_value, top_n);

    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }

  Eigen::SparseMatrix<double> out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}





