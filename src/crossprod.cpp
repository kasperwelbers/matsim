// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


bool comp ( const std::pair<double,int>& a, const std::pair<double,int>& b) {return a.first > b.first;}

std::vector<std::pair<double,int> > index_and_sort(std::vector<double>& x, int top_n, std::vector<bool>& use_pair) {
  // partial sort to get top_n highest scores, but after creating pair to remember the original indices.
  std::vector<std::pair<double,int> > xi(x.size());
  int res_i = 0;
  for (int i = 0; i < use_pair.size(); i++) {
    if (use_pair[i]) {
      xi[res_i] = xi[res_i] = std::pair<double,int>(x[res_i], i);
      res_i++;
    }
  }
  partial_sort(xi.begin(), xi.begin()+top_n, xi.end(), comp);
  return(xi);
}

void fill_triples(std::vector<Eigen::Triplet<double> >& tl, std::vector<double>& res, int i, double min_value, int top_n, std::vector<bool>& use_pair){
  if (top_n > 0 && top_n < res.size()) {
    std::vector<std::pair<double,int> > res_index = index_and_sort(res, top_n, use_pair);
    for (int res_i = 0; res_i < top_n; res_i++) {
      if (res_index[res_i].first > min_value)
        tl.push_back(Eigen::Triplet<double>(i, res_index[res_i].second, res_index[res_i].first));
    }

  } else {
    int res_i = 0;
    for (int pair_i = 0; pair_i < use_pair.size(); pair_i++) {
      if (use_pair[pair_i]) {
        if (res[res_i] > min_value)
           tl.push_back(Eigen::Triplet<double>(i, pair_i, res[res_i]));
        res_i++;
      }
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
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");

  std::vector<Eigen::Triplet<double> > tl;
  if (top_n > 0) tl.reserve(rows*top_n); else tl.reserve(rows*cols*0.05);  // for top_n we know max required capacity, otherwise guess (and update in loop)

  std::vector<bool> use_pair(cols);
  Progress p(rows, verbose);    // removed from cran
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

    fill_triples(tl, res, i, min_value, top_n, use_pair);

    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }

  Eigen::SparseMatrix<double> out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> tcrossprod_cpp(Eigen::SparseMatrix<double>& m1, Eigen::SparseMatrix<double>& m2,
                                           double min_value=0, int top_n=0, bool diag=true, bool only_upper=false,
                                           Rcpp::Nullable<Rcpp::StringVector> group1 = R_NilValue, Rcpp::Nullable<Rcpp::StringVector> group2 = R_NilValue, 
                                           Rcpp::Nullable<Rcpp::NumericVector> order1 = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> order2 = R_NilValue, int lwindow=0, int rwindow=0,
                                           bool verbose=false) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
  m1 = m1.transpose();
  
  int rows = m1.cols();
  int cols = m2.rows();
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  
  
  bool groupfilter = false;
  bool orderfilter = false;
  Rcpp::StringVector g1_s, g2_s;
  Rcpp::IntegerVector g1, g2;
  Rcpp::NumericVector o1, o2;
  if (group1.isNotNull() && group2.isNotNull()) {
     g1_s = Rcpp::StringVector(group1);
     g2_s = Rcpp::StringVector(group2);
     g1 = match(g1_s, sort_unique(g1_s));
     g2 = match(g2_s, sort_unique(g2_s));
     groupfilter = g1.size() == rows && g2.size() == cols;
  }
  if (order1.isNotNull() && order2.isNotNull()) {
    o1 = Rcpp::NumericVector(order1);
    o2 = Rcpp::NumericVector(order2);
    orderfilter = o1.size() == rows && o2.size() == cols;
  }
  

  std::vector<Eigen::Triplet<double> > tl;
  if (top_n > 0) {         // for top_n we know max required capacity, otherwise guess (and update in loop)
    tl.reserve(rows*top_n);
  } else {
    tl.reserve(rows*50);  
  }
  
  Progress p(rows, verbose);   
  std::vector<int> positions(cols);
  
  for (int i = 0; i < rows; i++) {
    int position_i = 0;
    std::vector<bool> use_pair(cols);
    
    int group_i;
    double order_i_l, order_i_r;
    if (groupfilter) group_i = g1[i];
    if (orderfilter) {
      order_i_l = o1[i] - lwindow;
      order_i_r = o1[i] + rwindow;
    }
    
    // make function to get j as vector based on binary search in group and order
    for (int j = 0; j < cols; j++) {
      if (groupfilter) if (group_i != g2[j]) continue;
      if (orderfilter) {
        if (order_i_l > o2[j]) continue; 
        if (order_i_r < o2[j]) continue;
      }
      if (!diag) if (i == j) continue;
      if (only_upper) if (i > j) continue;
      use_pair[j] = true;
      positions[j] = position_i;
      position_i++;
    }

    std::vector<double> res(position_i);
    if (tl.capacity() < tl.size() + position_i) {
      double pct_to_go = 1 - i/double(rows);
      tl.reserve(tl.capacity() * (1.2 + 0.8*pct_to_go));
    }
    
    for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
      for (Eigen::SparseMatrix<double>::InnerIterator it2(m2,it1.row()); it2; ++it2) {
        if (!use_pair[it2.row()]) continue;
        res[positions[it2.row()]] += it1.value() * it2.value();
      }
    }
    
    fill_triples(tl, res, i, min_value, top_n, use_pair);
    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }
  
  Eigen::SparseMatrix<double> out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}




/*** R
#m = Matrix::rsparsematrix(1000,1000,0.01)
#m2 = t(m)
*/
