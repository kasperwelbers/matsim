#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


void fill_triples(std::vector<Eigen::Triplet<double> >& tl, std::vector<double>& res, 
                  const std::vector<std::tuple<double,double,int> >& index1, const std::vector<std::tuple<double,double,int> >& index2,
                  int offset, int i, double min_value, int top_n){
  // fill the triples based on the results per row (i), using one of two options: top_n filtering or regular
  // the use_pair makes it somewhat complicated, but is necessary because res is a vector with only the results for pairs used for the current row. 
  // the positions in use_pair are the actual positions, and the true values match the positions in res.
  if (top_n > 0 && top_n < res.size()) {
    std::vector<std::pair<double,int> > res_index = index_and_sort_top_n<double>(res, top_n, offset);
    for (int res_i = 0; res_i < top_n; res_i++) {
      if (res_index[res_i].first > min_value)
        tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_index[res_i].second]), res_index[res_i].first));
    }
    
  } else {
    for (int res_i = 0; res_i < res.size(); res_i++) {
      if (res[res_i] > min_value)
        tl.push_back(Eigen::Triplet<double>(std::get<2>(index1[i]), std::get<2>(index2[res_i+offset]), res[res_i]));
    }
  }
}

std::vector<double> get_colsum(Eigen::SparseMatrix<double>& m) {
  // calculate column sums for an Eigen matrix
  std::vector<double> out(m.outerSize());

  for (int k=0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      out[k] += it.value();
    }
  }
  return(out);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> batched_tcrossprod_cpp(Eigen::SparseMatrix<double>& m1, Eigen::SparseMatrix<double>& m2,
                                                   IntegerVector group1, IntegerVector group2, 
                                                   NumericVector order1, NumericVector order2,
                                                   double min_value=0, int top_n=0, bool diag=true, bool only_upper=false, bool rowsum_div=false, std::string crossfun="prod",
                                                   int lwindow=0, int rwindow=0,
                                                   bool verbose=false, int batchsize = 10000) {
  if (m1.cols() != m2.cols()) stop("m1 and m2 need to have the same number of columns");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  std::vector<std::tuple<double,double,int> > index1, index2;
  index1 = create_index(group1,order1);
  index2 = create_index(group2,order2);
  
  // temp hack to prevent useless batches. Once everything is tested we should make a separate tcrossprod_cpp function without group and order
  bool not_batched = unique(group1).size() == 1 && unique(group2).size() == 1 && unique(order1).size() == 1 && unique(order2).size() == 1;
  if (not_batched) batchsize = m1.rows();
    
  m1 = sm_sort(m1, index1, true);   // sort and transpose m1
  m2 = sm_sort(m2, index2, false);  // only sort m2
  
  std::vector<double> rowsum;
  if (rowsum_div) rowsum = get_colsum(m1);

  int rows = m1.cols();
  int cols = m2.rows();
  if (rows != cols && only_upper) stop("using 'only_upper = true' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  if (rows != cols && !diag) stop("using 'diag = false' is only possible if m1 and m2 have the same number of rows (so output is symmetric)");
  
  bool fun_prod=false, fun_min=false;
  if (crossfun == "prod") fun_prod = true;
  if (crossfun == "min") fun_min = true;
  if (!fun_prod && !fun_min) stop("Not a valid crossfun (currently supports prod and min)");
  
  
  std::vector<Eigen::Triplet<double> > tl;
  if (top_n > 0) {         // for top_n we know max required capacity, otherwise guess (and update in loop)
    tl.reserve(rows*top_n);
  } else {
    tl.reserve(rows*10);  
  }
  
  Progress p(rows, verbose);   
  
  Eigen::SparseMatrix<double> m2_batch;
  
  std::vector<bool> use_pair;
  std::vector<int> positions;
  int i_end;
  int offset = 0;
  
  std::pair<int,int> i_indices;
  std::pair<int,int> batch_indices;
  
  for (int i = 0; i < rows; i++) {
    // CREATE BATCH
    if (i % batchsize == 0) {
      i_end = i + (batchsize - 1);
      if (i_end > rows) i_end = rows - 1;
      batch_indices = find_positions(index2, std::get<0>(index1[i]),
                                     std::get<0>(index1[i_end]),
                                     std::get<1>(index1[i]) - lwindow,
                                     std::get<1>(index1[i_end]) + rwindow);
      m2_batch = m2.middleRows(batch_indices.first, batch_indices.second - batch_indices.first);
      offset = batch_indices.first;
    } 
    std::vector<double> res(m2_batch.rows());
    
    
    // FIND ALL MATCHES FOR i IN BATCH
    std::vector<bool> use_pair(m2_batch.rows());
    double group1_val = std::get<0>(index1[i]);
    double order1_val = std::get<1>(index1[i]);
    double group2_val, order2_val;
    for (int j = 0; j < use_pair.size(); j++) {
      group2_val = std::get<0>(index2[j + offset]);
      if (group2_val != group1_val) continue;
      order2_val = std::get<1>(index2[j + offset]);
      if (order2_val < order1_val-lwindow) continue;
      if (order2_val > order1_val+rwindow) continue;
      if (!diag) if (i == j + offset) continue;
      if (only_upper) if (i > j + offset) continue;
      use_pair[j] = true;
    }
    
    // MANAGE TRIPLET VECTOR CAPACITY
    if (tl.capacity() < tl.size() + res.size()) {
      double pct_to_go = 1 - i/double(rows);
      tl.reserve(tl.capacity() * (1.2 + 0.8*pct_to_go));
    }
    
    // CROSSPROD
    for (Eigen::SparseMatrix<double>::InnerIterator it1(m1,i); it1; ++it1) {
      for (Eigen::SparseMatrix<double>::InnerIterator it2(m2_batch,it1.row()); it2; ++it2) {
        if (!use_pair[it2.row()]) continue;
        if (fun_prod) res[it2.row()] += it1.value() * it2.value();
        if (fun_min) res[it2.row()] += std::min(it1.value(), it2.value());
      }
    }
    if (rowsum_div) {
      for (int res_i = 0; res_i < res.size(); res_i++) {
        if (rowsum[res_i+offset] > 0) res[res_i] = res[res_i] / rowsum[res_i+offset];
      }
    }
    
    // SAVE VALUES WITH CORRECT POSITIONS
    fill_triples(tl, res, index1, index2, offset, i, min_value, top_n);
    if (Progress::check_abort())
      stop("Aborted");
    p.increment(1);
  }
  
  Eigen::SparseMatrix<double> out(rows,cols);
  out.setFromTriplets(tl.begin(), tl.end());
  return out;
}








