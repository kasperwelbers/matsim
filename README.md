Example
=======

    devtools::install_github('kasperwelbers/matsim')

    library(matsim)

    ## Loading required package: Matrix

    library(Matrix)
    library(microbenchmark)

Several options for filtering results. Filtering is performed on the
fly, after each iteration (rows in m).

    set.seed(1)
    m = Matrix::rsparsematrix(5,10,0.5)
    tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = T)

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                                            
    ## [1,] 0.1936 .        0.1100 0.2200 .       
    ## [2,] .      2.010149 .      .      .       
    ## [3,] 0.1100 .        2.7628 0.1657 .       
    ## [4,] 0.2200 .        0.1657 4.2504 .       
    ## [5,] .      .        .      .      2.189209

    tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = F)

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                            
    ## [1,] .    . 0.1100 0.2200 .
    ## [2,] .    . .      .      .
    ## [3,] 0.11 . .      0.1657 .
    ## [4,] 0.22 . 0.1657 .      .
    ## [5,] .    . .      .      .

    tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F)

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                       
    ## [1,] . . 0.11 0.2200 .
    ## [2,] . . .    .      .
    ## [3,] . . .    0.1657 .
    ## [4,] . . .    .      .
    ## [5,] . . .    .      .

    tcrossprod_sparse(m, min_value = 0.2, only_upper = T, diag = F)

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                  
    ## [1,] . . . 0.22 .
    ## [2,] . . . .    .
    ## [3,] . . . .    .
    ## [4,] . . . .    .
    ## [5,] . . . .    .

    tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F, top_n = 1)

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                    
    ## [1,] . . . 0.2200 .
    ## [2,] . . . .      .
    ## [3,] . . . 0.1657 .
    ## [4,] . . . .      .
    ## [5,] . . . .      .

    ## filtering by group/date. Documents have to be in the same group, or within the given date range
    m = Matrix::rsparsematrix(10,10,0.5)
    tcrossprod_sparse(m, group = c(1,1,1,2,2,2,3,3,3,3))

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                     
    ##  [1,] 3.029609 0.3356 .      .     .        .      .        .       
    ##  [2,] 0.335600 2.7838 .      .     .        .      .        .       
    ##  [3,] .        .      2.9841 .     .        .      .        .       
    ##  [4,] .        .      .      4.574 0.400000 1.0920 .        .       
    ##  [5,] .        .      .      0.400 2.984471 .      .        .       
    ##  [6,] .        .      .      1.092 .        4.8973 .        .       
    ##  [7,] .        .      .      .     .        .      1.594825 .       
    ##  [8,] .        .      .      .     .        .      .        0.000986
    ##  [9,] .        .      .      .     .        .      .        0.004560
    ## [10,] .        .      .      .     .        .      1.469000 0.028500
    ##                     
    ##  [1,] .       .     
    ##  [2,] .       .     
    ##  [3,] .       .     
    ##  [4,] .       .     
    ##  [5,] .       .     
    ##  [6,] .       .     
    ##  [7,] .       1.4690
    ##  [8,] 0.00456 0.0285
    ##  [9,] 4.00740 0.1076
    ## [10,] 0.10760 6.1493

    date = seq.Date(as.Date('2010-01-01'), as.Date('2010-01-10'), by=1)
    tcrossprod_sparse(m, date = date, lwindow = 1, rwindow = 1)

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                     
    ##  [1,] 3.029609 0.3356 .      .     .        .      .        .       
    ##  [2,] 0.335600 2.7838 .      .     .        .      .        .       
    ##  [3,] .        .      2.9841 1.289 .        .      .        .       
    ##  [4,] .        .      1.2890 4.574 0.400000 .      .        .       
    ##  [5,] .        .      .      0.400 2.984471 .      .        .       
    ##  [6,] .        .      .      .     .        4.8973 .        .       
    ##  [7,] .        .      .      .     .        .      1.594825 .       
    ##  [8,] .        .      .      .     .        .      .        0.000986
    ##  [9,] .        .      .      .     .        .      .        0.004560
    ## [10,] .        .      .      .     .        .      .        .       
    ##                     
    ##  [1,] .       .     
    ##  [2,] .       .     
    ##  [3,] .       .     
    ##  [4,] .       .     
    ##  [5,] .       .     
    ##  [6,] .       .     
    ##  [7,] .       .     
    ##  [8,] 0.00456 .     
    ##  [9,] 4.00740 0.1076
    ## [10,] 0.10760 6.1493

    tcrossprod_sparse(m, date = date, lwindow = 2, rwindow = 2)

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                      
    ##  [1,] 3.029609 0.3356 .       .     .        .      .        .       
    ##  [2,] 0.335600 2.7838 .       .     .        .      .        .       
    ##  [3,] .        .      2.98410 1.289 0.096360 .      .        .       
    ##  [4,] .        .      1.28900 4.574 0.400000 1.0920 .        .       
    ##  [5,] .        .      0.09636 0.400 2.984471 .      0.213400 .       
    ##  [6,] .        .      .       1.092 .        4.8973 .        .       
    ##  [7,] .        .      .       .     0.213400 .      1.594825 .       
    ##  [8,] .        .      .       .     .        .      .        0.000986
    ##  [9,] .        .      .       .     .        .      .        0.004560
    ## [10,] .        .      .       .     .        .      .        0.028500
    ##                     
    ##  [1,] .       .     
    ##  [2,] .       .     
    ##  [3,] .       .     
    ##  [4,] .       .     
    ##  [5,] .       .     
    ##  [6,] .       .     
    ##  [7,] .       .     
    ##  [8,] 0.00456 0.0285
    ##  [9,] 4.00740 0.1076
    ## [10,] 0.10760 6.1493

    ## using conditional probability
    m = abs(Matrix::rsparsematrix(10,10,0.5))
    m[m>0] = 1

    t(t(tcrossprod_sparse(m)) / rowSums(m))

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                    
    ##  [1,] 1.0000000 0.6666667 0.7142857 0.4 0.3333333 0.6 1.0000000 0.4
    ##  [2,] 0.6666667 1.0000000 0.7142857 0.4 .         0.4 0.6666667 0.6
    ##  [3,] 0.8333333 0.8333333 1.0000000 0.8 0.3333333 0.6 0.6666667 0.6
    ##  [4,] 0.3333333 0.3333333 0.5714286 1.0 0.6666667 0.4 .         0.6
    ##  [5,] 0.1666667 .         0.1428571 0.4 1.0000000 0.6 0.3333333 0.4
    ##  [6,] 0.5000000 0.3333333 0.4285714 0.4 1.0000000 1.0 1.0000000 0.6
    ##  [7,] 0.5000000 0.3333333 0.2857143 .   0.3333333 0.6 1.0000000 0.2
    ##  [8,] 0.3333333 0.5000000 0.4285714 0.6 0.6666667 0.6 0.3333333 1.0
    ##  [9,] 0.5000000 0.6666667 0.7142857 1.0 0.6666667 0.6 0.3333333 0.8
    ## [10,] .         0.3333333 0.1428571 0.4 0.3333333 0.2 .         0.4
    ##                          
    ##  [1,] 0.4285714 .        
    ##  [2,] 0.5714286 0.6666667
    ##  [3,] 0.7142857 0.3333333
    ##  [4,] 0.7142857 0.6666667
    ##  [5,] 0.2857143 0.3333333
    ##  [6,] 0.4285714 0.3333333
    ##  [7,] 0.1428571 .        
    ##  [8,] 0.5714286 0.6666667
    ##  [9,] 1.0000000 1.0000000
    ## [10,] 0.4285714 1.0000000

    tcrossprod_sparse(m, rowsum_div = T)

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                    
    ##  [1,] 1.0000000 0.6666667 0.7142857 0.4 0.3333333 0.6 1.0000000 0.4
    ##  [2,] 0.6666667 1.0000000 0.7142857 0.4 .         0.4 0.6666667 0.6
    ##  [3,] 0.8333333 0.8333333 1.0000000 0.8 0.3333333 0.6 0.6666667 0.6
    ##  [4,] 0.3333333 0.3333333 0.5714286 1.0 0.6666667 0.4 .         0.6
    ##  [5,] 0.1666667 .         0.1428571 0.4 1.0000000 0.6 0.3333333 0.4
    ##  [6,] 0.5000000 0.3333333 0.4285714 0.4 1.0000000 1.0 1.0000000 0.6
    ##  [7,] 0.5000000 0.3333333 0.2857143 .   0.3333333 0.6 1.0000000 0.2
    ##  [8,] 0.3333333 0.5000000 0.4285714 0.6 0.6666667 0.6 0.3333333 1.0
    ##  [9,] 0.5000000 0.6666667 0.7142857 1.0 0.6666667 0.6 0.3333333 0.8
    ## [10,] .         0.3333333 0.1428571 0.4 0.3333333 0.2 .         0.4
    ##                          
    ##  [1,] 0.4285714 .        
    ##  [2,] 0.5714286 0.6666667
    ##  [3,] 0.7142857 0.3333333
    ##  [4,] 0.7142857 0.6666667
    ##  [5,] 0.2857143 0.3333333
    ##  [6,] 0.4285714 0.3333333
    ##  [7,] 0.1428571 .        
    ##  [8,] 0.5714286 0.6666667
    ##  [9,] 1.0000000 1.0000000
    ## [10,] 0.4285714 1.0000000

    ## using the min function
    m = abs(Matrix::rsparsematrix(10,10,0.5))

    tcrossprod_sparse(m, crossfun = 'min')

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                      
    ##  [1,] 4.4956 0.352 0.2456 2.4656 0.490 0.540 2.192 1.1906 .     2.450
    ##  [2,] 0.3520 2.292 1.3450 0.7920 1.380 0.470 0.244 2.1920 0.280 0.942
    ##  [3,] 0.2456 1.345 3.6030 1.5730 1.973 0.153 0.240 2.5080 0.043 0.783
    ##  [4,] 2.4656 0.792 1.5730 7.3200 1.140 2.450 2.002 3.3450 0.890 2.893
    ##  [5,] 0.4900 1.380 1.9730 1.1400 3.190 0.610 0.490 1.5000 0.500 1.490
    ##  [6,] 0.5400 0.470 0.1530 2.4500 0.610 4.960 0.082 3.1000 1.900 1.083
    ##  [7,] 2.1920 0.244 0.2400 2.0020 0.490 0.082 2.192 0.1770 .     2.110
    ##  [8,] 1.1906 2.192 2.5080 3.3450 1.500 3.100 0.177 8.5450 1.500 1.735
    ##  [9,] .      0.280 0.0430 0.8900 0.500 1.900 .     1.5000 1.900 1.000
    ## [10,] 2.4500 0.942 0.7830 2.8930 1.490 1.083 2.110 1.7350 1.000 4.133

    tcrossprod_sparse(m, crossfun = 'min', rowsum_div = T) ## what percentage of 

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                     
    ##  [1,] 1.00000000 0.1535777 0.06816542 0.3368306 0.1536050 0.10887097
    ##  [2,] 0.07829878 1.0000000 0.37330003 0.1081967 0.4326019 0.09475806
    ##  [3,] 0.05463119 0.5868237 1.00000000 0.2148907 0.6184953 0.03084677
    ##  [4,] 0.54844737 0.3455497 0.43658063 1.0000000 0.3573668 0.49395161
    ##  [5,] 0.10899546 0.6020942 0.54759922 0.1557377 1.0000000 0.12298387
    ##  [6,] 0.12011745 0.2050611 0.04246461 0.3346995 0.1912226 1.00000000
    ##  [7,] 0.48758786 0.1064572 0.06661116 0.2734973 0.1536050 0.01653226
    ##  [8,] 0.26483673 0.9563700 0.69608659 0.4569672 0.4702194 0.62500000
    ##  [9,] .          0.1221640 0.01193450 0.1215847 0.1567398 0.38306452
    ## [10,] 0.54497731 0.4109948 0.21731890 0.3952186 0.4670846 0.21834677
    ##                                                 
    ##  [1,] 1.00000000 0.13933294 .          0.5927897
    ##  [2,] 0.11131387 0.25652428 0.14736842 0.2279216
    ##  [3,] 0.10948905 0.29350497 0.02263158 0.1894508
    ##  [4,] 0.91332117 0.39145699 0.46842105 0.6999758
    ##  [5,] 0.22354015 0.17554125 0.26315789 0.3605129
    ##  [6,] 0.03740876 0.36278525 1.00000000 0.2620373
    ##  [7,] 1.00000000 0.02071387 .          0.5105250
    ##  [8,] 0.08074818 1.00000000 0.78947368 0.4197919
    ##  [9,] .          0.17554125 1.00000000 0.2419550
    ## [10,] 0.96259124 0.20304272 0.52631579 1.0000000

Speed seems to be comparable to Matrix::crossprod(), with some speed-ups
if filtering is used. (This needs to be more properly tested)

    m = abs(Matrix::rsparsematrix(3000,10000,0.1))

    ## quick and dirty L2 norm (for cosine similarity)
    m = as(m, 'dgTMatrix')
    norm = sqrt(Matrix::colSums(m^2))
    m@x = m@x / norm[m@j+1]
    m = as(m, 'dgCMatrix')

    cp1 = tcrossprod_sparse(m)
    cp2 = tcrossprod(m)
    identical(as(cp1, 'dgCMatrix'), as(cp2, 'dgCMatrix'))

    ## [1] TRUE

    ## regular
    microbenchmark(tcrossprod_sparse(m),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                  expr      min       lq     mean   median       uq
    ##  tcrossprod_sparse(m) 2.446588 2.447164 2.462931 2.455663 2.482560
    ##       tcrossprod((m)) 2.247728 2.306535 2.316787 2.317991 2.342892
    ##       max neval cld
    ##  2.482679     5   b
    ##  2.368787     5  a

    ## only calculate upper triangle and do not calculate diagonal
    microbenchmark(tcrossprod_sparse(m, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                            expr      min       lq     mean
    ##  tcrossprod_sparse(m, only_upper = T, diag = F) 1.882685 1.899565 1.920442
    ##                                 tcrossprod((m)) 2.265994 2.288922 2.391366
    ##    median       uq      max neval cld
    ##  1.908928 1.909861 2.001169     5  a 
    ##  2.345652 2.466630 2.589634     5   b

    ## add minimum value
    microbenchmark(tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                                             expr      min
    ##  tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F) 1.642969
    ##                                                  tcrossprod((m)) 2.260224
    ##        lq     mean   median       uq      max neval cld
    ##  1.652664 1.668376 1.657611 1.686069 1.702566     5  a 
    ##  2.294994 2.326541 2.323772 2.361309 2.392406     5   b

Now, using a 'large' matrix (50.000 rows), that would normally crash
unless you have a good chunk of memory (worst case scenario: 50.000 \*
50.000 \* 8 bytes). Here we use the top\_n filter, for which the
required memory is at most nrow(m) \* top\_n \* 8 bytes.

    m = abs(Matrix::rsparsematrix(50000,20000,0.001))
    format(object.size(m), 'Mb')

    ## [1] "11.5 Mb"

    date = seq.Date(as.Date('1000-01-01'), as.Date('3010-01-10'), by=1)[1:nrow(m)]
    cp = tcrossprod_sparse(m, date=date, lwindow = 10, rwindow=10, verbose=T)

    cp = tcrossprod_sparse(m, top_n = 10, verbose=T)
    length(cp@x)

    ## [1] 500000
