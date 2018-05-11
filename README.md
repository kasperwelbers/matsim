Example
=======

    devtools::install_github('kasperwelbers/matsim')

    library(matsim)

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 3.3.2

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

    ## Warning in strptime(xx, f <- "%Y-%m-%d", tz = "GMT"): unknown timezone
    ## 'default/Europe/Amsterdam'

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
    ##  tcrossprod_sparse(m) 2.615479 2.666884 2.755037 2.753979 2.818006
    ##       tcrossprod((m)) 2.755363 2.816572 2.946422 2.883948 3.033160
    ##       max neval
    ##  2.920836     5
    ##  3.243070     5

    ## only calculate upper triangle and do not calculate diagonal
    microbenchmark(tcrossprod_sparse(m, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                            expr      min       lq     mean
    ##  tcrossprod_sparse(m, only_upper = T, diag = F) 1.857430 1.905484 1.919407
    ##                                 tcrossprod((m)) 2.711924 2.726251 2.763739
    ##    median       uq      max neval
    ##  1.929205 1.946529 1.958389     5
    ##  2.762128 2.775303 2.843090     5

    ## add minimum value
    microbenchmark(tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                                             expr      min
    ##  tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F) 1.594236
    ##                                                  tcrossprod((m)) 2.743505
    ##        lq     mean   median       uq      max neval
    ##  1.594396 1.652358 1.640681 1.687367 1.745109     5
    ##  2.784257 2.891470 2.797344 3.041113 3.091129     5

Now, using a 'large' matrix (50.000 rows), that would normally crash
unless you have a good chunk of memory (worst case scenario: 50.000 \*
50.000 \* 8 bytes). Here we use the top\_n filter, for which the
required memory is at most nrow(m) \* top\_n \* 8 bytes.

    m = abs(Matrix::rsparsematrix(50000,50000,0.001))
    format(object.size(m), 'Mb')

    ## [1] "28.8 Mb"

    cp = tcrossprod_sparse(m, top_n = 10)
    length(cp@x)

    ## [1] 500000
