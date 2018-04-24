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

Speed seems to be comparable to Matrix::crossprod(), with some speed-ups
if filtering is used. This needs to be more properly tested while taking
memory into account.

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
    ##                  expr     min       lq     mean   median       uq      max
    ##  tcrossprod_sparse(m) 1.79955 1.831394 1.846825 1.853715 1.859422 1.890043
    ##       tcrossprod((m)) 2.20692 2.237982 2.255434 2.259957 2.277980 2.294329
    ##  neval cld
    ##      5  a 
    ##      5   b

    ## only calculate upper triangle and do not calculate diagonal
    microbenchmark(tcrossprod_sparse(m, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                            expr      min       lq     mean
    ##  tcrossprod_sparse(m, only_upper = T, diag = F) 1.550285 1.575845 1.593503
    ##                                 tcrossprod((m)) 2.179167 2.227546 2.259518
    ##    median       uq      max neval cld
    ##  1.584861 1.616926 1.639599     5  a 
    ##  2.230339 2.317215 2.343324     5   b

    ## add minimum value
    microbenchmark(tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F),
                   tcrossprod((m)),
                   times=5)

    ## Unit: seconds
    ##                                                             expr      min
    ##  tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F) 1.418497
    ##                                                  tcrossprod((m)) 2.323160
    ##        lq     mean   median       uq      max neval cld
    ##  1.434765 1.435541 1.435815 1.436961 1.451665     5  a 
    ##  2.326524 2.335215 2.331723 2.345275 2.349396     5   b

Now, using a 'large' matrix (50.000 rows), that would normally crash
unless you have a good chunk of memory (worst case scenario: 50.000 \*
50.000 \* 8 bytes). Here we use the top\_n filter, for which the
required memory is at most nrow(m) \* top\_n.

It will still take a while though, since you still have (at worst)
50.000 \* 50.000 combinations. So this is where the verbose argument is
usefull.

    m = abs(Matrix::rsparsematrix(50000,10000,0.02))
    format(object.size(m), 'Mb')

    ## [1] "114.5 Mb"

    cp = tcrossprod_sparse(m, top_n = 10, verbose=T)
    length(cp@x)

    ## [1] 500000
