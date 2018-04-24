Example
============

Several options for filtering results. Filtering is performed on the fly, after each iteration (rows in m).

```{r}
set.seed(1)
m = Matrix::rsparsematrix(5,10,0.5)
tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = T)
tcrossprod_sparse(m, min_value = 0, only_upper = F, diag = F)
tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F)
tcrossprod_sparse(m, min_value = 0.2, only_upper = T, diag = F)
tcrossprod_sparse(m, min_value = 0, only_upper = T, diag = F, top_n = 1)
```

Speed seems to be comparable to Matrix::crossprod(), with some speed-ups if filtering is used.
This needs to be more properly tested while taking memory into account.

```{r}
m = abs(Matrix::rsparsematrix(3000,10000,0.1))

## quick and dirty L2 norm (for cosine similarity)
m = as(m, 'dgTMatrix')
norm = sqrt(Matrix::colSums(m^2))
m@x = m@x / norm[m@j+1]
m = as(m, 'dgCMatrix')

library(Matrix)
cp1 = tcrossprod_sparse(m)
cp2 = tcrossprod(m)
identical(as(cp1, 'dgCMatrix'), as(cp2, 'dgCMatrix'))

## regular
microbenchmark(tcrossprod_sparse(m),
               tcrossprod((m)),
               times=5)

## only calculate upper triangle and do not calculate diagonal
microbenchmark(tcrossprod_sparse(m, only_upper = T, diag = F),
               tcrossprod((m)),
               times=5)

## add minimum value
microbenchmark(tcrossprod_sparse(m, min_value = 0.5, only_upper = T, diag = F),
               tcrossprod((m)),
               times=5)
```

Now, using a 'large' matrix (50.000 rows), that would normally crash unless you have a good chunk of memory (worst case scenario: 500.000 * 50.000 * 8 bytes).

It will still take a while though, since you still have (at worst) 50.000 * 50.000 combinations. So this is where the verbose argument is usefull.

```{r}
m = abs(Matrix::rsparsematrix(50000,10000,0.1))
format(object.size(m), 'Mb')

cp = tcrossprod_sparse(m, only_upper = T, diag = F, verbose=T, top_n = 10, verbose=T)
length(cp@x)
```
