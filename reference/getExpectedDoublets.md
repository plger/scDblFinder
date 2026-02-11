# getExpectedDoublets

getExpectedDoublets

## Usage

``` r
getExpectedDoublets(x, dbr = NULL, only.heterotypic = TRUE, dbr.per1k = 0.008)
```

## Arguments

- x:

  A vector of cluster labels for each cell

- dbr:

  The expected doublet rate.

- only.heterotypic:

  Logical; whether to return expectations only for heterotypic doublets

- dbr.per1k:

  The expected proportion of doublets per 1000 cells.

## Value

The expected number of doublets of each combination of clusters

## Examples

``` r
# random cluster labels
cl <- sample(head(LETTERS,4), size=2000, prob=c(.4,.2,.2,.2), replace=TRUE)
getExpectedDoublets(cl)
#>      A+B      A+C      B+C      A+D      B+D      C+D 
#> 5.808704 4.766784 2.611776 4.870976 2.668864 2.190144 
```
