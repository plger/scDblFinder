# getCellPairs

Given a vector of cluster labels, returns pairs of cross-cluster cells

## Usage

``` r
getCellPairs(
  clusters,
  n = 1000,
  ls = NULL,
  q = c(0.1, 0.9),
  selMode = "proportional",
  soft.min = 5
)
```

## Arguments

- clusters:

  A vector of cluster labels for each cell, or a list containing
  metacells and graph

- n:

  The number of cell pairs to obtain

- ls:

  Optional library sizes

- q:

  Library size quantiles between which to include cells (ignored if
  \`ls\` is NULL)

- selMode:

  How to decide the number of pairs of each kind to produce. Either
  'proportional' (default, proportional to the abundance of the
  underlying clusters), 'uniform' or 'sqrt'.

- soft.min:

  Minimum number of pairs of a given type.

## Value

A data.frame with the columns

## Examples

``` r
# create random labels
x <- sample(head(LETTERS), 100, replace=TRUE)
getCellPairs(x, n=6)
#>    cell1 cell2 orig.clusters
#> 1     14    15           A+B
#> 2     62    39           A+C
#> 3     76    30           B+C
#> 4    100    11           A+D
#> 5     43    41           B+D
#> 6     90    64           C+D
#> 7     70    55           A+E
#> 8     25    83           B+E
#> 9     49    80           C+E
#> 10    19     6           D+E
#> 11    14    44           A+F
#> 12    98    45           B+F
#> 13    49    44           C+F
#> 14    56    58           D+F
#> 15    92    42           E+F
```
