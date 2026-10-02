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
#> 1     71    84           A+B
#> 2     68    76           A+C
#> 3     25    22           B+C
#> 4     59    58           A+D
#> 5     25    30           B+D
#> 6      7    41           C+D
#> 7     83    21           A+E
#> 8     85    46           B+E
#> 9     40    74           C+E
#> 10    39    31           D+E
#> 11    35    54           A+F
#> 12    93    13           B+F
#> 13     9    63           C+F
#> 14    51    63           D+F
#> 15    95    47           E+F
```
