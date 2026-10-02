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
#> 1     22    28           A+B
#> 2     64    83           A+C
#> 3     25    79           B+C
#> 4     57    26           A+D
#> 5     12    17           B+D
#> 6     53     3           C+D
#> 7     46    92           A+E
#> 8     34    72           B+E
#> 9      6    54           C+E
#> 10    17    29           D+E
#> 11     7    86           A+F
#> 12     9    38           B+F
#> 13   100    27           C+F
#> 14     3    35           D+F
#> 15    93    99           E+F
```
