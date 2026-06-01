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
#> 1     48    23           A+B
#> 2     20    58           A+C
#> 3     92    62           B+C
#> 4     48    39           A+D
#> 5     23    36           B+D
#> 6     57    36           C+D
#> 7     54    30           A+E
#> 8     91    37           B+E
#> 9     17    26           C+E
#> 10    46    22           D+E
#> 11    89    83           A+F
#> 12    95    96           B+F
#> 13    59    80           C+F
#> 14     5    83           D+F
#> 15    99    71           E+F
```
