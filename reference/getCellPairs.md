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
#> 1     88    36           A+B
#> 2     63    41           A+C
#> 3     97    20           B+C
#> 4     45    33           A+D
#> 5     40    80           B+D
#> 6      1    80           C+D
#> 7     63    56           A+E
#> 8     11     5           B+E
#> 9     53    14           C+E
#> 10    54     5           D+E
#> 11    24    61           A+F
#> 12    89    65           B+F
#> 13    20    61           C+F
#> 14    66    65           D+F
#> 15    57    16           E+F
```
