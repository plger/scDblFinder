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
#> 1     32    46           A+B
#> 2     43    85           A+C
#> 3     30    13           B+C
#> 4     47    72           A+D
#> 5     99    61           B+D
#> 6      7    42           C+D
#> 7     54    24           A+E
#> 8     46    91           B+E
#> 9     13     5           C+E
#> 10    42    69           D+E
#> 11    38    79           A+F
#> 12    55    15           B+F
#> 13    96    36           C+F
#> 14    34    36           D+F
#> 15    74    78           E+F
```
