# propHomotypic

Computes the proportion of pairs expected to be made of elements from
the same cluster.

## Usage

``` r
propHomotypic(clusters)
```

## Arguments

- clusters:

  A vector of cluster labels

## Value

A numeric value between 0 and 1.

## Examples

``` r
clusters <- sample(LETTERS[1:5], 100, replace=TRUE)
propHomotypic(clusters)
#> [1] 0.2056
```
