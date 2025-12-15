# TFIDF

The Term Frequency - Inverse Document Frequency (TF-IDF) normalization,
as implemented in Stuart & Butler et al. 2019.

## Usage

``` r
TFIDF(x, sf = 10000)
```

## Arguments

- x:

  The matrix of occurrences

- sf:

  Scaling factor

## Value

An array of same dimensions as \`x\`

## Examples

``` r
m <- matrix(rpois(500,1),nrow=50)
m <- TFIDF(m)
```
