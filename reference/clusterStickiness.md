# clusterStickiness

Tests for enrichment of doublets created from each cluster (i.e.
cluster's stickiness). Only applicable with \>=4 clusters. Note that
when applied to an multisample object, this functions assumes that the
cluster labels match across samples.

## Usage

``` r
clusterStickiness(
  x,
  type = c("quasibinomial", "nbinom", "binomial", "poisson"),
  inclDiff = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  A table of double statistics, or a SingleCellExperiment on which
  [scDblFinder](https://plger.github.io/scDblFinder/reference/scDblFinder.md)
  was run using the cluster-based approach.

- type:

  The type of test to use (quasibinomial recommended).

- inclDiff:

  Logical; whether to include the difficulty in the model. If NULL, will
  be used only if there is a significant trend with the enrichment.

- verbose:

  Logical; whether to print additional running information.

## Value

A table of test results for each cluster.

## Examples

``` r
sce <- mockDoubletSCE(rep(200,5), dbl.rate=0.2)
sce <- scDblFinder(sce, clusters=TRUE, artificialDoublets=500)
#> Warning: Some cells in `sce` have an extremely low read counts; note that these could trigger errors and might best be filtered out
#> Clustering cells...
#> 5 clusters
#> Creating ~500 artificial doublets...
#> Dimensional reduction
#> Evaluating kNN...
#> Training model...
#> iter=0, 40 cells excluded from training.
#> iter=1, 43 cells excluded from training.
#> iter=2, 52 cells excluded from training.
#> Threshold found:0.753
#> 44 (3.8%) doublets called
clusterStickiness(sce)
#>      Estimate Std. Error     t value    p.value       FDR
#> 5  0.66860293  0.2988549  2.23721585 0.07547797 0.3773898
#> 2 -0.43630836  0.4405344 -0.99040691 0.36745252 1.0000000
#> 1  0.31394198  0.3417711  0.91857374 0.40046109 1.0000000
#> 4 -0.20239823  0.3490364 -0.57987706 0.58714505 1.0000000
#> 3 -0.02531098  0.3371028 -0.07508386 0.94305952 1.0000000
```
