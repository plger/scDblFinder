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
#> iter=1, 41 cells excluded from training.
#> iter=2, 34 cells excluded from training.
#> Threshold found:0.618
#> 45 (3.9%) doublets called
clusterStickiness(sce)
#>     Estimate Std. Error    t value   p.value       FDR
#> 5  0.6565350  0.3648642  1.7993952 0.1318602 0.6593011
#> 1  0.4033159  0.4195934  0.9612066 0.3805937 1.0000000
#> 2 -0.2903152  0.4781229 -0.6071980 0.5702351 1.0000000
#> 4 -0.1490259  0.3994685 -0.3730606 0.7243962 1.0000000
#> 3  0.1296982  0.3850223  0.3368588 0.7499025 1.0000000
```
