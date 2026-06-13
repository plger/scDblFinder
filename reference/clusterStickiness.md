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
#> iter=0, 38 cells excluded from training.
#> iter=1, 40 cells excluded from training.
#> iter=2, 44 cells excluded from training.
#> Threshold found:0.734
#> 45 (3.9%) doublets called
clusterStickiness(sce)
#>     Estimate Std. Error    t value    p.value       FDR
#> 2  0.4852933  0.2323662  2.0884854 0.09107487 0.4553743
#> 3 -0.3071994  0.2804440 -1.0954038 0.32327665 1.0000000
#> 1  0.1633938  0.2659148  0.6144592 0.56579277 1.0000000
#> 4  0.1031945  0.2647532  0.3897761 0.71275204 1.0000000
#> 5 -0.1129954  0.3171997 -0.3562279 0.73620809 1.0000000
```
