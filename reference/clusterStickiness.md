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
#> iter=1, 43 cells excluded from training.
#> iter=2, 41 cells excluded from training.
#> Threshold found:0.657
#> 44 (3.8%) doublets called
clusterStickiness(sce)
#>     Estimate Std. Error    t value   p.value FDR
#> 2  0.4475349  0.3374576  1.3261958 0.2421237   1
#> 1  0.3410064  0.3572148  0.9546257 0.3836075   1
#> 3 -0.3727388  0.4072857 -0.9151778 0.4020786   1
#> 4  0.1818591  0.3587658  0.5069021 0.6337812   1
#> 5 -0.1333212  0.4380387 -0.3043595 0.7731142   1
```
