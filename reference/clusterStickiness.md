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
#> iter=1, 37 cells excluded from training.
#> iter=2, 37 cells excluded from training.
#> Threshold found:0.699
#> 41 (3.5%) doublets called
clusterStickiness(sce)
#>       Estimate Std. Error     t value    p.value       FDR
#> 1  0.570949414  0.2367866  2.41124002 0.06077094 0.3038547
#> 4  0.537158468  0.2379533  2.25741087 0.07359007 0.3038547
#> 5 -1.061269307  0.5524059 -1.92117660 0.11275844 0.3382753
#> 3 -0.463794435  0.3297274 -1.40659974 0.21855358 0.4371072
#> 2  0.004126825  0.2783859  0.01482412 0.98874583 0.9887458
```
