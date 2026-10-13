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
#> iter=0, 43 cells excluded from training.
#> iter=1, 46 cells excluded from training.
#> iter=2, 43 cells excluded from training.
#> Threshold found:0.658
#> 48 (4.1%) doublets called
clusterStickiness(sce)
#>      Estimate Std. Error    t value    p.value       FDR
#> 5  0.46478892  0.1847161  2.5162333 0.05341973 0.2670987
#> 2 -0.13701315  0.2532653 -0.5409867 0.61173819 1.0000000
#> 3 -0.10894496  0.2168956 -0.5022922 0.63679654 1.0000000
#> 1 -0.08931500  0.2379699 -0.3753205 0.72281691 1.0000000
#> 4  0.03613934  0.2128837  0.1697609 0.87185298 1.0000000
```
