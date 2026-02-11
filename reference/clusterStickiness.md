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
#> iter=2, 44 cells excluded from training.
#> Threshold found:0.688
#> 45 (3.8%) doublets called
clusterStickiness(sce)
#>     Estimate Std. Error    t value    p.value       FDR
#> 4  0.5669779  0.2048421  2.7678773 0.03946444 0.1973222
#> 5 -1.1509962  0.5242631 -2.1954552 0.07954797 0.3181919
#> 1  0.4464573  0.2077003  2.1495267 0.08429545 0.3181919
#> 3 -0.3525247  0.2725346 -1.2935046 0.25236718 0.5047344
#> 2  0.1105706  0.2378989  0.4647799 0.66162650 0.6616265
```
