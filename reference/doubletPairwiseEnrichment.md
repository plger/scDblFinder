# doubletPairwiseEnrichment

Calculates enrichment in any type of doublet (i.e. specific combination
of clusters) over random expectation. Note that when applied to an
multisample object, this functions assumes that the cluster labels match
across samples.

## Usage

``` r
doubletPairwiseEnrichment(
  x,
  lower.tail = FALSE,
  sampleWise = FALSE,
  type = c("poisson", "binomial", "nbinom", "chisq"),
  inclDiff = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  A table of double statistics, or a SingleCellExperiment on which
  scDblFinder was run using the cluster-based approach.

- lower.tail:

  Logical; defaults to FALSE to test enrichment (instead of depletion).

- sampleWise:

  Logical; whether to perform tests sample-wise in multi-sample
  datasets. If FALSE (default), will aggregate counts before testing.

- type:

  Type of test to use.

- inclDiff:

  Logical; whether to regress out any effect of the identification
  difficulty in calculating expected counts

- verbose:

  Logical; whether to output eventual warnings/notes

## Value

A table of significances for each combination.

## Examples

``` r
sce <- mockDoubletSCE()
sce <- scDblFinder(sce, clusters=TRUE, artificialDoublets=500)
#> Clustering cells...
#> 4 clusters
#> Creating ~500 artificial doublets...
#> Dimensional reduction
#> Evaluating kNN...
#> Training model...
#> iter=0, 19 cells excluded from training.
#> iter=1, 19 cells excluded from training.
#> iter=2, 22 cells excluded from training.
#> Threshold found:0.827
#> 23 (4.3%) doublets called
doubletPairwiseEnrichment(sce)
#> theta=0.0522644304279938
#>   combination log2enrich      p.value          FDR
#> 1         1+2  2.1712385 9.877137e-11 5.926282e-10
#> 3         2+3 -1.1429288 1.000000e+00 1.000000e+00
#> 2         1+3 -0.9147971 1.000000e+00 1.000000e+00
#> 4         1+4 -0.8701050 1.000000e+00 1.000000e+00
#> 5         2+4 -0.8328621 1.000000e+00 1.000000e+00
#> 6         3+4 -0.6037463 1.000000e+00 1.000000e+00
```
