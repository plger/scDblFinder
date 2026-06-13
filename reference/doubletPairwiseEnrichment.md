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
#> Warning: Some cells in `sce` have an extremely low read counts; note that these could trigger errors and might best be filtered out
#> Clustering cells...
#> 3 clusters
#> Creating ~500 artificial doublets...
#> Dimensional reduction
#> Evaluating kNN...
#> Training model...
#> iter=0, 30 cells excluded from training.
#> iter=1, 0 cells excluded from training.
#> iter=2, 0 cells excluded from training.
#> Threshold found:0.995
#> 30 (5.7%) doublets called
doubletPairwiseEnrichment(sce)
#> theta=0.0329346479541671
#>   combination log2enrich      p.value          FDR
#> 1         1+2   1.475139 1.098949e-07 3.296848e-07
#> 3         2+3  -1.363649 1.000000e+00 1.000000e+00
#> 2         1+3  -1.037319 1.000000e+00 1.000000e+00
```
