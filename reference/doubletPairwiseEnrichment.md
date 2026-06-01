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
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> 4 clusters
#> Creating ~500 artificial doublets...
#> Dimensional reduction
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Evaluating kNN...
#> Training model...
#> iter=0, 19 cells excluded from training.
#> iter=1, 19 cells excluded from training.
#> iter=2, 19 cells excluded from training.
#> Threshold found:0.743
#> 19 (3.6%) doublets called
doubletPairwiseEnrichment(sce)
#> theta=0.0940531381130913
#>   combination log2enrich      p.value          FDR
#> 3         2+3  1.6771530 7.920289e-06 4.752173e-05
#> 6         3+4  0.1424823 1.954982e-01 9.774910e-01
#> 1         1+2 -1.2070349 1.000000e+00 1.000000e+00
#> 4         1+4 -1.0426637 1.000000e+00 1.000000e+00
#> 2         1+3 -0.7659583 1.000000e+00 1.000000e+00
#> 5         2+4 -0.6495687 1.000000e+00 1.000000e+00
```
