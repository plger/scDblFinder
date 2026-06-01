# selFeatures

Selects features based on cluster-wise expression or marker detection,
or a combination.

## Usage

``` r
selFeatures(
  sce,
  clusters = NULL,
  nfeatures = 1000,
  propMarkers = 0,
  FDR.max = 0.05
)
```

## Arguments

- sce:

  A
  [`SummarizedExperiment-class`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  [`SingleCellExperiment-class`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  with a 'counts' assay.

- clusters:

  Optional cluster assignments. Should either be a vector of labels for
  each cell.

- nfeatures:

  The number of features to select.

- propMarkers:

  The proportion of features to select from markers (rather than on the
  basis of high expression). Ignored if \`clusters\` isn't given.

- FDR.max:

  The maximum marker binom FDR to be included in the selection. (see
  [`findMarkers`](https://rdrr.io/pkg/scran/man/findMarkers.html)).

## Value

A vector of feature (i.e. row) names.

## Examples

``` r
sce <- mockDoubletSCE()
selFeatures(sce, clusters=sce$cluster, nfeatures=5)
#> Warning: 'sumCountsAcrossCells' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' instead.
#> See help("Deprecated")
#> Warning: 'summarizeAssayByGroup' is deprecated.
#> Use 'scrapper::aggregateAcrossCells' or 'beachmat::tatami.sums.by.group' instead.
#> [1] "gene13"  "gene149" "gene70"  "gene84"  "gene171"
```
