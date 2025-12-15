# plotDoubletMap

Plots a heatmap of observed versus expected doublets. Requires the
\`ComplexHeatmap\` package.

## Usage

``` r
plotDoubletMap(
  sce,
  colorBy = "enrichment",
  labelBy = "observed",
  addSizes = TRUE,
  col = NULL,
  column_title = "Clusters",
  row_title = "Clusters",
  column_title_side = "bottom",
  na_col = "white",
  ...
)
```

## Arguments

- sce:

  A SingleCellExperiment object on which \`scDblFinder\` has been run
  with the cluster-based approach.

- colorBy:

  Determines the color mapping. Either "enrichment" (for log2-enrichment
  over expectation) or any column of
  \`metadata(sce)\$scDblFinder.stats\`

- labelBy:

  Determines the cell labels. Either "enrichment" (for log2-enrichment
  over expectation) or any column of
  \`metadata(sce)\$scDblFinder.stats\`

- addSizes:

  Logical; whether to add the sizes of clusters to labels

- col:

  The colors scale to use (passed to \`ComplexHeatmap::Heatmap\`)

- column_title:

  passed to \`ComplexHeatmap::Heatmap\`

- row_title:

  passed to \`ComplexHeatmap::Heatmap\`

- column_title_side:

  passed to \`ComplexHeatmap::Heatmap\`

- na_col:

  color for NA cells

- ...:

  passed to \`ComplexHeatmap::Heatmap\`

## Value

a Heatmap object
