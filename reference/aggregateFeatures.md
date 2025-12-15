# aggregateFeatures

Aggregates similar features (rows).

## Usage

``` r
aggregateFeatures(
  x,
  dims.use = seq(2L, 12L),
  k = 1000,
  num_init = 3,
  use.mbk = NULL,
  use.subset = 20000,
  minCount = 1L,
  norm.fn = TFIDF,
  twoPass = FALSE,
  ...
)
```

## Arguments

- x:

  A integer/numeric (sparse) matrix, or a \`SingleCellExperiment\`
  including a \`counts\` assay.

- dims.use:

  The PCA dimensions to use for clustering rows.

- k:

  The approximate number of meta-features desired

- num_init:

  The number of initializations used for k-means clustering.

- use.mbk:

  Logical; whether to use minibatch k-means (see
  [`mbkmeans`](https://rdrr.io/pkg/mbkmeans/man/mbkmeans.html)). If
  NULL, the minibatch approach will be used if there are more than 30000
  features.

- use.subset:

  How many cells (columns) to use to cluster the features.

- minCount:

  The minimum number of counts for a region to be included.

- norm.fn:

  The normalization function to use on the un-clustered data (a function
  taking a count matrix as a single argument and returning a matrix of
  the same dimensions).
  [TFIDF](https://plger.github.io/scDblFinder/reference/TFIDF.md) by
  default.

- twoPass:

  Logical; whether to perform the procedure twice, so in the second
  round cells are aggregated based on the meta-features of the first
  round, before re-clustering the features. Ignored if the dataset has
  fewer than \`use.subset\` cells.

- ...:

  Passed to
  [`mbkmeans`](https://rdrr.io/pkg/mbkmeans/man/mbkmeans.html). Can for
  instance be used to pass the \`BPPARAM\` argument for multithreading.

## Value

An aggregated version of \`x\` (either an array or a
\`SingleCellExperiment\`, depending on the input). If \`x\` is a
\`SingleCellExperiment\`, the feature clusters will also be stored in
\`metadata(x)\$featureGroups\`
