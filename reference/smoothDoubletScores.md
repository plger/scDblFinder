# Smooth Doublet Scores

This function applies a distance-weighted kNN smoothing to single-cell
doublet probabilities, amplifying high-probability doublet clusters
while suppressing isolated noisy predictions.

## Usage

``` r
smoothDoubletScores(
  x,
  knn = NULL,
  coords = "PCA",
  k = 20,
  alpha = 0.5,
  gamma = 1,
  weights = TRUE,
  decayByDistance = c("linear", "exp", "none"),
  scoreColumn = "scDblFinder.score",
  outColumn = "smoothedDoubletScore",
  ...
)
```

## Arguments

- x:

  A `SingleCellExperiment` object (containing the 'scDblFinder.score'
  colData column), or a numeric vector of doublet probabilities.

- knn:

  Optional k nearest neighbors. A list containing `index` and `distance`
  matrices, typically the output of
  [`findKNN`](https://rdrr.io/pkg/BiocNeighbors/man/findKNN.html). If
  `NULL`, the kNN graph will be computed.

- coords:

  Either a character scalar indicating the name of the reduced dimension
  to use for kNN computation, or a matrix of such reduced dimensions,
  with cells as rows and dimensions as columns. Ignored if `knn` is
  given. Default "PCA".

- k:

  Integer. The number of nearest neighbors to compute if `knn` is
  `NULL`. Defaults to 30.

- alpha:

  Numeric \[0, 1\]. The blending factor between the cell's own score (0)
  and the neighborhood consensus (1). Defaults to 0.5.

- gamma:

  Numeric \>= 1.0. The non-linear amplification exponent. Values above 1
  exponentially amplify clusters of high-probability cells. Defaults to
  2.0.

- weights:

  Logical. If `TRUE` (default), uses adaptive Gaussian kernel weighting
  based on continuous distance. If `FALSE`, applies uniform weights
  across all neighbors.

- decayByDistance:

  Whether (and how) to reduce the importance of the neighborhood for
  isolated cells (based on the first neighbor distance), forcing them to
  rely primarily on their own raw score. Default 'linear'.

- scoreColumn:

  Character. The column name in `colData(x)` containing the raw doublet
  scores. Defaults to `"scDblFinder.score"`.

- outColumn:

  Character. The column name to store the smoothed scores if `x` is a
  `SingleCellExperiment`. Defaults to `"smoothedDoubletScore"`.

- ...:

  Passed to
  [`findKNN`](https://rdrr.io/pkg/BiocNeighbors/man/findKNN.html).

## Value

If `x` is a `SingleCellExperiment`, returns the object with an added
`colData` column containing the smoothed scores. If `x` is a numeric
vector, returns a numeric vector of smoothed scores.
