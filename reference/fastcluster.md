# fastcluster

Performs a fast two-step clustering: first clusters using k-means with a
very large k, then uses louvain clustering of the k cluster averages and
reports back the cluster labels.

## Usage

``` r
fastcluster(
  x,
  k = NULL,
  rdname = "PCA",
  nstart = 3,
  iter.max = 50,
  ndims = NULL,
  nfeatures = 1000,
  verbose = TRUE,
  returnType = c("clusters", "preclusters", "metacells", "graph"),
  ...
)
```

## Arguments

- x:

  An object of class SCE

- k:

  The number of k-means clusters to use in the primary step (should be
  much higher than the number of expected clusters). Defaults to 1/10th
  of the number of cells with a maximum of 3000.

- rdname:

  The name of the dimensionality reduction to use.

- nstart:

  Number of starts for k-means clustering

- iter.max:

  Number of iterations for k-means clustering

- ndims:

  Number of dimensions to use

- nfeatures:

  Number of features to use (ignored if \`rdname\` is given and the
  corresponding dimensional reduction exists in \`sce\`)

- verbose:

  Logical; whether to output progress messages

- returnType:

  See return.

- ...:

  Arguments passed to \`scater::runPCA\` (e.g. BPPARAM or BSPARAM) if
  \`x\` does not have \`rdname\`.

## Value

By default, a vector of cluster labels. If \`returnType='preclusters'\`,
returns the k-means pre-clusters. If \`returnType='metacells'\`, returns
the metacells aggretated by pre-clusters and the corresponding cell
indexes. If \`returnType='graph'\`, returns the graph of (meta-)cells
and the corresponding cell indexes.

## Examples

``` r
sce <- mockDoubletSCE()
sce$cluster <- fastcluster(sce)
#> Reduced dimension not found - running PCA...
#> Building KNN graph and clustering
```
