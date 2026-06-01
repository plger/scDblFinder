# Recover intra-sample doublets

Recover intra-sample doublets that are neighbors to known inter-sample
doublets in a multiplexed experiment.

## Usage

``` r
recoverDoublets(x, ...)

# S4 method for class 'ANY'
recoverDoublets(
  x,
  doublets,
  samples,
  k = 50,
  transposed = FALSE,
  subset.row = NULL,
  BNPARAM = KmknnParam(),
  BPPARAM = SerialParam()
)

# S4 method for class 'SummarizedExperiment'
recoverDoublets(x, ..., assay.type = "logcounts")

# S4 method for class 'SingleCellExperiment'
recoverDoublets(x, ..., use.dimred = NULL)
```

## Arguments

- x:

  A log-expression matrix for all cells (including doublets) in columns
  and genes in rows. If `transposed=TRUE`, this should be a matrix of
  low-dimensional coordinates where each row corresponds to a cell.

  Alternatively, a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  or
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  containing (i) a log-expression matrix in the
  [`assays`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  as specified by `assay.type`, or (ii) a matrix of reduced dimensions
  in the
  [`reducedDims`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
  as specified by `use.dimred`.

- ...:

  For the generic, additional arguments to pass to specific methods.

  For the SummarizedExperiment method, additional arguments to pass to
  the ANY method.

  For the SingleCellExperiment method, additional arguments to pass to
  the SummarizedExperiment method.

- doublets:

  A logical, integer or character vector specifying which cells in `x`
  are known (inter-sample) doublets.

- samples:

  A numeric vector containing the relative proportions of cells from
  each sample, used to determine how many cells are to be considered as
  intra-sample doublets.

- k:

  Integer scalar specifying the number of nearest neighbors to use for
  computing the local doublet proportions.

- transposed:

  Logical scalar indicating whether `x` is transposed, i.e., cells in
  the rows.

- subset.row:

  A logical, integer or character vector specifying the genes to use for
  the neighbor search. Only used when `transposed=FALSE`.

- BNPARAM:

  A BiocNeighborParam object specifying the algorithm to use for the
  nearest neighbor search.

- BPPARAM:

  A
  [BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying the parallelization to use for the nearest neighbor
  search.

- assay.type:

  A string specifying which assay values contain the log-expression
  matrix.

- use.dimred:

  A string specifying whether existing values in
  [`reducedDims`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)`(x)`
  should be used.

## Value

A DataFrame containing one row per cell and the following fields:

- `proportion`, a numeric field containing the proportion of neighbors
  that are doublets.

- `known`, a logical field indicating whether this cell is a known
  inter-sample doublet.

- `predicted`, a logical field indicating whether this cell is a
  predicted intra-sample doublet.

The `metadata` contains `intra`, a numeric scalar containing the
expected number of intra-sample doublets.

## Details

In multiplexed single-cell experiments, we can detect doublets as
libraries with labels for multiple samples. However, this approach fails
to identify doublets consisting of two cells with the same label. Such
cells may be problematic if they are still sufficiently abundant to
drive formation of spurious clusters.

This function identifies intra-sample doublets based on the similarity
in expression profiles to known inter-sample doublets. For each cell, we
compute the proportion of the `k` neighbors that are known doublets. Of
the “unmarked” cells that are not known doublets, those with top \\X\\
largest proportions are considered to be intra-sample doublets. We use
`samples` to obtain a reasonable estimate for \\X\\, see the vignette
for details.

A larger value of `k` provides more stable estimates of the doublet
proportion in each cell. However, this comes at the cost of assuming
that each cell actually has `k` neighboring cells of the same state. For
example, if a doublet cluster has fewer than `k` members, its doublet
proportions will be “diluted” by inclusion of unmarked cells in the
next-closest cluster.

## See also

`doubletCells` and `doubletCluster`, for alternative methods of doublet
detection when no prior doublet information is available.

`hashedDrops` from the DropletUtils package, to identify doublets from
cell hashing experiments.

More detail on the mathematical background of this function is provided
in the corresponding vignette at
[`vignette("recoverDoublets", package="scDblFinder")`](https://plger.github.io/scDblFinder/articles/recoverDoublets.md).

## Author

Aaron Lun

## Examples

``` r
# Mocking up an example.
set.seed(100)
ngenes <- 1000
mu1 <- 2^rnorm(ngenes, sd=2)
mu2 <- 2^rnorm(ngenes, sd=2)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes) # Pure type 1
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes) # Pure type 2
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes) # Doublets (1 & 2)
all.counts <- cbind(counts.1, counts.2, counts.m)
lcounts <- scuttle::normalizeCounts(all.counts)
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")

# Pretending that half of the doublets are known. Also pretending that 
# the experiment involved two samples of equal size.
known <- 200 + seq_len(10) 
out <- recoverDoublets(lcounts, doublets=known, k=10, samples=c(1, 1))
out
#> DataFrame with 220 rows and 3 columns
#>     proportion     known predicted
#>      <numeric> <logical> <logical>
#> 1            0     FALSE     FALSE
#> 2            0     FALSE     FALSE
#> 3            0     FALSE     FALSE
#> 4            0     FALSE     FALSE
#> 5            0     FALSE     FALSE
#> ...        ...       ...       ...
#> 216        0.5     FALSE      TRUE
#> 217        0.5     FALSE      TRUE
#> 218        0.4     FALSE      TRUE
#> 219        0.3     FALSE      TRUE
#> 220        0.4     FALSE      TRUE
```
