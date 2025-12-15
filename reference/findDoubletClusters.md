# Detect doublet clusters

Identify potential clusters of doublet cells based on whether they have
intermediate expression profiles, i.e., their profiles lie between two
other “source” clusters.

## Usage

``` r
findDoubletClusters(x, ...)

# S4 method for class 'ANY'
findDoubletClusters(
  x,
  clusters,
  subset.row = NULL,
  threshold = 0.05,
  get.all.pairs = FALSE,
  ...
)

# S4 method for class 'SummarizedExperiment'
findDoubletClusters(x, ..., assay.type = "counts")

# S4 method for class 'SingleCellExperiment'
findDoubletClusters(x, clusters = colLabels(x, onAbsence = "error"), ...)
```

## Arguments

- x:

  A numeric matrix-like object of count values, where each column
  corresponds to a cell and each row corresponds to an endogenous gene.

  Alternatively, a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  or
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing such a matrix.

- ...:

  For the generic, additional arguments to pass to specific methods.

  For the ANY method, additional arguments to pass to `findMarkers`.

  For the SummarizedExperiment method, additional arguments to pass to
  the ANY method.

  For the SingleCellExperiment method, additional arguments to pass to
  the SummarizedExperiment method.

- clusters:

  A vector of length equal to `ncol(x)`, containing cluster identities
  for all cells. If `x` is a SingleCellExperiment, this is taken from
  [`colLabels`](https://rdrr.io/pkg/SingleCellExperiment/man/colLabels.html)`(x)`
  by default.

- subset.row:

  See `?"scran-gene-selection"`.

- threshold:

  A numeric scalar specifying the FDR threshold with which to identify
  significant genes.

- get.all.pairs:

  Logical scalar indicating whether statistics for all possible source
  pairings should be returned.

- assay.type:

  A string specifying which assay values to use, e.g., `"counts"` or
  `"logcounts"`.

## Value

A DataFrame containing one row per query cluster with the following
fields:

- `source1`::

  String specifying the identity of the first source cluster.

- `source2`::

  String specifying the identity of the second source cluster.

- `num.de`::

  Integer, number of genes that are significantly non-intermediate in
  the query cluster compared to the two putative source clusters.

- `median.de`::

  Integer, median number of genes that are significantly
  non-intermediate in the query cluster across all possible source
  cluster pairings.

- `best`::

  String specifying the identify of the top gene with the lowest p-value
  against the doublet hypothesis for this combination of query and
  source clusters.

- `p.value`::

  Numeric, containing the adjusted p-value for the `best` gene.

- `lib.size1`::

  Numeric, ratio of the median library sizes for the first source
  cluster to the query cluster.

- `lib.size2`::

  Numeric, ratio of the median library sizes for the second source
  cluster to the query cluster.

- `prop`::

  Numeric, proportion of cells in the query cluster.

- `all.pairs`::

  A SimpleList object containing the above statistics for every pair of
  potential source clusters, if `get.all.pairs=TRUE`.

Each row is named according to its query cluster.

## Details

This function detects clusters of doublet cells in a manner similar to
the method used by Bach et al. (2017). For each “query” cluster, we
examine all possible pairs of “source” clusters, hypothesizing that the
query consists of doublets formed from the two sources. If so, gene
expression in the query cluster should be strictly intermediate between
the two sources after library size normalization.

We apply pairwise t-tests to the normalized log-expression profiles to
reject this null hypothesis. This is done by identifying genes that are
consistently up- or down-regulated in the query compared to *both*
sources. We count the number of genes that reject the null hypothesis at
the specified FDR `threshold`. For each query cluster, the most likely
pair of source clusters is that which minimizes the number of
significant genes.

Potential doublet clusters are identified using the following
characteristics, in order of importance:

- Low number of significant genes (i.e., `num.de`). Ideally, `median.de`
  is also high to indicate that the absence of strong DE is not due to a
  lack of power.

- A reasonable proportion of cells in the cluster, i.e., `prop`. This
  requires some expectation of the doublet rate in the experimental
  protocol.

- Library sizes of the source clusters that are below that of the query
  cluster, i.e., `lib.size*` values below unity. This assumes that the
  doublet cluster will contain more RNA and have more counts than either
  of the two source clusters.

For each query cluster, the function will only report the pair of source
clusters with the lowest `num.de`. Setting `get.all.pairs=TRUE` will
retrieve statistics for all pairs of potential source clusters. This can
be helpful for diagnostics to identify relationships between specific
clusters.

The reported `p.value` is of little use in a statistical sense, and is
only provided for inspection. Technically, it could be treated as the
Simes combined p-value against the doublet hypothesis for the query
cluster. However, this does not account for the multiple testing across
all pairs of clusters for each chosen cluster, especially as we are
chosing the pair that is most concordant with the doublet null
hypothesis.

We use library size normalization (via
[`librarySizeFactors`](https://rdrr.io/pkg/scuttle/man/librarySizeFactors.html))
even if existing size factors are present. This is because intermediate
expression of the doublet cluster is not guaranteed for arbitrary size
factors. For example, expression in the doublet cluster will be higher
than that in the source clusters if normalization was performed with
spike-in size factors.

## References

Bach K, Pensa S, Grzelak M, Hadfield J, Adams DJ, Marioni JC and Khaled
WT (2017). Differentiation dynamics of mammary epithelial cells revealed
by single-cell RNA sequencing. *Nat Commun.* 8, 1:2128.

## See also

`findMarkers`, to detect DE genes between clusters.

## Author

Aaron Lun

## Examples

``` r
# Mocking up an example.
library(SingleCellExperiment)
sce <- mockDoubletSCE(c(200,300,200))

# Compute doublet-ness of each cluster:
dbl <- findDoubletClusters(counts(sce), sce$cluster)
dbl
#> DataFrame with 3 rows and 9 columns
#>              source1     source2    num.de median.de        best      p.value
#>          <character> <character> <integer> <integer> <character>    <numeric>
#> cluster3    cluster2    cluster1       100       100      gene12 1.72469e-139
#> cluster1    cluster3    cluster2       105       105      gene29 8.90974e-102
#> cluster2    cluster3    cluster1       110       110     gene171 1.56657e-102
#>          lib.size1 lib.size2      prop
#>          <numeric> <numeric> <numeric>
#> cluster3  0.981618  1.014706  0.268817
#> cluster1  0.985507  0.967391  0.305108
#> cluster2  1.018727  1.033708  0.426075

# Narrow this down to clusters with very low 'N':
library(scuttle)
isOutlier(dbl$num.de, log=TRUE, type="lower")
#> [1] FALSE FALSE FALSE
#> attr(,"class")
#> [1] "outlier.filter" "logical"       
#> attr(,"thresholds")
#>   lower  higher 
#> 85.3746     Inf 

# Get help from "lib.size" below 1.
dbl$lib.size1 < 1 & dbl$lib.size2 < 1
#> [1] FALSE  TRUE FALSE
```
