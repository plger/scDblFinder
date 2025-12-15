# cxds2

Calculates a coexpression-based doublet score using the method developed
by [Bais and Kostka
2020](https://doi.org/10.1093/bioinformatics/btz698). This is the
original implementation from the
[scds](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
package, but enabling scores to be calculated for all cells while the
gene coexpression is based only on a subset (i.e. excluding
known/artificial doublets) and making it robust to low sparsity.

## Usage

``` r
cxds2(x, whichDbls = c(), ntop = 500, binThresh = NULL)
```

## Arguments

- x:

  A matrix of counts, or a \`SingleCellExperiment\` containing a
  'counts'

- whichDbls:

  The columns of \`x\` which are known doublets.

- ntop:

  The number of top features to keep.

- binThresh:

  The count threshold to be considered expressed.

## Value

A cxds score or, if \`x\` is a \`SingleCellExperiment\`, \`x\` with an
added \`cxds_score\` colData column.

## References

<https://doi.org/10.1093/bioinformatics/btz698>

## Examples

``` r
sce <- mockDoubletSCE()
sce <- cxds2(sce)
# which is equivalent to
# sce$cxds_score <- cxds2(counts(sce))
```
