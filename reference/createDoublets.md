# createDoublets

Creates artificial doublet cells by combining given pairs of cells

## Usage

``` r
createDoublets(
  x,
  dbl.idx,
  clusters = NULL,
  resamp = 0.5,
  halfSize = 0.5,
  adjustSize = FALSE,
  prefix = "dbl."
)
```

## Arguments

- x:

  A count matrix of real cells

- dbl.idx:

  A matrix or data.frame with pairs of cell indexes stored in the first
  two columns.

- clusters:

  An optional vector of cluster labels (for each column of \`x\`)

- resamp:

  Logical; whether to resample the doublets using the poisson
  distribution. Alternatively, if a proportion between 0 and 1, the
  proportion of doublets to resample.

- halfSize:

  Logical; whether to half the library size of doublets (instead of just
  summing up the cells). Alternatively, a number between 0 and 1 can be
  given, determining the proportion of the doublets for which to perform
  the size adjustment. Ignored if not resampling.

- adjustSize:

  Logical; whether to adjust the size of the doublets using the median
  sizes per cluster of the originating cells. Requires \`clusters\` to
  be given. Alternatively to a logical value, a number between 0 and 1
  can be given, determining the proportion of the doublets for which to
  perform the size adjustment.

- prefix:

  Prefix for the colnames generated.

## Value

A matrix of artificial doublets.

## Examples

``` r
sce <- mockDoubletSCE()
idx <- getCellPairs(sce$cluster, n=200)
art.dbls <- createDoublets(sce, idx)
```
