# addDoublets

Adds artificial doublets to an existing dataset

## Usage

``` r
addDoublets(
  x,
  clusters,
  dbr = (0.01 * ncol(x)/1000),
  only.heterotypic = TRUE,
  adjustSize = FALSE,
  prefix = "doublet.",
  ...
)
```

## Arguments

- x:

  A count matrix of singlets, or a
  [`SummarizedExperiment-class`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- clusters:

  A vector of cluster labels for each column of \`x\`

- dbr:

  The doublet rate

- only.heterotypic:

  Whether to add only heterotypic doublets.

- adjustSize:

  Whether to adjust the library sizes of the doublets.

- prefix:

  Prefix for the colnames generated.

- ...:

  Any further arguments to
  [`createDoublets`](https://plger.github.io/scDblFinder/reference/createDoublets.md).

## Value

A \`SingleCellExperiment\` with the colData columns \`cluster\` and
\`type\` (indicating whether the cell is a singlet or doublet).

## Examples

``` r
sce <- mockDoubletSCE(dbl.rate=0)
sce <- addDoublets(sce, clusters=sce$cluster)
```
