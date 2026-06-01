# amuletFromCounts

A reimplementation of the Amulet doublet detection method for
single-cell ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021),
based on tile/peak counts. Note that this is only a fast approximation
to the original Amulet method, and \*performs considerably worse\*; for
an equivalent implementation, see
[`amulet`](https://plger.github.io/scDblFinder/reference/amulet.md).

## Usage

``` r
amuletFromCounts(x, maxWidth = 500L, exclude = c("chrM", "M", "Mt"))
```

## Arguments

- x:

  A \`SingleCellExperiment\` object, or a matrix of counts with cells as
  columns. If the rows represent peaks, it is recommended to limite
  their width (see details).

- maxWidth:

  the maximum width for a feature to be included. This is ignored unless
  \`x\` is a \`SingleCellExperiment\` with \`rowRanges\`.

- exclude:

  an optional \`GRanges\` of regions to be excluded. This is ignored
  unless \`x\` is a \`SingleCellExperiment\` with \`rowRanges\`.

## Value

If \`x\` is a \`SingleCellExperiment\`, returns the object with an
additional \`amuletFromCounts.q\` colData column. Otherwise returns a
vector of the amulet doublet q-values for each cell.

## Details

The rationale for the amulet method is that a single diploid cell should
not have more than two reads covering a single genomic location, and the
method looks for cells enriched with sites covered by more than two
reads. If the method is applied on a peak-level count matrix, however,
larger peaks can however contain multiple reads even though no single
nucleotide is covered more than once. Therefore, in such case we
recommend to limit the width of the peaks used for this analysis,
ideally to maximum twice the upper bound of the fragment size. For
example, with a mean fragment size of 250bp and standard deviation of
125bp, peaks larger than 500bp are very likely to contain
non-overlapping fragments, and should therefore be excluded using the
\`maxWidth\` argument.

## See also

[`amulet`](https://plger.github.io/scDblFinder/reference/amulet.md)

## Examples

``` r
x <- mockDoubletSCE()
x <- amuletFromCounts(x)
table(call=x$amuletFromCounts.q<0.05, truth=x$type)
#>        truth
#> call    singlet doublet
#>   FALSE     500       0
#>   TRUE        0      22
```
