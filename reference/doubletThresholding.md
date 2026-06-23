# doubletThresholding

Sets the doublet scores threshold; typically called by
[`scDblFinder`](https://plger.github.io/scDblFinder/reference/scDblFinder.md).

## Usage

``` r
doubletThresholding(
  d,
  dbr = NULL,
  dbr.sd = NULL,
  dbr.per1k = 0.008,
  stringency = 0.5,
  p = 0.1,
  method = c("auto", "optim", "dbr", "griffiths"),
  perSample = TRUE,
  returnType = c("threshold", "call")
)
```

## Arguments

- d:

  A data.frame of cell properties, with each row representing a cell, as
  produced by \`scDblFinder(..., returnType="table")\`, or minimally
  containing a \`score\` column.

- dbr:

  The expected (mean) doublet rate. If \`d\` contains a \`cluster\`
  column, the doublet rate will be adjusted for homotypic doublets.

- dbr.sd:

  The standard deviation of the doublet rate, representing the
  uncertainty in the estimate. Ignored if \`method!="optim"\`.

- dbr.per1k:

  The expected proportion of doublets per 1000 cells.

- stringency:

  A numeric value \>0 and \<1 which controls the relative weight of
  false positives (i.e. real cells) and false negatives (artificial
  doublets) in setting the threshold. A value of 0.5 gives equal weight
  to both; a higher value (e.g. 0.7) gives higher weight to the false
  positives, and a lower to artificial doublets. Ignored if
  \`method!="optim"\`.

- p:

  The p-value threshold determining the deviation in doublet score.

- method:

  The thresholding method to use, either 'auto' (default, automatic
  selection depending on the available fields), 'optim' (optimization of
  misclassification rate and deviation from expected doublet rate),
  'dbr' (strictly based on the expected doublet rate), or 'griffiths'
  (cluster-wise number of median absolute deviation in doublet score).

- perSample:

  Logical; whether to perform thresholding individually for each sample.

- returnType:

  The type of value to return, either doublet calls (\`call\`) or
  thresholds (\`threshold\`).

## Value

A vector of doublet calls if \`returnType=="call"\`, or a threshold (or
vector of thresholds) if \`returnType=="threshold"\`.

## Examples

``` r
sce <- mockDoubletSCE()
d <- scDblFinder(sce, verbose=FALSE, returnType="table")
th <- doubletThresholding(d, dbr=0.05)
th
#> [1] 0.6180837
```
