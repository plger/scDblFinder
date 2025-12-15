# plotThresholds

Plots scores used for thresholding.

## Usage

``` r
plotThresholds(d, ths = (0:100)/100, dbr = NULL, dbr.sd = NULL, do.plot = TRUE)
```

## Arguments

- d:

  A data.frame of cell properties, with each row representing a cell, as
  produced by \`scDblFinder(..., returnType="table")\`.

- ths:

  A vector of thresholds between 0 and 1 at which to plot values.

- dbr:

  The expected (mean) doublet rate.

- dbr.sd:

  The standard deviation of the doublet rate, representing the
  uncertainty in the estimate.

- do.plot:

  Logical; whether to plot the data (otherwise will return the
  underlying data.frame).

## Value

A ggplot, or a data.frame if \`do.plot==FALSE\`.
