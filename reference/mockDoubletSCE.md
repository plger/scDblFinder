# mockDoubletSCE

Creates a mock random single-cell experiment object with doublets

## Usage

``` r
mockDoubletSCE(
  ncells = c(200, 300),
  ngenes = 200,
  mus = NULL,
  dbl.rate = 0.1,
  only.heterotypic = TRUE
)
```

## Arguments

- ncells:

  A positive integer vector indicating the number of cells per cluster
  (min 2 clusters)

- ngenes:

  The number of genes to simulate. Ignored if \`mus\` is given.

- mus:

  A list of cluster averages.

- dbl.rate:

  The doublet rate

- only.heterotypic:

  Whether to create only heterotypic doublets

## Value

A SingleCellExperiment object, with the colData columns \`type\`
indicating whether the cell is a singlet or doublet, and \`cluster\`
indicating from which cluster (or cluster combination) it was simulated.

## Examples

``` r
sce <- mockDoubletSCE()
```
