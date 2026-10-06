# predictDoubletOrigins : run an origins classifier on new cells

predictDoubletOrigins : run an origins classifier on new cells

## Usage

``` r
predictDoubletOrigins(model, doublets, ret = c("call", "probs"))
```

## Arguments

- model:

  The output of
  [`identifyDoubletOrigins`](https://plger.github.io/scDblFinder/reference/identifyDoubletOrigins.md).

- doublets:

  A \`SingleCellExperiment\` or counts matrix of doublets.

- ret:

  Either 'call' (origin call, default) or 'probs' (per-class
  probabilities).

## Value

Either a factor (for \`ret="call"\`) or a matrix of probabilities.

## Examples

``` r
# we generate a random dataset
sce <- mockDoubletSCE(ncells = c(20,30,40), ngenes = 500)
# to have the example run fast, we set a low number of artificial doublets 
# and a low maximum learning rounds
clf <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=10)
#> `doublets` not specified, and `sce` does not include doublet annotations. We will train themodel but not run predictions
#> Generating 100 artificial doublets.
#> Warning: scDblFinder might not work well with very low numbers of cells.
#> Running cross-validation...
#> Multiple eval metrics are present. Will use test_mlogloss for early stopping.
#> Will train until test_mlogloss hasn't improved in 2 rounds.
#> 
#> [1]  train-mlogloss:1.059803±0.004535    test-mlogloss:1.091313±0.015154 
#> [2]  train-mlogloss:1.024378±0.002622    test-mlogloss:1.059191±0.018629 
#> [3]  train-mlogloss:0.987741±0.003865    test-mlogloss:1.024226±0.017201 
#> [4]  train-mlogloss:0.955461±0.002877    test-mlogloss:0.994742±0.017810 
#> [5]  train-mlogloss:0.925710±0.004317    test-mlogloss:0.968028±0.014986 
#> [6]  train-mlogloss:0.895066±0.003824    test-mlogloss:0.938593±0.017145 
#> [7]  train-mlogloss:0.865324±0.004363    test-mlogloss:0.912177±0.019816 
#> [8]  train-mlogloss:0.840360±0.006295    test-mlogloss:0.889653±0.019729 
#> [9]  train-mlogloss:0.816509±0.005654    test-mlogloss:0.869626±0.019209 
#> [10] train-mlogloss:0.793831±0.005581    test-mlogloss:0.849244±0.019835 
#> Will use 10 rounds
#> Training final model...
#> Accuracy on artifical doublets:0.99
# we run on a new set of doublets:
isDoublet <- which(sce$type=="doublet")
res <- predictDoubletOrigins(clf, sce[,isDoublet])
table(true=sce$origin[isDoublet], res)
#>                    res
#> true                cluster1+cluster2 cluster1+cluster3 cluster2+cluster3
#>   cluster1+cluster3                 1                 2                 0
#>   cluster2+cluster3                 2                 0                 1
```
