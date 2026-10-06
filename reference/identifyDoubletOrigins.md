# identifyDoubletOrigins

Trains a classifier based on artificial doublets to identify the origins
of doublets.

## Usage

``` r
identifyDoubletOrigins(
  sce,
  clusters,
  samples = NULL,
  doublets = NULL,
  balance = TRUE,
  nArtificial = NULL,
  verbose = TRUE,
  xgb.param = list(booster = "gbtree", objective = "multi:softprob", eval_metric =
    "mlogloss", subsample = 0.5, colsample_bytree = 0.4, eta = 0.5, lambda = 100, alpha =
    1),
  max_rounds = 300,
  nthread = 1
)
```

## Arguments

- sce:

  A SingleCellExperiment object with a 'counts' assay.

- clusters:

  A vector of cluster labels for each column of \`sce\`, or the name of
  a colData column of \`sce\` containing such labels.

- samples:

  An optional vector of sample labels for each column of \`sce\`, or the
  name of a colData column of \`sce\` containing such labels. If
  provided, artificial doublets will be generated within-sample (but the
  classifier trained across samples).

- doublets:

  The doublets to be identified. If NULL, doublets will be taken from
  \`sce\$scDblFinder.class\` if present; if not, only training is
  performed.

- balance:

  Logical; whether to balance doublet types (default TRUE).

- nArtificial:

  Number of artificial doublets. If omitted, 100 per cluster combination
  will be used, up to a maximum of 20000.

- verbose:

  Logical; whether to output progress messages.

- xgb.param:

  A named list of parameters passed to xgboost.

- max_rounds:

  The maximum training round during cross-validation.

- nthread:

  The number of threads.

## Value

A list with: - model : the xgboost model - train_contigency : the
contingency matrix on the training data - predictions : the per-class
probabilities on \`doublets\` (if given) - calls : the origin calls on
\`doublets\` (if given) - features : the ordered features (i.e. genes)
needed to run the model.

## Examples

``` r
# we generate a random dataset
sce <- mockDoubletSCE(ncells = c(20,30,40), ngenes = 500)
# to have the example run fast, we set a low number of artificial doublets 
# and a low maximum learning rounds
res <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=10)
#> `doublets` not specified, and `sce` does not include doublet annotations. We will train themodel but not run predictions
#> Generating 100 artificial doublets.
#> Warning: scDblFinder might not work well with very low numbers of cells.
#> Running cross-validation...
#> Multiple eval metrics are present. Will use test_mlogloss for early stopping.
#> Will train until test_mlogloss hasn't improved in 2 rounds.
#> 
#> [1]  train-mlogloss:1.060435±0.002905    test-mlogloss:1.106107±0.027793 
#> [2]  train-mlogloss:1.026786±0.002313    test-mlogloss:1.079765±0.030362 
#> [3]  train-mlogloss:0.992117±0.002862    test-mlogloss:1.048461±0.028363 
#> [4]  train-mlogloss:0.960264±0.003601    test-mlogloss:1.023692±0.030042 
#> [5]  train-mlogloss:0.929879±0.005930    test-mlogloss:0.993621±0.027248 
#> [6]  train-mlogloss:0.902678±0.007541    test-mlogloss:0.968996±0.026746 
#> [7]  train-mlogloss:0.873568±0.006825    test-mlogloss:0.943397±0.026365 
#> [8]  train-mlogloss:0.847350±0.005659    test-mlogloss:0.918052±0.025004 
#> [9]  train-mlogloss:0.820982±0.006390    test-mlogloss:0.893774±0.025520 
#> [10] train-mlogloss:0.796431±0.006969    test-mlogloss:0.870935±0.024708 
#> Will use 9 rounds
#> Training final model...
#> Accuracy on artifical doublets:0.9899
# if desired, we could then re-run the same classifier on a new sample using
# predictDoubletOrigins()
```
