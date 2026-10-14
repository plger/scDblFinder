#' identifyDoubletOrigins
#' 
#' Trains a classifier based on artificial doublets to identify the origins of 
#' doublets.
#'
#' @param sce A SingleCellExperiment object with a 'counts' assay.
#' @param clusters A vector of cluster labels for each column of `sce`, or the 
#'   name of a colData column of `sce` containing such labels.
#' @param samples An optional vector of sample labels for each column of `sce`,
#'   or the name of a colData column of `sce` containing such labels. If 
#'   provided, artificial doublets will be generated within-sample (but the 
#'   classifier trained across samples).
#' @param doublets The doublets to be identified. If NULL, doublets will be 
#'   taken from `sce$scDblFinder.class` if present; if not, only training is
#'   performed.
#' @param balance Logical; whether to balance doublet types (default TRUE).
#' @param nArtificial Number of artificial doublets. If omitted, 100 per 
#'   cluster combination will be used, up to a maximum of 20000.
#' @param verbose Logical; whether to output progress messages.
#' @param xgb.param A named list of parameters passed to xgboost.
#' @param max_rounds The maximum training round during cross-validation.
#' @param nthread The number of threads.
#'
#' @returns A list with:
#'    - model : the xgboost model
#'    - train_contigency : the contingency matrix on the training data
#'    - predictions : the per-class probabilities on `doublets` (if given)
#'    - calls : the origin calls on `doublets` (if given)
#'    - stats : the observed and expected number of doublets of each type (and
#'      other related statistics), updated using the classifier's calls (if 
#'      `doublets` is given and `sce` contains `scDblFinder` statistics). When 
#'      multiple samples are present, this is aggregated across samples.
#'    - stats.per.sample : same as `stats`, but with a row per combination and
#'      sample (identified in the `sample` column). Only available when the 
#'      statistics of `sce` are sample-specific.
#'    - features : the ordered features (i.e. genes) needed to run the model.
#' @export
#' @importFrom xgboost xgb.train
#'
#' @examples
#' # we generate a random dataset
#' sce <- mockDoubletSCE(ncells = c(20,30,40), ngenes = 500)
#' # to have the example run fast, we set a low number of artificial doublets 
#' # and a low maximum learning rounds
#' res <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=10)
#' # if desired, we could then re-run the same classifier on a new sample using
#' # predictDoubletOrigins()
identifyDoubletOrigins <- function(sce, clusters, samples=NULL, doublets=NULL,
                                   balance=TRUE, nArtificial=NULL, verbose=TRUE,
                                   xgb.param=list(
                                     booster = "gbtree",
                                     objective = "multi:softprob",
                                     eval_metric = "mlogloss",
                                     subsample = 0.5,
                                     colsample_bytree = 0.4,
                                     eta = 0.5,
                                     lambda = 100,
                                     alpha = 1
                                   ), max_rounds=300, nthread=1){
  
  stopifnot(inherits(sce, "SingleCellExperiment"))
  if(is.null(xgb.param$nthread)) xgb.param$nthread <- nthread
  
  if(is.null(doublets)){
    if(is.null(sce$scDblFinder.class)){
      if(verbose) message("`doublets` not specified, and `sce` does not ",
                          "include doublet annotations. We will train the",
                          "model but not run predictions")
    }else{
      doublets <- which(sce$scDblFinder.class=="doublet")
    }
  }
  allSamples <- .checkColArg(sce, samples)
  dblSamples <- NULL
  if(!is.null(doublets)){
    if(inherits(doublets, "SingleCellExperiment")){
      doublets <- assay(doublets, "counts")
    }
    if(!is.array(doublets)){
      # `doublets` indexes the columns of `sce`, so we know their samples
      if(!is.null(allSamples)){
        idx <- if(is.character(doublets)) match(doublets, colnames(sce)) else
          seq_len(ncol(sce))[doublets]
        dblSamples <- as.character(allSamples[idx])
      }
      doublets <- assay(sce, "counts")[,doublets]
    }
  }
  if(!is.null(sce$scDblFinder.class)){
    w <- which(sce$scDblFinder.class!="doublet")
  }else{
    w <- seq_len(ncol(sce))
  }
  
  clusters <- droplevels(as.factor(.checkColArg(sce, clusters)[w]))
  samples <- .checkColArg(sce, samples)[w]
  if(is.null(samples)){
    nSamples <- 1L
  }else{
    samples <- droplevels(as.factor(samples))
    nSamples <- length(unique(samples))
  }

  if(is.null(nArtificial))
    nArtificial <- min((100/nSamples)*length(levels(clusters))^2,
                       20000/nSamples)
  
  if(verbose) message("Generating ", nArtificial*nSamples,
                      " artificial doublets.")
  
  out <- scDblFinder(sce[,w], clusters=clusters, samples=samples,
                     returnType="counts", artificialDoublets=nArtificial, 
                     nfeatures=1000, propMarkers=0.5, verbose=FALSE,
                     selMode=ifelse(balance, "uniform", "proportional"),
                     meta.triplets=FALSE, propRandom=0)
  if(!is.null(doublets)) doublets <- doublets[row.names(out),]
  
  w <- which(out$type=="doublet" & !is.na(out$origin))
  ad <- assay(out)[,w]
  label <- droplevels(as.factor(out$origin[w]))
  ad <- t(ad)/colSums(ad)
  xgb.param$num_class <- length(unique(label))
  xgb.param$max_depth <- pmax(3,pmin(round(sqrt(xgb.param$num_class)), 7))
  
  dtrain <- xgb.DMatrix(data = ad, label = as.integer(label) - 1L)
  
  if(verbose) message("Running cross-validation...")
  cv <- xgb.cv(
    params = xgb.param,
    data = dtrain,
    nrounds = max_rounds,
    nfold = 5,
    early_stopping_rounds = 2,
    eval_set = 0.2,
    verbose = verbose
  )
  best_nrounds <- cv$early_stop$best_iteration
  
  e <- cv$evaluation_log
  testm <- grep("test.+mean", colnames(e))
  best <- which.min(e[,testm])
  ac <- e[[testm]][best] + e[[grep("test.+std", colnames(e))]][best]
  best_nrounds <- min(which(e[[grep("test.+mean", colnames(e))]] <= ac))
  
  if(is.null(best_nrounds)){
    warning("Cross-validation did not reach plateau, using max rounds")
    best_nrounds <- max_rounds
  }else if(verbose){
    message("Will use ", best_nrounds, " rounds")
  }
  
  if(verbose) message("Training final model...")
  model <- xgb.train(
    params = xgb.param,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = verbose,
  )
  
  preds1 <- .xgpreds(predict(model, ad, type="class"), levels(label))
  tt <- unclass(table(label, factor(apply(preds1, 1, which.max),
                                    seq_len(ncol(preds1)),
                                    levels(label))))
  ac <- sum(diag(tt))/sum(tt)
  if(verbose) message("Accuracy on artifical doublets:", round(ac,4))
  if(ac<.9) warning("Low classifier accuracy on training data!")
  
  stats <- stats.per.sample <- calls <- preds2 <- NULL
  if(!is.null(doublets)){
    if(verbose) message("Predicting origins of real doublets")
    doublets <- t(doublets)/colSums(doublets)
    preds2 <- .xgpreds(predict(model, doublets), levels(label))
    calls <- factor(apply(preds2, 1, which.max),
                    seq_len(ncol(preds2)), colnames(preds2))
    row.names(preds2) <- names(calls) <- row.names(doublets)
    if(!is.null(metadata(sce)$scDblFinder.stats)){
      st <- .updateDoubletOriginsStats(sce, calls, dblSamples)
      stats <- st$stats
      stats.per.sample <- st$per.sample
    }
  }
  
  list(
    model=model,
    train_contigency=tt,
    predictions=preds2,
    calls=calls,
    stats=stats,
    stats.per.sample=stats.per.sample,
    features=colnames(ad))
  
}

.xgpreds <- function(pred, clnames){
  # no cells to predict (e.g. no doublets called): xgboost returns an empty object
  if(length(pred)==0L)
    return(matrix(numeric(0), nrow=0, ncol=length(clnames),
                  dimnames=list(NULL, clnames)))
  if(is.vector(pred)) {
    pred <- matrix(pred, ncol=length(clnames), byrow=TRUE)
  }
  stopifnot(ncol(pred) == length(clnames))
  colnames(pred) <- clnames
  pred
}

#' predictDoubletOrigins : run an origins classifier on new cells
#'
#' @param model The output of \code{\link{identifyDoubletOrigins}}.
#' @param doublets A `SingleCellExperiment` or counts matrix of doublets.
#' @param ret Either 'call' (origin call, default) or 'probs' (per-class 
#'   probabilities).
#'
#' @returns Either a factor (for `ret="call"`) or a matrix of probabilities.
#' @export
#'
#' @examples
#' # we generate a random dataset
#' sce <- mockDoubletSCE(ncells = c(20,30,40), ngenes = 500)
#' # to have the example run fast, we set a low number of artificial doublets 
#' # and a low maximum learning rounds
#' clf <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=10)
#' # we run on a new set of doublets:
#' isDoublet <- which(sce$type=="doublet")
#' res <- predictDoubletOrigins(clf, sce[,isDoublet])
#' table(true=sce$origin[isDoublet], res)
predictDoubletOrigins <- function(model, doublets, ret=c("call","probs")){
  ret <- match.arg(ret)
  stopifnot(is.list(model) && all(c("model", "features") %in% names(model)))
  if(is(doublets, "SingleCellExperiment")) doublets <- counts(doublets)
  doublets <- t(doublets[model$features,])
  res <- predict(model$model, newdata=doublets)
  colnames(res) <- colnames(model$train_contigency)
  if(ret=="probs") return(res)
  factor(apply(res, 1, which.max), seq_len(ncol(res)), colnames(res))
}


# updates the doublet statistics of `sce` using the multi-label classifier 
# results. Returns a list with `stats` (aggregated across samples)
# and `per.sample` (if relevant).
.updateDoubletOriginsStats <- function(s, calls, samples=NULL) {
  if(is(s, "SingleCellExperiment")) s <- metadata(s)$scDblFinder.stats
  s <- as.data.frame(s)
  
  .update_one <- function(df, calls){
    obs <- tabulate(calls, nlevels(calls))
    names(obs) <- levels(calls)
    df$observed <- as.numeric(obs[as.character(df$combination)])
    df$deviation <- abs(df$expected - df$observed)
    df$prop.deviation <- df$deviation / sum(df$expected)
    df
  }
  
  if(is.null(s$sample)){
    return(list(stats=.update_one(s, calls), per.sample=NULL))
  }
  
  cols <- setdiff(colnames(s), "sample")
  s$sample <- as.character(s$sample)
  ss <- split(s, factor(s$sample, unique(s$sample)))
  per.sample <- NULL
  if(!is.null(samples)){
    samples <- as.character(samples)
    if(!all(samples %in% names(ss))){
      warning("The sample labels do not match those of the existing doublet ",
              "statistics; per-sample statistics were not updated.")
    }else{
      per.sample <- do.call(rbind, lapply(names(ss), FUN=function(n){
        data.frame(sample=n, .update_one(ss[[n]][,colnames(s)!="sample"],
                                         calls[which(samples==n)]),
                   row.names=NULL)
      }))
      rownames(per.sample) <- NULL
    }
  }
  
  other <- setdiff(colnames(s), c("sample","combination","observed","expected",
                                  "deviation","prop.deviation"))
  s <- do.call(rbind, lapply(split(s, factor(s$combination, 
                                              unique(s$combination))), 
                             FUN=function(x){
    w <- if(sum(x$expected)>0) x$expected else rep(1,nrow(x))
    y <- data.frame(combination=x$combination[1], observed=NA_real_,
                    expected=sum(x$expected))
    for(f in other) y[[f]] <- weighted.mean(x[[f]], w, na.rm=TRUE)
    y
  }))
  s <- .update_one(s, calls)
  s <- s[,intersect(cols, colnames(s))]
  rownames(s) <- NULL
  list(stats=s, per.sample=per.sample)
}
