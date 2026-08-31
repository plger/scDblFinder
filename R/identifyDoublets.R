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
#' @param nArtificial Number of artificial doublets. If omitted, 100 per cluster
#'   combination will be used, up to a maximum of 20000.
#' @param verbose Logical; whether to output progress messages.
#' @param xgb.param A named list of parameters passed to xgboost.
#' @param max_rounds The maximum training round during cross-validation.
#' @param nthread The number of threads.
#'
#' @returns A list with:
#'    - model : the xgboost model,
#'    - train_contigency : the contingency matrix on the training data
#'    - predictions : the per-class probabilities on `doublets`
#'    - features : the ordered features (i.e. genes) needed to run the model.
#' @export
#'
#' @examples
#' # we generate a random dataset
#' sce <- mockDoubletSCE(ncells = c(20,30,40), ngenes = 500)
#' # to have the example run fast, we set a low number of artificial doublets 
#' # and a low maximum learning rounds
#' res <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=10)
identifyDoubletOrigins <- function(sce, clusters, samples=NULL, doublets=NULL,
                                   balance=TRUE, nArtificial=NULL, verbose=TRUE,
                                   xgb.param=list(
                                     booster = "gbtree",
                                     objective = "multi:softprob",
                                     eval_metric = "mlogloss",
                                     subsample = 0.8,
                                     colsample_bytree = 0.7,
                                     eta = 0.1,
                                     max_depth = 8
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
  if(!is.null(doublets)){
    if(inherits(doublets, "SingleCellExperiment")){
      doublets <- assay(doublets, "counts")
    }
    if(!is.array(doublets)){
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
    nArtificial <- min((100/nSamples)*length(levels(clusters))^2, 20000/nSamples)
  
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
  
  dtrain <- xgb.DMatrix(data = ad, label = as.integer(label) - 1L)
  
  if(verbose) message("Running cross-validation...")
  cv <- xgb.cv(
    params = xgb.param,
    data = dtrain,
    nrounds = max_rounds,
    nfold = 5,
    early_stopping_rounds = 2,
    verbose = verbose
  )
  best_nrounds <- cv$early_stop$best_iteration
  if(is.null(best_nrounds)){
    warning("Cross-validation did not reach plateau, using max rounds")
    best_nrounds <- max_rounds
  }
  
  if(verbose) message("Training final model...")
  model <- xgb.train(
    params = xgb.param,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = verbose,
  )
  
  preds1 <- predict(model, ad, type="class")
  tt <- unclass(table(label, apply(preds1, 1, which.max)))
  colnames(tt) <- levels(out$origin)
  ac <- sum(diag(tt))/sum(tt)
  if(verbose) message("Accuracy on artifical doublets:", round(ac,4))
  if(ac<.9) warning("Low classifier accuracy on training data!")
  
  colnames(tt) <- levels(out$origin)
  
  preds2 <- NULL
  if(!is.null(doublets)){
    if(verbose) message("Predicting origins of real doublets")
    doublets <- t(doublets)/colSums(doublets)
    preds2 <- predict(model, doublets, type="class")
    colnames(preds2) <- levels(out$origin)
  }
  
  list(
    model=model,
    train_contigency=tt,
    predictions=preds2,
    features=colnames(ad))
  
}
