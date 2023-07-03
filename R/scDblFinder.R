#' scDblFinder
#'
#' Identification of heterotypic (or neotypic) doublets in single-cell RNAseq
#' using cluster-based generation of artificial doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment-class}}, or array of
#' counts.
#' @param artificialDoublets The approximate number of artificial doublets to
#' create. If \code{NULL}, will be the maximum of the number of cells or
#' \code{5*nbClusters^2} (with a minimum of 1500).
#' @param clusters The optional cluster assignments. This is used to make
#' doublets more efficiently. \code{clusters} should either be a vector of
#' labels for each cell, or the name of a colData column of \code{sce}.
#' Alternatively, if `clusters=TRUE`, fast clustering will be performed. If
#' `clusters` is a single integer, it will determine how many clusters to
#' create (using k-means clustering). If `clusters` is NULL or FALSE, purely
#' random artificial doublets will be generated.
#' @param samples A vector of the same length as cells (or the name of a column
#' of \code{colData(x)}), indicating to which sample each cell belongs. Here, a
#' sample is understood as being processed independently. If omitted, doublets
#' will be searched for with all cells together. If given, doublets will be
#' searched for independently for each sample, which is preferable if they
#' represent different captures. If your samples were multiplexed using cell
#' hashes, want you want to give here are the different batches/wells (i.e.
#' independent captures, since doublets cannot arise across them) rather
#' than biological samples.
#' @param multiSampleMode Either "split" (recommended if there is
#' heterogeneity across samples), "singleModel", "singleModelSplitThres", or
#' "asOne" (see details below).
#' @param knownDoublets An optional logical vector of known doublets (e.g.
#' through cell barcodes), or the name of a colData column of `sce` containing
#' that information. The way these are used depends on the `knownUse` argument.
#' @param knownUse The way to use known doublets, either 'discard' (they are
#' discarded for the purpose of training, but counted as positive for 
#' thresholding) or 'positive' (they are used as positive doublets for training 
#' - usually leads to a mild decrease in accuracy due to the fact that known 
#' doublets typically include a sizeable fraction of homotypic doublets). Note
#' that `scDblFinder` does *not* enforce that the knownDoublets be necessarily
#' called as doublets in the final classification, if they are not predicted as 
#' such.
#' @param nfeatures The number of top features to use. Alternatively, a 
#'   character vectors of feature names (e.g. highly-variable genes) to use.
#' @param dims The number of dimensions used.
#' @param dbr The expected doublet rate. By default this is assumed to be 1\%
#' per thousand cells captured (so 4\% among 4000 thousand cells), which is
#' appropriate for 10x datasets. Corrections for homeotypic doublets will be
#' performed on the given rate.
#' @param dbr.sd The uncertainty range in the doublet rate, interpreted as
#' a +/- around `dbr`. During thresholding, deviation from the expected doublet
#' rate will be calculated from these boundaries, and will be considered null
#' within these boundaries. If NULL, will be 40\% of `dbr`. Set to `dbr.sd=0` to
#'  disable the uncertainty around the doublet rate, or to `dbr.sd=1` to disable
#'  any expectation of the number of doublets (thus letting the thresholding be
#'  entirely driven by the misclassification of artificial doublets).
#' @param k Number of nearest neighbors (for KNN graph). If more than one value
#' is given, the doublet density will be calculated at each k (and other values
#' at the highest k), and all the information will be used by the classifier.
#' If omitted, a reasonable set of values is used.
#' @param clustCor Include Spearman correlations to cell type averages in the
#' predictors. If `clustCor` is a matrix of cell type marker expressions (with
#' features as rows and cell types as columns), the subset of these which are
#' present in the selected features will be correlated to each cell to produce
#' additional predictors (i.e. one per cell type). Alternatively, if `clustCor`
#' is a positive integer, this number of inter-cluster markers will be selected
#' and used for correlation (se `clustCor=Inf` to use all available genes).
#' @param removeUnidentifiable Logical; whether to remove artificial doublets of
#'  a combination that is generally found to be unidentifiable.
#' @param includePCs The index of principal components to include in the
#' predictors (e.g. `includePCs=1:2`), or the number of top components to use
#' (e.g. `includePCs=10`, equivalent to 1:10).
#' @param propRandom The proportion of the artificial doublets which
#' should be made of random cells (as opposed to inter-cluster combinations).
#' If clusters is FALSE or NULL, this is ignored (and set to 1).
#' @param propMarkers The proportion of features to select based on marker
#' identification.
#' @param trainingFeatures The features to use for training (defaults to an
#' optimal pre-selection based on benchmark datasets). To exclude features
#' (rather than list those to be included), prefix them with a "-".
#' @param unident.th The score threshold below which artificial doublets will be
#' considered unidentifiable.
#' @param processing Counts (real and artificial) processing before KNN. Either
#' 'default' (normal \code{scater}-based normalization and PCA), "rawPCA" (PCA
#' without normalization), "rawFeatures" (no normalization/dimensional
#' reduction), "normFeatures" (uses normalized features, without PCA) or a
#' custom function with (at least) arguments `e` (the matrix of counts) and
#' `dims` (the desired number of dimensions), returning a named matrix with
#' cells as rows and components as columns.
#' @param returnType Either "sce" (default), "table" (to return the table of
#' cell attributes including artificial doublets), or "full" (returns an SCE
#' object containing both the real and artificial cells).
#' @param score Score to use for final classification.
#' @param metric Error metric to optimize during training (e.g. 'merror',
#' 'logloss', 'auc', 'aucpr').
#' @param nrounds Maximum rounds of boosting. If NULL, will be determined
#' through cross-validation. If a number <=1, will used the best
#' cross-validation round minus `nrounds` times the standard deviation of the
#' classification error.
#' @param max_depth Maximum depths of each tree.
#' @param iter A positive integer indicating the number of scoring iterations
#' (ignored if `score` isn't based on classifiers). At each iteration, real
#' cells that would be called as doublets are excluding from the training, and
#' new scores are calculated. Recommended values are 1 or 2.
#' @param threshold Logical; whether to threshold scores into binary doublet
#' calls
#' @param aggregateFeatures Whether to perform feature aggregation (recommended
#'  for ATAC). Can also be a positive integer, in which case this will indicate
#'  the number of components to use for feature aggregation (if TRUE, `dims`
#'  will be used.)
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when
#' `samples!=NULL`); otherwise passed to eventual PCA and K/SNN calculations.
#' @param ... further arguments passed to \code{\link{getArtificialDoublets}}.
#'
#' @return The \code{sce} object with several additional colData columns, in
#' particular `scDblFinder.score` (the final score used) and `scDblFinder.class`
#' (whether the cell is called as 'doublet' or 'singlet'). See
#' \code{vignette("scDblFinder")} for more details; for alternative return
#' values, see the `returnType` argument.
#'
#' @details
#' This function generates artificial doublets from real cells, evaluates their
#' prevalence in the neighborhood of each cells, and uses this along with
#' additional cell-level features to classify doublets. The approach is
#' complementary to doublets identified via cell hashes and SNPs in multiplexed
#' samples: the latter can identify doublets formed by cells of the same type
#' from two samples, which are nearly undistinguishable from real cells
#' transcriptionally, but cannot identify doublets made by cells of the
#' same sample. See \code{vignette("scDblFinder")} for more details on the
#' method.
#'
#' The `clusters` and `propRandom` argument determines whether the artificial
#' doublets are generated between clusters or randomly.
#'
#' When multiple samples/captures are present, they should be specified using
#' the \code{samples} argument. In this case, we recommend the use of
#' \code{BPPARAM} to perform several of the steps in parallel. Artificial
#' doublets and kNN networks will be computed separately; then the behavior will
#' then depend on the `multiSampleMode` argument:
#'
#' \itemize{
#'   \item \emph{split}: the whole process is split by sample. This is the
#'   default and recommended mode, because it is the most robust (e.g. to
#'   heterogeneity between samples, also for instance in the number of cells),
#'   and in practice we have not seen major gains in sharing information across
#'   samples;
#'   \item \emph{singleModel}: the doublets are generated on a per-sample basis,
#'   but the classifier and thresholding will be trained globally;
#'   \item \emph{singleModelSplitThres}: the doublets are generated on a
#'   per-sample basis, the classifier is trained globally, but the final
#'   thresholding is per-sample;
#'   \item \emph{asOne}: the doublet rate (if not given) is calculated as the
#'   weighted average of sample-specific doublet rates, and all samples are
#'   otherwise run as if they were one sample. This can get computationally
#'   more intensive, and can lead to biases if there are batch effects.
#' }
#'
#' When inter-sample doublets are available, they can be provided to
#' `scDblFinder` through the \code{knownDoublets} argument to improve the
#' identification of further doublets. How exactly these are used depends on the
#' `knownUse` argument: with 'discard' (default), the known doublets are
#' excluded from the training step, but counted as positives. With 'positive',
#' they are included and treated as positive doublets for the training step.
#' Note that because known doublets can in practice include a lot of homotypic
#' doublets, this second approach can often lead to a slight decrease in the
#' accuracy of detecting heterotypic doublets.
#'
#' Finally, for some types of data, such as single-cell ATAC-seq, selecting a
#' number of top features is ineffective due to the high sparsity of the signal.
#' In such contexts, rather than _selecting_ features we recommend to use the
#' alternative approach of _aggregating_ similar features (with
#' `aggregateFeatures=TRUE`), which strongly improves accuracy. See the
#' vignette for more detail.
#'
#' @import SingleCellExperiment BiocParallel
#' @importFrom SummarizedExperiment colData<- assayNames
#' @importFrom scuttle normalizeCounts
#' @importFrom scater runPCA
#' @importFrom methods is
#' @importFrom DelayedArray as.matrix
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocSingular IrlbaParam
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- mockDoubletSCE()
#' sce <- scDblFinder(sce)
#' table(truth=sce$type, call=sce$scDblFinder.class)
#'
#' @export
#' @rdname scDblFinder
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom BiocParallel SerialParam bpnworkers
scDblFinder <- function(
  sce, clusters=NULL, samples=NULL, clustCor=NULL, artificialDoublets=NULL,
  knownDoublets=NULL, knownUse=c("discard","positive"), dbr=NULL, dbr.sd=NULL, 
  nfeatures=1352, dims=20, k=NULL, removeUnidentifiable=TRUE, includePCs=19, 
  propRandom=0, propMarkers=0, aggregateFeatures=FALSE,
  returnType=c("sce","table","full","counts"),
  score=c("xgb","weighted","ratio"), processing="default", metric="logloss",
  nrounds=0.25, max_depth=4, iter=3, trainingFeatures=NULL, unident.th=NULL, 
  multiSampleMode=c("split","singleModel","singleModelSplitThres","asOne"),
  threshold=TRUE, verbose=is.null(samples), BPPARAM=SerialParam(), ...){

  multiSampleMode <- match.arg(multiSampleMode)

  ## check arguments
  sce <- .checkSCE(sce)
  score <- match.arg(score)
  knownUse <- match.arg(knownUse)
  if(!is.null(clustCor)){
    if(is.null(dim(clustCor)) && (!is.numeric(clustCor) || clustCor<0))
      stop("`clustCor` should be either a matrix of marker expression per cell",
      " types, or a positive integer indicating the number of markers to use.")
  }
  returnType <- match.arg(returnType)
  if(!is.null(clusters) && (!is.logical(clusters))){
    if(length(clusters)>1 || !is.numeric(clusters))
      clusters <- .checkColArg(sce, clusters)
    if(is.factor(clusters)) clusters <- droplevels(clusters)
  }
  if(is.null(unident.th))
    unident.th <- ifelse(is.null(clusters) || isFALSE(clusters), 0.2, 0)
  knownDoublets <- .checkColArg(sce, knownDoublets)
  samples <- .checkColArg(sce, samples)
  if(!is.null(samples)) samples <- as.factor(samples)
  .checkPropArg(propMarkers)
  .checkPropArg(propRandom)
  .checkPropArg(dbr.sd)
  .checkPropArg(dbr, acceptNull=TRUE)
  processing <- .checkProcArg(processing)

  if(!bpisup(BPPARAM)){
    ## pre-start params for independent seeds between bplapply calls
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM))
  }

  if(length(nfeatures)>1){
    if(!all(nfeatures %in% row.names(sce)))
      stop("'nfeatures' has a length >1, which is interpreted as feature (i.e.",
           " row) names to use, but not all of the features specified are ",
           "found in the object. ")
    sel_features <- nfeatures
    nfeatures <- length(sel_features)
  }else{
    ## if clusters are given, it's more efficient to do feature selection before
    ## eventually splitting the dataset
    if(!is.null(clusters) && length(clusters)>1 && !aggregateFeatures){
      sel_features <- selFeatures(sce, clusters, nfeatures=nfeatures,
                                  propMarkers=propMarkers)
    }else{
      sel_features <- row.names(sce)
    }
  }

  if(!is.null(samples) && multiSampleMode=="asOne"){
    if(is.null(dbr)){
      tt <- as.numeric(table(samples))
      dbr <- weighted.mean(tt/100000, tt)
    }
    samples <- NULL
  }
  if(!is.null(samples)){
    ## splitting by samples
    if(!(isSplitMode <- multiSampleMode=="split")) includePCs <- c()
    if(returnType=="full")
      warning("`returnType='full'` ignored when splitting by samples")
    cs <- split(seq_along(samples), samples, drop=TRUE)
    names(nn) <- nn <- names(cs)
    ## run scDblFinder individually
    d <- bplapply(nn, BPPARAM=BPPARAM, FUN=function(n){
      x <- cs[[n]]
      if(!is.null(clusters) && length(clusters)>1) clusters <- clusters[x]
      if(!is.null(knownDoublets) && length(knownDoublets)>1){
        knownDoublets <- knownDoublets[x]
        if(!any(knownDoublets)) knownDoublets <- NULL
      }
      out <- tryCatch(
        scDblFinder(sce[sel_features,x], clusters=clusters, dims=dims, dbr=dbr,
                    dbr.sd=dbr.sd, clustCor=clustCor, unident.th=unident.th,
                    knownDoublets=knownDoublets, knownUse=knownUse,
                    artificialDoublets=artificialDoublets, k=k,
                    processing=processing, nfeatures=nfeatures,
                    propRandom=propRandom, includePCs=includePCs,
                    propMarkers=propMarkers, trainingFeatures=trainingFeatures,
                    returnType=ifelse(returnType=="counts","counts","table"),
                    threshold=isSplitMode, score=ifelse(isSplitMode,score,"weighted"),
                    removeUnidentifiable=removeUnidentifiable, verbose=FALSE,
                    aggregateFeatures=aggregateFeatures, ...),
               error=function(e){
                 stop("An error occured while processing sample '",n,"':\n", e)
               })
      if(!is.matrix(out)) out$sample <- n
      out
    })
    if(returnType=="counts") return(do.call(cbind, d))

    ## aggregate the property tables
    d <- .aggResultsTable(d)
    if(multiSampleMode!="split"){
      ## score and thresholding
      d <- .scDblscore(d, scoreType=score, threshold=threshold, dbr=dbr,
                       dbr.sd=dbr.sd, max_depth=max_depth, nrounds=nrounds,
                       iter=iter, BPPARAM=BPPARAM, verbose=verbose,
                       features=trainingFeatures, unident.th=unident.th,
                       metric=metric, filterUnidentifiable=removeUnidentifiable,
                       perSample=multiSampleMode=="singleModelSplitThres",
                       includeSamples=TRUE)
    }
    if(returnType=="table") return(d)
    return(.scDblAddCD(sce, d))
  }

  ## Handling a single sample

  if(ncol(sce)<100)
    warning("scDblFinder might not work well with very low numbers of cells.")
  if(verbose && ncol(sce)>25000 && multiSampleMode!="asOne")
    warning("You are trying to run scDblFinder on a very large number of ",
            "cells. If these are from different captures, please specify this",
            " using the `samples` argument.", immediate=TRUE)

  k <- .defaultKnnKs(k, ncol(sce))

  orig <- sce
  wDbl <- c()

  ## if known doublets are given, we need to treat them separately
  if(!is.null(knownDoublets) && length(wDbl <- which(knownDoublets))>0){
    sce$knownDoublet <- knownDoublets
    sce.dbl <- sce[,wDbl,drop=FALSE]
    sce <- sce[,-wDbl,drop=FALSE]
    if(!is.null(clusters) && length(clusters)>1){
      clusters.dbl <- clusters[wDbl]
      clusters <- clusters[-wDbl]
      if(is.factor(clusters)) clusters <- droplevels(clusters)
    }
  }

  if(aggregateFeatures){
    if(verbose) message("Aggregating features...")
    if(is.numeric(aggregateFeatures)){
      fdims <- aggregateFeatures
    }else{
      fdims <- dims
    }
    if(length(fdims)==1) fdims <- seq_len(dims)[-1]
    sce <- aggregateFeatures(sce, dims.use=fdims, k=nfeatures)
    sel_features <- row.names(sce)
  }

  ## clustering (if required)
  if(isFALSE(clusters)) clusters <- NULL
  if(!is.null(clusters)){
    if(!is.null(clusters) && length(clusters)==1 && !isFALSE(clusters)){
      if(verbose) message("Clustering cells...")
      if(isTRUE(clusters)) clusters <- NULL
      if(!is.null(clusters)){
        clusters <- fastcluster(sce, ndims=dims, k=clusters, nfeatures=nfeatures,
                                returnType="preclusters",
                                BPPARAM=BPPARAM, verbose=FALSE)
      }else{
        clusters <- fastcluster(sce, ndims=dims, nfeatures=nfeatures,
                                BPPARAM=BPPARAM, verbose=FALSE)
      }
    }
    nc <- length(unique(clusters))
    if(nc==1) stop("Only one cluster generated. Consider specifying `cluster` ",
                   "(e.g. `cluster=10`)")
    if(verbose) message(nc, " clusters")
  }else{
    characterize <- FALSE
  }
  cl <- clusters

  ## feature selection
  if(length(sel_features)>nfeatures)
    sel_features <- selFeatures(sce[sel_features,], cl, nfeatures=nfeatures,
                                propMarkers=propMarkers)
  sce <- sce[sel_features,]
  if(length(wDbl)>0) sce.dbl <- sce.dbl[sel_features,]

  ## get the artificial doublets
  if(is.null(artificialDoublets))
    artificialDoublets <- min( 25000, max(1500,
                                          ceiling(ncol(sce)*0.8),
                                          10*length(unique(cl))^2 ) )
  if(artificialDoublets<=2)
    artificialDoublets <- min(ceiling(artificialDoublets*ncol(sce)),25000)

  if(verbose)
    message("Creating ~", artificialDoublets, " artificial doublets...")
  ad <- getArtificialDoublets(counts(sce), n=artificialDoublets,
                              clusters=clusters, propRandom=propRandom, ...)

  gc(verbose=FALSE)

  ado <- ad$origins
  ad <- ad$counts

  no <- ncol(sce) + length(wDbl)
  ado2 <- as.factor(c(rep(NA, no), as.character(ado)))
  src <- factor( rep(1:2, c(no,ncol(ad))), labels = c("real","artificial"))
  ctype <- factor( rep(c(1L,ifelse(knownUse=="positive",2L,1L),2L),
                       c(ncol(sce),length(wDbl),ncol(ad))),
                   labels=c("real","doublet") )
  inclInTrain <- rep(c(TRUE,ifelse(knownUse=="positive",TRUE,FALSE),TRUE),
                     c(ncol(sce),length(wDbl),ncol(ad)))

  e <- counts(sce)
  if(!is.null(wDbl)) e <- cbind(e, counts(sce.dbl))
  e <- cbind(e, ad[row.names(sce),])

  # evaluate by library size and non-zero features
  lsizes <- Matrix::colSums(e)
  cxds_score <- cxds2(e, whichDbls=which(ctype=="doublet" | !inclInTrain))
  nfeatures <- Matrix::colSums(e>0L)
  nAbove2 <- Matrix::colSums(e>2L)

  if(returnType=="counts"){
    sce_out <- SingleCellExperiment(list(
      counts=cbind(counts(sce), ad[row.names(sce),])))
    sce_out$type <- ctype
    sce_out$src <- src
    sce_out$origin <- ado2
    sce_out$cluster <- NA
    if(!is.null(clusters)) colData(sce_out)[colnames(sce),"cluster"] <- clusters
    sce_out$cxds_score <- cxds_score
    return(sce_out)
  }

  if(verbose) message("Dimensional reduction")

  if(!is.null(clustCor) && !is.null(clusters)){
    if(!is.null(dim(clustCor))){
      clustCor <- .clustSpearman(e, clustCor)
    }else{
      clustCor <- .clustSpearman(e, clusters, nMarkers=clustCor)
    }
  }

  if(is.character(processing)){
    pca <- switch(processing,
                  default=.defaultProcessing(e, dims=dims),
                  rawPCA=.defaultProcessing(e, dims=dims, doNorm=FALSE),
                  rawFeatures=t(e),
                  atac=.atacProcessing(e, dims=dims),
                  normFeatures=t(normalizeCounts(e)),
                  stop("Unknown processing function.")
    )
  }else{
    pca <- processing(e, dims=dims)
    stopifnot(identical(row.names(pca),colnames(e)))
  }

  ex <- NULL
  if(!is.null(clusters)) ex <- getExpectedDoublets(clusters, dbr)

  if(verbose) message("Evaluating kNN...")
  d <- .evaluateKNN(pca, ctype, ado2, expected=ex, k=k)

  #if(characterize) knn <- d$knn   ## experimental
  d <- d$d
  if(!is.null(clusters)){
    d$cluster <- NA
    d[colnames(sce),"cluster"] <- clusters
  }else{
    d$cluster <- NULL
  }
  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$nAbove2 <- nAbove2
  d$src <- src
  d$cxds_score <- cxds_score
  d$include.in.training <- inclInTrain

  if(!is.null(clustCor)) d <- cbind(d, clustCor)

  ## classify
  if(length(includePCs)==1) includePCs <- seq_len(includePCs)
  includePCs <- includePCs[includePCs<ncol(pca)]
  d <- .scDblscore(d, scoreType=score, addVals=pca[,includePCs,drop=FALSE],
                   threshold=threshold, dbr=dbr, dbr.sd=dbr.sd, nrounds=nrounds,
                   max_depth=max_depth, iter=iter, BPPARAM=BPPARAM,
                   features=trainingFeatures, verbose=verbose, metric=metric,
                   filterUnidentifiable=removeUnidentifiable,
                   unident.th=unident.th)

  #if(characterize) d <- .callDblType(d, pca, knn=knn, origins=ado2)
  if(returnType=="table") return(d)
  if(returnType=="full"){
    sce_out <- SingleCellExperiment(list(
      counts=cbind(counts(sce), ad[row.names(sce),])), colData=d)
    reducedDim(sce_out, "PCA") <- pca
    if(is(d,"DataFrame") && !is.null(metadata(d)$scDblFinder.stats))
        metadata(sce_out)$scDblFinder.stats <- metadata(d)$scDblFinder.stats
    return(sce_out)
  }
  rowData(orig)$scDblFinder.selected <- row.names(orig) %in% sel_features
  .scDblAddCD(orig, d)
}

#' @importFrom BiocNeighbors AnnoyParam
.evaluateKNN <- function(pca, ctype, origins, expected=NULL, k){
  knn <- suppressWarnings(findKNN(pca, max(k), BNPARAM=AnnoyParam()))
  hasOrigins <- length(unique(origins))>1
  knn$type <- matrix(as.integer(ctype)[knn$index]-1L, nrow=nrow(knn$index))
  if(hasOrigins) knn$orig <- matrix(origins[knn$index], nrow=nrow(knn[[1]]))
  if(any(w <- knn$distance==0))
    knn$distance[w] <- min(knn$distance[knn$distance[,1]>0,1])

  md <- max(knn$distance[,1])
  dr <- t(vapply(seq_len(nrow(knn$distance)), FUN.VALUE=numeric(2L),
                 FUN=function(x){
                   w <- knn$type[x,]==1
                   dA <- ifelse(length(wA <- which(w))==0, 2*md,
                                knn$distance[x,wA[1]])
                   dB <- ifelse(length(wB <- which(!w))==0, 2*md,
                                knn$distance[x,wB[1]])
                   c(dA,dB)
                 }))
  dw <- sqrt(max(k)-seq_len(max(k))) * 1/knn$distance
  dw <- dw/rowSums(dw)
  d <- data.frame( row.names=row.names(pca), type=ctype, cluster=NA,
                   weighted=rowSums(knn$type*dw),
                   distanceToNearest=knn$distance[,1],
                   distanceToNearestDoublet=dr[,1],
                   distanceToNearestReal=dr[,2],
                   nearestClass=knn$type[,1] )
  if(hasOrigins) d <- cbind(d, .getMostLikelyOrigins(knn, origins))

  for(ki in k)
    d[[paste0("ratio.k",ki)]] <- rowSums(knn$type[,seq_len(ki)])/ki

  if(hasOrigins && !is.null(expected)){
    w <- which(d$type=="doublet")
    class.weighted <- vapply( split(d$weighted[w], d$mostLikelyOrigin[w]),
                              FUN.VALUE=numeric(1L), FUN=mean )

    d$difficulty <- 1
    w <- which(!is.na(d$mostLikelyOrigin))
    d$difficulty[w] <- 1-class.weighted[d$mostLikelyOrigin[w]]
    #d$difficulty <- .knnSmooth(knn, d$difficulty, use.distance=FALSE)

    d$expected <- expected[d$mostLikelyOrigin]
    ob <- table(d$mostLikelyOrigin)
    d$observed <- ob[d$mostLikelyOrigin]
    w <- which(is.na(d$mostLikelyOrigin))
    d$observed[w] <- d$expected[w] <- 0
  }
  list(knn=knn, d=d)
}

#' @importFrom stats quantile weighted.mean
.knnSmooth <- function(knn, score, use.distance=TRUE, type=NULL){
  w <- seq_len(ncol(knn$index))
  if(use.distance){
    mind <- quantile(knn$distance[,1], probs=0.1)
    if(mind==0) mind <- 0.5
  }
  vapply(seq_len(nrow(knn$index)), FUN.VALUE=numeric(1L), FUN=function(i){
    x <- knn$index[i,]
    if(!is.null(type)){
      w <- knn$type[i,]==type
    }
    if(sum(w)==0) return(score[i])
    x <- x[w]
    if(use.distance){
      weights <- mind+c(0,knn$distance[i,][w])
      weights <- 1/sqrt(weights)
    }else{
      weights <- 1/seq_len(1+length(x))
    }
    weighted.mean(c(score[i],score[x]),weights)
  })
}

#' @importFrom S4Vectors DataFrame metadata
#' @importFrom stats predict quantile
.scDblscore <- function(d, scoreType="xgb", nrounds=NULL, max_depth=5, iter=2,
                        threshold=TRUE, verbose=TRUE, dbr=NULL, dbr.sd=NULL,
                        features=NULL, filterUnidentifiable=TRUE, addVals=NULL,
                        metric="logloss", eta=0.3, BPPARAM=SerialParam(),
                        includeSamples=FALSE, perSample=TRUE, unident.th=0.1, ...){
  gdbr <- .gdbr(d, dbr)
  if(!is.null(d$sample) && length(unique(d$sample))==1) d$sample <- NULL
  if(is.null(dbr.sd)) dbr.sd <- 0.3*gdbr+0.025
  if(scoreType=="xgb"){
    if(verbose) message("Training model...")
    d$score <- NULL
    if(is.null(features)){
      prds <- .defTrainFeatures(d)
    }else{
      if("ratio.k*" %in% features)
        features <- c(features[features!="ratio.k*"],
                      grep("^ratio\\.k",colnames(d),value=TRUE))
      toExclude <- grep("^-",features)
      if(length(toExclude)==0){
        if(length(mis <- setdiff(features, colnames(d)))>0)
          warning("The following features were not found: ",
                  paste(mis,collapse=", "))
        prds <- intersect(features, colnames(d))
      }else{
        if(length(toExclude)!=length(features))
          stop("Mixture of included/excluded features - use only either.")
        prds <- setdiff(.defTrainFeatures(d), gsub("^-","",features))
      }
      prds <- setdiff(prds,c("type","src","class","cluster"))
    }
    if(!is.null(features)) message(paste("Features used for training:\n",
                                   paste(prds,collapse=", ")))
    preds <- as(as.matrix(d[,prds,drop=FALSE]), "CsparseMatrix")


    if(includeSamples && !is.null(d$sample))
      preds <- cbind(preds, as(stats::model.matrix(~d$sample)[,-1,drop=FALSE],
                               "CsparseMatrix"))

    if(!is.null(addVals)){
      stopifnot(nrow(addVals)==nrow(preds))
      preds <- cbind(preds, as(addVals, "CsparseMatrix"))
      rm(addVals)
    }
    w <- which(d$type=="real")
    ratio <- rev(grep("^ratio\\.k",colnames(d)))[1]
    if(!is.null(d$sample) && !perSample){
      tt <- table(d$type, d$sample)
      expected.ratio <- (1+tt["doublet",])/(1+tt["real",])
      d$adjusted.ratio <- d[[ratio]]/expected.ratio[d$sample]
      d <- .rescaleSampleScores(d, TRUE, what="cxds_score",
                                newName="adjusted.cxds")
      d$score <- (d$adjusted.ratio + d$adjusted.cxds)/2
    }else{
      d$score <- (d$cxds_score + d[[ratio]]/max(d[[ratio]]))/2
    }
    max.iter <-  iter
    while(iter>0){
      # remove cells with a high chance of being doublets from the training,
      # as well as unidentifiable artificial doublets
        w <- which( (d$type=="real" &
          doubletThresholding(d, dbr=dbr, dbr.sd=dbr.sd, stringency=0.7,
                              perSample=perSample,
                              returnType="call")=="doublet") |
            (d$type=="doublet" & d$score<unident.th & filterUnidentifiable) |
            !d$include.in.training )
      if(verbose) message("iter=",max.iter-iter,", ", length(w),
                          " cells excluded from training.")
      d$score <- tryCatch({
        fit <- .xgbtrain(preds[-w,], d$type[-w], nrounds, metric=metric,
                         max_depth=max_depth, eta=eta, #base_score=gdbr,
                         nthreads=BiocParallel::bpnworkers(BPPARAM))
        predict(fit, as.matrix(preds))
      }, error=function(e) d$score)
      if(!is.null(d$mostLikelyOrigin)){
        wO <- which(d$type!="real" & !is.na(d$mostLikelyOrigin))
        class.diff <- vapply( split(d$score[wO], d$mostLikelyOrigin[wO]),
                              FUN.VALUE=numeric(1L), FUN=mean )
        d$difficulty <- mean(class.diff)
        wO <- which(!is.na(d$mostLikelyOrigin))
        d$difficulty[wO] <- 1-class.diff[d$mostLikelyOrigin[wO]]
        if(filterUnidentifiable && iter==max.iter)
          d <- .filterUnrecognizableDoublets(d)
      }
      iter <- iter-1
    }
    d$include.in.training[w] <- FALSE
    ########################
    # Uncomment and use with scDblFinder(..., returnType="table") to extract
    # variable importance
    # return(xgb.importance(model=fit))
    #######################
  }else{
    if(scoreType=="ratio"){
      d$score <- d$ratio
    }else{
      d$score <- d$weighted
    }
  }
  d <- DataFrame(d)
  if(threshold){
    th <- doubletThresholding( d, dbr=dbr, dbr.sd=dbr.sd, perSample=perSample,
                               ... )
    if(!is.null(d$sample) && length(th)>1){
      d$class <- ifelse(d$score >= th[d$sample], "doublet", "singlet")
    }else{
      d$class <- ifelse(d$score >= th, "doublet", "singlet")
    }
    if(verbose) message("Threshold found:", paste(round(th,3), collapse=" "))
    ## set class of known (i.e. inputted) doublets:
    d$class[d$src=="real" & d$type=="doublet"] <- "doublet"
    if(!is.null(d$mostLikelyOrigin)){
      th.stats <- .getDoubletStats(d, th, dbr, dbr.sd)
      metadata(d)$scDblFinder.stats <- th.stats
    }
    metadata(d)$scDblFinder.threshold <- th
    d$nearestClass <- factor(d$nearestClass, levels = 0:1,
                             labels=c("cell","artificialDoublet"))
    dbr <- sum(d$class=="doublet" & d$src=="real")/sum(d$src=="real")
    if(verbose) message(sum(d$class=="doublet" & d$src=="real"), " (",
                        round(100*dbr,1),"%) doublets called")
  }
  d
}

.defTrainFeatures <- function(d){
  setdiff(colnames(d), c("mostLikelyOrigin","originAmbiguous",
                         "distanceToNearestDoublet", "type",
                         "src","distanceToNearest","class",
                         "nearestClass","cluster","sample","expected",
                         "include.in.training","observed"))
}

#' @importFrom xgboost xgb.cv xgboost
.xgbtrain <- function(d2, ctype, nrounds=NULL, max_depth=6, nfold=5,
                      tree_method="exact", subsample=0.75, nthreads=1,
                      metric="logloss", ...){
  if(!is.integer(ctype)) ctype <- as.integer(ctype)-1
  d2 <- as.matrix(d2)
  if(is.null(nrounds)) nrounds <- 0L
  if(!is.numeric(nrounds) || nrounds<0)
    stop("If given, `nrounds` must be a positive number!")
  if(nrounds<=1){
    # use cross-validation
    res <- xgb.cv(data=d2, label=ctype, nrounds=200, max_depth=max_depth,
                  objective="binary:logistic", eval_metric=metric,
                  early_stopping_rounds=2, tree_method=tree_method, nfold=nfold,
                  subsample=subsample, nthread=nthreads, verbose=FALSE, ...)
    best <- res$best_iteration
    if(nrounds==0){
      nrounds <- best
    }else{
      e <- res$evaluation_log
      ac <- e[[grep("test.+mean",colnames(e))]][best] +
        nrounds * e[[grep("test.+std",colnames(e))]][best]
      nrounds <- min(which(e[[grep("test.+mean",colnames(e))]] <= ac))
    }
    #message("Best iteration: ", best, "; selected nrounds: ", nrounds)
  }
  xgboost( d2, ctype, nrounds=nrounds, eval_metric=metric,
           objective="binary:logistic", tree_method=tree_method,
           max_depth=max_depth, early_stopping_rounds=2, verbose=FALSE,
           nthread=nthreads, ... )
}

.aggResultsTable <- function(d, keep.col=NULL){
  ths <- sapply(d, FUN=function(x){
    if(!is(x,"DFrame") || is.null(th <- metadata(x)$scDblFinder.threshold))
      return(NULL)
    th
  })
  cn <- table(unlist(lapply(d, colnames)))
  cn <- c(names(cn)[cn==length(d)], "total.prop.real")
  if(!is.null(keep.col)) cn <- intersect(cn, keep.col)
  d <- do.call(rbind, lapply(d, FUN=function(x){
    x$total.prop.real <- sum(x$type=="real",na.rm=TRUE)/nrow(x)
    if(!is.null(x$cluster)) x$cluster <- as.character(x$cluster)
    x[,cn]
  }))
  if(!is.null(d$cluster)) d$cluster <- as.factor(d$cluster)
  if(!is.null(d$sample)) d$sample <- as.factor(d$sample)
  metadata(d)$scDblFinder.threshold <- ths
  d
}

# add the relevant fields of the scDblFinder results table to the SCE
#' @importFrom stats relevel
.scDblAddCD <- function(sce, d){
  fields <- c("sample","cluster","class","score","ratio","weighted",
              "difficulty","cxds_score","mostLikelyOrigin","originAmbiguous",
              "origin.prob", "origin.call", "origin.2ndBest")
  if(!is.data.frame(d) && is.list(d)) d <- .aggResultsTable(d, fields)
  d <- d[colnames(sce),]
  for(f in fields){
    if(!is.null(d[[f]])) sce[[paste0("scDblFinder.",f)]] <- d[[f]]
  }
  if(!is.null(sce$scDblFinder.class)) sce$scDblFinder.class <-
      relevel(as.factor(sce$scDblFinder.class),"singlet")
  if(is(d,"DataFrame")){
    if(!is.null(metadata(d)$scDblFinder.stats))
      metadata(sce)$scDblFinder.stats <- metadata(d)$scDblFinder.stats
    metadata(sce)$scDblFinder.threshold <- metadata(d)$scDblFinder.threshold
  }
  sce
}


## sets a reasonable set of ks (for KNN)
.defaultKnnKs <- function(k=NULL, n){
  if(!is.null(dim(n))) n <- ncol(n)
  if(!is.null(k)) return(k[k<=ceiling(n/2)])
  kmax <- max(ceiling(sqrt(n/2)),25)
  k <- c(3,10,15,20,25,50,kmax)
  unique(k[k<=kmax])
}
