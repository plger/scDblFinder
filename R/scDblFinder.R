#' scDblFinder
#' 
#' Identification of doublets in single-cell RNAseq directly from counts using 
#' overclustering-based generation of artifical doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param artificialDoublets The approximate number of artificial doublets to 
#' create. If NULL, will be the maximum of the number of cells or 
#' `5*nbClusters^2`.
#' @param clusters The optional cluster assignments (if omitted, will run 
#' clustering). This is used to make doublets more efficiently. `clusters` 
#' should either be a vector of labels for each cell, or the name of a colData 
#' column of `sce`.
#' @param clust.method The clustering method if `clusters` is not given.
#' @param samples A vector of the same length as cells (or the name of a column 
#' of `colData(sce)`), indicating to which sample each cell belongs. Here, a 
#' sample is understood as being processed independently. If omitted, doublets 
#' will be searched for with all cells together. If given, doublets will be 
#' searched for independently for each sample, which is preferable if they 
#' represent different captures. If your samples were multiplexed using cell
#' hashes, want you want to give here are the different batches/wells (i.e. 
#' independent captures, since doublets cannot arise across them) rather
#' than biological samples.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if `clusters` is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if 
#' `clusters` is given or if `clust.method`!='overcluster'.
#' @param nfeatures The number of top features to use (default 1000)
#' @param dims The number of dimensions used to build the network (default 20)
#' @param dbr The expected doublet rate. By default this is assumed to be 1\% 
#' per thousand cells captured (so 4\% among 4000 thousand cells), which is 
#' appropriate for 10x datasets.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.015.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param fullTable Logical; whether to return the full table including 
#' artificial doublets, rather than the table for real cells only (default).
#' @param score Score to use for final classification.
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when 
#' `samples!=NULL`); otherwise passed to eventual PCA and K/SNN calculations.
#'
#' @return The `sce` object with several additional colData columns, in 
#' particular `scDblFinder.score` (the final score used) and `scDblFinder.class` 
#' (whether the cell is called as 'doublet' or 'singlet'). Alternatively, if 
#' `fullTable=TRUE`, a data.frame will be returned with information about real 
#' and artificial cells.
#' 
#' @import SingleCellExperiment Matrix BiocParallel xgboost
#' @importFrom SummarizedExperiment colData<- assayNames
#' @importFrom scater normalizeCounts runPCA
#' @importFrom dplyr bind_rows
#' @importFrom methods is
#' @importFrom DelayedArray as.matrix
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocSingular IrlbaParam
#' 
#' @examples
#' library(SingleCellExperiment)
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' sce <- SingleCellExperiment( list(counts=m) )
#' sce <- scDblFinder(sce, verbose=FALSE)
#' table(sce$scDblFinder.class)
#' 
#' @export
scDblFinder <- function( sce, artificialDoublets=NULL, clusters=NULL,
                         clust.method=c("fastcluster","overcluster"),
                         samples=NULL, minClusSize=min(50,ncol(sce)/5), 
                         maxClusSize=NULL, nfeatures=1000, dims=20, dbr=NULL, 
                         dbr.sd=0.015, k=20, fullTable=FALSE,
                         score=c("xgb.local.optim","xgb","weighted","ratio"),
                         verbose=is.null(samples), BPPARAM=SerialParam()
                        ){
  if(!is(sce, "SingleCellExperiment"))
      stop("`sce` should be a SingleCellExperiment")
  clust.method <- match.arg(clust.method)
  if(!is.null(clusters) && is.character(clusters) && length(clusters)==1){
      if(!(clusters %in% colnames(colData(sce)))) 
          stop("Could not find `clusters` column in colData!")
      clusters <- colData(sce)[[clusters]]
  }
  if(!is.null(score)) score <- match.arg(score)
  stopifnot(is(sce, "SingleCellExperiment"))
  if( !("counts" %in% assayNames(sce)) ) 
      stop("`sce` should have an assay named 'counts'")
  if(!is.null(samples)){
    # splitting by samples
    if(length(samples)==1 && samples %in% colnames(colData(sce)))
        samples <- colData(sce)[[samples]]
    if(ncol(sce)!=length(samples))
        stop("Samples vector matches number of cells")
    if(fullTable) warning("`fullTable` param ignored when splitting by samples")
    if(verbose) message("`verbose` param ignored when splitting by samples.")
    cs <- split(seq_along(samples), samples, drop=TRUE)
    tmp <- bplapply(cs, BPPARAM=BPPARAM, FUN=function(x){ 
        if(!is.null(clusters) && length(clusters)>1) clusters <- clusters[x]
        x <- colData( scDblFinder(sce[,x], artificialDoublets=artificialDoublets, 
                                  clusters=clusters[x], minClusSize=minClusSize, 
                                  maxClusSize=maxClusSize, dims=dims, dbr=dbr, 
                                  dbr.sd=dbr.sd, k=k, clust.method=clust.method,
                                  score=score, nfeatures=nfeatures, 
                                  verbose=FALSE ) )
        fields <- paste0("scDblFinder.",c("weighted","ratio","class","score"))
        as.data.frame(x[,fields])
    })
    CD <- bind_rows(tmp)
    row.names(CD) <- unlist(lapply(tmp, row.names))
    for(f in colnames(CD)) colData(sce)[[f]] <- CD[unlist(cs),f]
    return(sce)
  }
  if(ncol(sce)<100) warning("scDblFinder might not work well with very low 
numbers of cells.")

  if(is.null(dbr)){
      ## dbr estimated as for chromium data, 1% per 1000 cells captured:
      dbr <- (0.01*ncol(sce)/1000)
  }
  
  orig <- sce
  if(nrow(sce)>nfeatures){
    if(verbose) message("Identifying top genes...")
    sce <- sce[.clusterTopG(sce,nfeatures=nfeatures,clusters),]
  }  
  if(is.null(clusters)){
      if(verbose) message("Clustering cells...")
      if(clust.method=="overcluster"){
        clusters <- overcluster(sce, min.size=minClusSize, max.size=maxClusSize,
                                ndims=ndims)
      }else{
        clusters <- fastcluster(sce, ndims=ndims)
      }
  }
  if(length(unique(clusters)) == 1) stop("Only one cluster generated")

  maxSameDoublets <- length(clusters)*(dbr+2*dbr.sd) * 
      prod(sort(table(clusters), decreasing=TRUE)[1:2] / length(clusters))
  if(min(table(clusters)) <= maxSameDoublets){
    warning("In light of the expected rate of doublets given, and of the size ",
            "of the clusters, it is possible that some of the smaller clusters",
            " are composed of doublets of the same type.")
  }
  cli <- split(seq_len(ncol(sce)), clusters)
  
  if(is.null(colnames(sce)))
      colnames(orig) <- colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  if(is.null(row.names(sce)))
      row.names(sce) <- paste0("f",seq_len(nrow(sce)))

  # get the artificial doublets
  if(is.null(artificialDoublets)){
    artificialDoublets <- max( ncol(sce), 
                               min(40000, 5*length(unique(clusters))^2))
  }
  
  if(verbose){
    message("Creating ~", artificialDoublets, " artifical doublets...")
    ad <- getArtificialDoublets( as.matrix(counts(sce)), n=artificialDoublets, 
                                 clusters=clusters )
  }else{
    ad <- suppressWarnings( getArtificialDoublets(as.matrix(counts(sce)), 
                                                  n=artificialDoublets, 
                                                  clusters=clusters ) )
  }
  ado <- ad$origins
  ad <- ad$counts
  ado2 <- as.factor(c(rep(NA, ncol(sce)), as.character(ado)))

  e <- cbind(as.matrix(counts(sce)), ad[row.names(sce),])
  e <- normalizeCounts(e)
  pca <- tryCatch({
            scater::calculatePCA(e, dims, subset_row=seq_len(nrow(e)),
                                 BSPARAM=BiocSingular::IrlbaParam())
        }, error=function(msg){
            reducedDim( scater::runPCA( SingleCellExperiment(list(logcounts=e)), 
                                        ncomponents=dims, ntop=nrow(e),
                                        BSPARAM=BiocSingular::IrlbaParam()) )
        })
  if(is.list(pca)) pca <- pca$x
  row.names(pca) <- colnames(e)

  # evaluate by library size and non-zero features
  lsizes <- c(colSums(counts(sce)),colSums(ad))
  nfeatures <- c(colSums(counts(sce)>0), colSums(ad>0))
  
  
  ctype <- factor( rep(1:2, c(ncol(sce),ncol(ad))), 
                   labels = c("real","artificial"))
  knn <- .evaluateKNN(pca, ctype, ado2, k=k, dbr=dbr, BPPARAM=BPPARAM, 
                      verbose=verbose)
  d <- knn$d
  knn <- knn$knn
  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$type <- ctype

  if(is.null(score)) return(d)

  if(score %in% c("xgb.local.optim","xgb")){
    if(verbose) message("Training model...")
    prds <- setdiff(colnames(d), c("mostLikelyOrigin","originAmbiguous","type",
                                   "distanceToNearest"))
    d2  <- as.matrix(d[,prds])
    nd <- round(sum(ex))+sum(d$type=="artificial")
    fit <- xgboost(d2, as.integer(d$type)-1, nrounds=50, 
                   max_depth=6, objective="binary:logistic", 
                   #scale_pos_weight=sqrt(nd/(nrow(d)-nd)),
                   early_stopping_rounds=5, verbose=FALSE )
    d$score <- predict(fit, d2)
  }else{
    if(score=="ratio"){
      d$score <- d$ratio
    }else{
      d$score <- d$weighted
    }
  }
  
  d$cluster <- NA
  d$cluster[d$type=="real"] <- clusters
  d$smoothed.score <- .knnSmooth(knn, d$score, type=0)
  rm(knn)
  
  if(verbose) message("Finding threshold...")
  th <- doubletThresholding( d, dbr=dbr, dbr.sd=dbr.sd, 
                             local=score=="xgb.local.optim" )
  if(score=="xgb.local.optim"){
    d$score.global <- d$score
    d$score <- th$finalScores
  }
  metadata(orig)$scDblFinder.stats <- th$stats
  th <- th$th
  if(verbose) message("Threshold found:", round(th,3))

  d$nearestClass <- factor(d$nearestClass, levels = 0:1, 
                           labels=c("cell","artificialDoublet"))
  d$class <- ifelse(d$score >= th, "doublet", "singlet")
  if(fullTable) return(d)
  
  d <- d[seq_len(ncol(orig)),]
  d$type <- NULL
  row.names(d) <- colnames(orig)
  d <- d[colnames(orig),]
  
  orig$scDblFinder.cluster <- clusters
  for(f in c("cluster","distanceToNearest","nearestClass","difficulty","ratio",
             "weighted","score.global","score","smoothed.score","class",
             "mostLikelyOrigin","originAmbiguous")){
    if(!is.null(d[[f]])) orig[[paste0("scDblFinder.",f)]] <- d[[f]]
  }
  orig
}

.evaluateKNN <- function(pca, ctype, origins, k, dbr, BPPARAM=SerialParam(), 
                         verbose=TRUE){
  if(verbose) message("Finding KNN...")
  knn <- suppressWarnings(findKNN(pca, k, BPPARAM=BPPARAM))
  
  if(verbose) message("Evaluating cell neighborhoods...")
  knn$type <- matrix(as.numeric(ctype)[knn$index]-1, nrow=nrow(knn$index))
  knn$orig <- matrix(origins[knn$index], nrow=nrow(knn[[1]]))
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
  
  dw <- 1/knn$distance
  dw <- dw/rowSums(dw)
  d <- data.frame( type=ctype,
                   weighted=rowSums(knn$type*dw),
                   distanceToNearest=knn$distance[,1],
                   distanceToNearestDoublet=dr[,1],
                   distanceToNearestReal=dr[,2],
                   nearestClass=knn$type[,1],
                   ratio=rowSums(knn$type)/k,
                   .getMostLikelyOrigins(knn, origins) )
  w <- which(d$type=="artificial")
  class.weighted <- vapply( split(d$weighted[w], d$mostLikelyOrigin[w]), 
                            FUN.VALUE=numeric(1L), FUN=mean )
  d$difficulty <- 1
  w <- which(!is.na(d$mostLikelyOrigin))
  d$difficulty[w] <- 1-class.weighted[d$mostLikelyOrigin[w]]
  d$difficulty <- .knnSmooth(knn, d$difficulty)
  
  ex <- getExpectedDoublets(clusters, dbr)
  d$expected <- ex[d$mostLikelyOrigin]
  ob <- table(d$mostLikelyOrigin)
  d$observed <- ob[d$mostLikelyOrigin]
  w <- which(is.na(d$mostLikelyOrigin))
  d$observed[w] <- d$expected[w] <- 0
  list(knn=knn, d=d)
}

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