#' scDblFinder
#' 
#' Identification of heterotypic (or neotypic) doublets in single-cell RNAseq 
#' using cluster-based generation of artifical doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment-class}}, or array of
#' counts.
#' @param artificialDoublets The approximate number of artificial doublets to 
#' create. If \code{NULL}, will be the maximum of the number of cells or 
#' \code{5*nbClusters^2}.
#' @param clusters The optional cluster assignments (if omitted, will run 
#' clustering). This is used to make doublets more efficiently. \code{clusters} 
#' should either be a vector of labels for each cell, or the name of a colData 
#' column of \code{sce}.
#' @param clust.method The clustering method if \code{clusters} is not given.
#' @param samples A vector of the same length as cells (or the name of a column 
#' of \code{colData(x)}), indicating to which sample each cell belongs. Here, a 
#' sample is understood as being processed independently. If omitted, doublets 
#' will be searched for with all cells together. If given, doublets will be 
#' searched for independently for each sample, which is preferable if they 
#' represent different captures. If your samples were multiplexed using cell
#' hashes, want you want to give here are the different batches/wells (i.e. 
#' independent captures, since doublets cannot arise across them) rather
#' than biological samples.
#' @param knownDoublets An optional logical vector of known doublets (e.g. 
#' through cell barcodes), or the name of a colData column of `sce` containing
#' that information.
#' @param use.cxds Logical; whether to use `scds::cxds` scores in addition to 
#' information from artificial/known doublets as part of the predictors.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if \code{clusters} is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if 
#' \code{clusters} is given or if \code{clust.method!='overcluster'}.
#' @param nfeatures The number of top features to use (default 1000)
#' @param dims The number of dimensions used to build the network (default 20)
#' @param dbr The expected doublet rate. By default this is assumed to be 1\% 
#' per thousand cells captured (so 4\% among 4000 thousand cells), which is 
#' appropriate for 10x datasets.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.015.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param returnType Either "sce" (default), "table" (to return the table of 
#' cell attributes including artificial doublets), or "full" (returns an SCE
#' object containing both the real and artificial cells.
#' @param score Score to use for final classification.
#' @param threshold Logical; whether to threshold scores into binary doublet 
#' calls
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when 
#' `samples!=NULL`); otherwise passed to eventual PCA and K/SNN calculations.
#'
#' @return The \code{sce} object with several additional colData columns, in 
#' particular `scDblFinder.score` (the final score used) and `scDblFinder.class` 
#' (whether the cell is called as 'doublet' or 'singlet'). See 
#' \code{vignette("scDblFinder")} for more details; for alternative return 
#' values, see the `returnType` argument.
#' 
#' @details
#' This function generates artificial doublets from clusters of real cells, 
#' evaluates their prevalence in the neighborhood of each cells, and uses this 
#' along with additional features to classify doublets. The approach is 
#' complementary to doublets identified via cell hashes and SNPs in multiplexed
#' samples: the latter can identify doublets formed by cells of the same type
#' from two samples, which are nearly undistinguishable from real cells 
#' transcriptionally, but cannot identify doublets made by cells of the 
#' same sample. See \code{vignette("scDblFinder")} for more details on the 
#' method.
#' 
#' When inter-sample doublets are available, they can be provided to 
#' `scDblFinder` through the \code{knownDoublets} argument to improve the 
#' identification of further doublets.
#' 
#' When multiple samples/captures are present, they can be specified using the 
#' \code{samples} argument and will be processed separately (in parallel if 
#' \code{BPPARAM} is given).
#' 
#' @import SingleCellExperiment BiocParallel xgboost
#' @importFrom SummarizedExperiment colData<- assayNames
#' @importFrom scuttle normalizeCounts
#' @importFrom scater runPCA
#' @importFrom scds cxds
#' @importFrom methods is
#' @importFrom DelayedArray as.matrix
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocSingular IrlbaParam
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- mockDoubletSCE()
#' sce <- scDblFinder(sce, verbose=FALSE)
#' table(sce$scDblFinder.class)
#' 
#' @export
#' @rdname scDblFinder
#' @import SingleCellExperiment
scDblFinder <- function( sce, clusters=NULL, samples=NULL, 
                         artificialDoublets=NULL, knownDoublets=NULL,
                         clust.method=c("fastcluster","overcluster"),
                         use.cxds=TRUE, minClusSize=min(50,ncol(sce)/5),
                         maxClusSize=NULL, nfeatures=1000, dims=20, dbr=NULL, 
                         dbr.sd=0.015, k=20, returnType=c("sce","table","full"),
                         score=c("xgb","xgb.local.optim","weighted","ratio"),
                         threshold=TRUE, verbose=is.null(samples), 
                         BPPARAM=SerialParam()
                        ){
  sce <- .checkSCE(sce)
  score <- match.arg(score)
  clust.method <- match.arg(clust.method)
  returnType <- match.arg(returnType)
  clusters <- .checkColArg(sce, clusters)
  knownDoublets <- .checkColArg(sce, knownDoublets)
  samples <- .checkColArg(sce, samples)
  

  if(!is.null(samples)){
    # splitting by samples
    if(returnType=="full") 
        warning("`returnType='full'` ignored when splitting by samples")
    cs <- split(seq_along(samples), samples, drop=TRUE)
    d <- bplapply(cs, BPPARAM=BPPARAM, FUN=function(x){ 
        if(!is.null(clusters) && length(clusters)>1) clusters <- clusters[x]
        scDblFinder(sce[,x], artificialDoublets=artificialDoublets, 
                    clusters=clusters[x], minClusSize=minClusSize, 
                    maxClusSize=maxClusSize, dims=dims, dbr=dbr, 
                    dbr.sd=dbr.sd, k=k, clust.method=clust.method,
                    score="weighted", nfeatures=nfeatures, returnType="table",
                    knownDoublets=knownDoublets, threshold=FALSE, verbose=FALSE)
    })
    ss <- factor(rep(seq_along(names(d)),vapply(d,nrow,integer(1))), 
                 levels=seq_along(names(d)), labels=names(d))
    d <- do.call(rbind, d)
    d$sample <- ss
    d <- .scDblscore(d, scoreType=score, threshold=FALSE, verbose=verbose,
                     dbr=dbr, dbr.sd=dbr.sd)
    if(returnType=="table") return(d)
    return(.scDblAddCD(sce, d))
  }
  if(ncol(sce)<100)
    warning("scDblFinder might not work well with very low numbers of cells.")
  if(ncol(sce)>25000)
    warning("You are trying to run scDblFinder on a very large number of ",
            "cells. If these are from different captures, please specify this",
            " using the `samples` argument.", immediate=TRUE)

  if(is.null(dbr)){
      ## dbr estimated as for chromium data, 1% per 1000 cells captured:
      dbr <- (0.01*ncol(sce)/1000)
  }
  
  orig <- sce
  if(nrow(sce)>nfeatures){
    if(verbose) message("Identifying top genes...")
    sce <- sce[.clusterTopG(sce,nfeatures=nfeatures,clusters),]
  }

  wDbl <- c()
  if(!is.null(knownDoublets) && length(wDbl <- which(knownDoublets))>0){
    sce$knownDoublet <- knownDoublets
    sce.dbl <- sce[,wDbl,drop=FALSE]
    sce <- sce[,-wDbl,drop=FALSE]
    if(!is.null(clusters)){
      clusters.dbl <- clusters[wDbl]
      clusters <- clusters[-wDbl]
      if(is.factor(clusters)) clusters <- droplevels(clusters)
    }
  }
  if(is.null(clusters)){
      if(verbose) message("Clustering cells...")
      if(clust.method=="overcluster"){
        clusters <- overcluster(sce, min.size=minClusSize, max.size=maxClusSize,
                                ndims=dims)
      }else{
        clusters <- fastcluster(sce, ndims=dims)
      }
  }
  if((nc <- length(unique(clusters))) == 1) stop("Only one cluster generated")
  if(verbose) message(nc, " clusters")

  maxSameDoublets <- length(clusters)*(dbr+2*dbr.sd) * 
      prod(sort(table(clusters), decreasing=TRUE)[1:2] / length(clusters))
  if(min(table(clusters)) <= maxSameDoublets){
    warning("In light of the expected rate of doublets given, and of the size ",
            "of the clusters, it is possible that some of the smaller clusters",
            " are composed of doublets of the same type.")
  }

  # get the artificial doublets
  if(is.null(artificialDoublets))
    artificialDoublets <- min(max(ncol(sce), 5*length(unique(clusters))^2),
                              25000)
  
  if(verbose){
    message("Creating ~", artificialDoublets, " artifical doublets...")
    ad <- getArtificialDoublets( as.matrix(counts(sce)), n=artificialDoublets, 
                                 clusters=clusters )
  }else{
    ad <- suppressWarnings( getArtificialDoublets(as.matrix(counts(sce)), 
                                                  n=artificialDoublets, 
                                                  clusters=clusters ) )
  }
  gc(verbose=FALSE)
  
  ado <- ad$origins
  ad <- ad$counts
  
  no <- ncol(sce) + length(wDbl)
  ado2 <- as.factor(c(rep(NA, no), as.character(ado)))
  src <- factor( rep(1:2, c(no,ncol(ad))), labels = c("real","artificial"))
  ctype <- factor( rep(c(1,2,2), c(ncol(sce),length(wDbl),ncol(ad))), 
                   labels=c("real","doublet") )
  
  if(verbose) message("Dimensional reduction")
  
  e <- counts(sce)
  if(!is.null(wDbl)) e <- cbind(e, counts(sce.dbl))
  e <- cbind(e, ad[row.names(sce),])
  
  # evaluate by library size and non-zero features
  lsizes <- Matrix::colSums(e)
  cxds_score <- NULL
  if(use.cxds)
    cxds_score <- scds::cxds(SingleCellExperiment(list(counts=e)))$cxds_score
  nfeatures <- Matrix::colSums(e>0)
  
  # skip normalization if data is too large
  if(ncol(e)<=25000) e <- normalizeCounts(e)
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
  
  ex <- getExpectedDoublets(clusters, dbr)
  d <- .evaluateKNN(pca, ctype, ado2, expected=ex, k=k, BPPARAM=BPPARAM, 
                      verbose=verbose)$d
  d[colnames(sce),"cluster"] <- clusters
  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$src <- src
  if(use.cxds) d$cxds_score <- cxds_score
  
  d <- .scDblscore(cbind(d, pca[,1:2]), scoreType=score, 
                   threshold=threshold, verbose=verbose, 
                   dbr=dbr, dbr.sd=dbr.sd)
  if(returnType=="table") return(d)
  if(returnType=="full"){
      sce_out <- SingleCellExperiment(list(
          counts=cbind(counts(sce), ad[row.names(sce),])), colData=d)
      reducedDim(sce_out, "PCA") <- pca
      if(is(d,"DataFrame") && !is.null(metadata(d)$scDblFinder.stats)) 
        metadata(sce_out)$scDblFinder.stats <- metadata(d)$scDblFinder.stats
      return(sce_out)
  }
  .scDblAddCD(orig, d)
}

.evaluateKNN <- function(pca, ctype, origins, expected, k, 
                         BPPARAM=SerialParam(), verbose=TRUE){
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
  d <- data.frame( row.names=row.names(pca), type=ctype, cluster=NA, 
                   weighted=rowSums(knn$type*dw),
                   distanceToNearest=knn$distance[,1],
                   distanceToNearestDoublet=dr[,1],
                   distanceToNearestReal=dr[,2],
                   nearestClass=knn$type[,1],
                   ratio=rowSums(knn$type)/k,
                   .getMostLikelyOrigins(knn, origins) )
  w <- which(d$type=="doublet")
  class.weighted <- vapply( split(d$weighted[w], d$mostLikelyOrigin[w]), 
                            FUN.VALUE=numeric(1L), FUN=mean )
  
  d$difficulty <- 1
  w <- which(!is.na(d$mostLikelyOrigin))
  d$difficulty[w] <- 1-class.weighted[d$mostLikelyOrigin[w]]
  d$difficulty <- .knnSmooth(knn, d$difficulty)
  
  d$expected <- expected[d$mostLikelyOrigin]
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

#' @import xgboost
#' @importFrom S4Vectors DataFrame metadata
.scDblscore <- function(d, scoreType="xgb", nrounds=NULL, threshold=TRUE, 
                        verbose=TRUE, dbr=NULL, ...){
    if(scoreType %in% c("xgb.local.optim","xgb")){
        if(verbose) message("Training model...")
        d$score <- NULL
        prds <- setdiff(colnames(d), c("mostLikelyOrigin","originAmbiguous",
                                       "type","src","distanceToNearest","class",
                                       "nearestClass","cluster","sample"))
        d2 <- d[,prds]
        
        # remove cells with a high chance of being doublets from the training
        w <- which(d$type=="real" & d$ratio>=0.8)
        if(!is.null(d$cxds_score))
            w <- union(w, which(d$type=="real" & 
                                d$cxds_score>=quantile(d$cxds_score, 
                                                       prob=1-(dbr/2))))
        if(length(w)>0){
            ctype <- d$type[-w]
            d2 <- d2[-w,]
        }
        ctype <- as.integer(ctype)-1
        d2 <- as.matrix(d2)
        if(is.null(nrounds)){
            # use cross-validation
            res <- xgb.cv(data=d2, label=ctype, nrounds=300, 
                          objective="binary:logistic", nfold=5, 
                          early_stopping_rounds=3, tree_method="hist", 
                          metrics=list("error"), subsample=0.6, verbose=FALSE)
            ni = res$best_iteration
            ac = res$evaluation_log$test_error_mean[ni] + 
                1 * res$evaluation_log$test_error_std[ni]
            nrounds = min(which(res$evaluation_log$test_error_mean <= ac))
        }
        fit <- xgboost(d2, ctype, nrounds=nrounds, 
                       objective="binary:logistic", max_depth=6, 
                       early_stopping_rounds=5, verbose=FALSE )
        d$score <- predict(fit, as.matrix(d[,prds]))
    }else{
        if(scoreType=="ratio"){
            d$score <- d$ratio
        }else{
            d$score <- d$weighted
        }
    }
    d <- DataFrame(d)
    if(threshold){
        if(verbose) message("Finding threshold...")
        if(!is.null(d$sample) && is.null(dbr) && scoreType!="xgb.local.optim"){
            # per-sample thresholding
            th <- lapply(split(d[,c("cluster","src","type","mostLikelyOrigin",
                                    "originAmbiguous","score")], FUN=fun),
                         FUN=function(x){
                             dbr <- 0.01*sum(d$src=="real",na.rm=TRUE)/1000
                             doubletThresholding(x, local=FALSE, dbr=dbr, ...)
                         })
            th.stats <- lapply(th, FUN=function(x) x$stats)
            th <- vapply(th, FUN=function(x) x$th, FUN.VALUE=numeric(1))
            d$class <- ifelse(d$score >= th[d$sample], "doublet", "singlet")
            if(verbose) message("Thresholds found:\n", 
                                paste(paste(names(th),round(th,3),sep="="),
                                      collapse=", "))
        }else{
            th <- doubletThresholding( d, local=scoreType=="xgb.local.optim", 
                                       dbr=dbr, ... )
            if(scoreType=="xgb.local.optim"){
                d$score.global <- d$score
                d$score <- th$finalScores
            }
            th.stats <- th$stats
            th <- th$th
            d$class <- ifelse(d$score >= th, "doublet", "singlet")
            if(verbose) message("Threshold found:", round(th,3))
        }
        metadata(d)$scDblFinder.stats <- th.stats
        metadata(d)$scDblFinder.threshold <- th
        d$nearestClass <- factor(d$nearestClass, levels = 0:1, 
                                 labels=c("cell","artificialDoublet"))
        dbr <- sum(d$class=="doublet" & d$src=="real")/sum(d$src=="real")
        if(verbose) message(sum(d$class=="doublet" & d$src=="real"), " (", 
                            round(100*dbr,1),"%) doublets called")
    }
    d
}

# add the relevant fields of the scDblFinder results table to the SCE
.scDblAddCD <- function(sce, d){
  d <- d[colnames(sce),]
  for(f in c("sample","cluster","distanceToNearest","nearestClass","difficulty",
             "ratio","cxds_score","weighted","score.global","score",
             "class","mostLikelyOrigin","originAmbiguous")){
    if(!is.null(d[[f]])) sce[[paste0("scDblFinder.",f)]] <- d[[f]]
  }
  if(!is.null(sce$scDblFinder.class)) sce$scDblFinder.class <- 
    relevel(as.factor(sce$scDblFinder.class),"singlet")
  if(is(d,"DataFrame") && !is.null(metadata(d)$scDblFinder.stats))
      metadata(sce)$scDblFinder.stats <- metadata(d)$scDblFinder.stats
  sce
}



.checkSCE <- function(sce){
  if(is(sce, "SummarizedExperiment")){
    sce <- as(sce, "SingleCellExperiment")
  }else if(!is(sce, "SingleCellExperiment")){
    if(is.null(dim(sce)) || any(sce<0))
      stop("`sce` should be a SingleCellExperiment, a SummarizedExperiment, ",
           "or an array (i.e. matrix, sparse matric, etc.) of counts.")
    message("Assuming the input to be a matrix of counts or expected counts.")
    sce <- SingleCellExperiment(list(counts=sce))
  }
  if( !("counts" %in% assayNames(sce)) ) 
      stop("`sce` should have an assay named 'counts'")
  counts(sce) <- as(counts(sce),"dgCMatrix")
  if(is.null(colnames(sce)))
      colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  if(is.null(row.names(sce)))
      row.names(sce) <- paste0("f",seq_len(nrow(sce)))
  sce
}