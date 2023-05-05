#' fastcluster
#'
#' Performs a fast two-step clustering: first clusters using k-means with a very
#' large k, then uses louvain clustering of the k cluster averages and reports
#' back the cluster labels.
#'
#' @param x An object of class SCE
#' @param k The number of k-means clusters to use in the primary step (should
#' be much higher than the number of expected clusters). Defaults to 1/10th of
#' the number of cells with a maximum of 3000.
#' @param rdname The name of the dimensionality reduction to use.
#' @param nstart Number of starts for k-means clustering
#' @param iter.max Number of iterations for k-means clustering
#' @param ndims Number of dimensions to use
#' @param nfeatures Number of features to use (ignored if `rdname` is given and
#' the corresponding dimensional reduction exists in `sce`)
#' @param verbose Logical; whether to output progress messages
#' @param returnType See return.
#' @param ... Arguments passed to `scater::runPCA` (e.g. BPPARAM or BSPARAM) if
#' `x` does not have `rdname`.
#'
#' @return By default, a vector of cluster labels. If
#' `returnType='preclusters'`, returns the k-means pre-clusters. If
#' `returnType='metacells'`, returns the metacells aggretated by pre-clusters
#' and the corresponding cell indexes. If `returnType='graph'`, returns the
#' graph of (meta-)cells and the corresponding cell indexes.
#'
#' @importFrom igraph cluster_louvain membership
#' @importFrom scran buildKNNGraph
#' @importFrom stats kmeans
#'
#' @examples
#' sce <- mockDoubletSCE()
#' sce$cluster <- fastcluster(sce)
#'
#' @export
#' @importFrom bluster makeKNNGraph
#' @importFrom igraph membership cluster_louvain
#' @importFrom DelayedArray rowsum
fastcluster <- function( x, k=NULL, rdname="PCA", nstart=3, iter.max=50,
                         ndims=NULL, nfeatures=1000, verbose=TRUE,
                         returnType=c("clusters","preclusters","metacells",
                                      "graph"), ...){
  returnType <- match.arg(returnType)
  x <- .getDR(x, ndims=ndims, nfeatures=nfeatures, rdname=rdname,
              verbose=verbose, ...)
  if(is.null(k)) k <- min(2500, floor(nrow(x)/10))
  if((returnType != "clusters" || nrow(x)>1000) && nrow(x)>k){
    if(verbose) message("Building meta-cells")
    k <- kmeans(x, k, iter.max=iter.max, nstart=nstart)$cluster
    if(returnType=="preclusters") return(k)
    x <- rowsum(x, k)
    x <- x/as.integer(table(k)[rownames(x)])
    if(returnType=="metacells") return(list(meta=x,idx=k))
  }else{
    k <- seq_len(nrow(x))
  }
  if(verbose) message("Building KNN graph and clustering")
  x <- makeKNNGraph(x, k=min(max(2,floor(sqrt(length(unique(k))))-1),10))
  if(returnType=="graph") return(list(k=k, graph=x))
  cl <- membership(cluster_louvain(x))
  cl[k]
}

#' @importFrom scater runPCA
#' @importFrom scuttle logNormCounts librarySizeFactors computeLibraryFactors
#' @importFrom BiocSingular IrlbaParam
#' @import SingleCellExperiment
.prepSCE <- function(sce, ndims=30, nfeatures=1000, ...){
    if(!("logcounts" %in% assayNames(sce))){
        if(is.null(librarySizeFactors(sce)))
            sce <- computeLibraryFactors(sce)
        ls <- librarySizeFactors(sce)
        if(any(is.na(ls) | ls==0))
            stop("Some of the size factors are invalid. Consider removing",
                 "cells with sizeFactors of zero, or filling in the",
                 "`logcounts' assay yourself.")
        sce <- logNormCounts(sce)
    }
    if(!("PCA" %in% reducedDimNames(sce))){
        sce <- runPCA(sce, ncomponents=ifelse(is.null(ndims),30,ndims),
                      ntop=min(nfeatures,nrow(sce)),
                      BSPARAM=IrlbaParam(), ...)
    }
    sce
}

.getDR <- function(x, ndims=30, nfeatures=1000, rdname="PCA", verbose=TRUE, ...){
  if(!(rdname %in% reducedDimNames(x))){
    if(verbose) message("Reduced dimension not found - running PCA...")
    x <- .prepSCE(x, ndims=ndims, nfeatures=nfeatures, ...)
  }
  x <- reducedDim(x, rdname)
  if(is.null(ndims)) dims <- 20
  x[,seq_len(min(ncol(x),as.integer(ndims)))]
}

.getMetaGraph <- function(x, clusters, BPPARAM=SerialParam()){
  x <- rowsum(x, clusters)
  x <- x/as.integer(table(clusters)[rownames(x)])
  makeKNNGraph(x, k=min(max(2,floor(sqrt(length(unique(clusters))))-1),10),
               BPPARAM=BPPARAM)
}
