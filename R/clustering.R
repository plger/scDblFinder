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
#' @param returnType See return.
#'
#' @return By default, a vector of cluster labels. If 
#' `returnType='preclusters'`, returns the k-means pre-clusters. If
#' `returnType='metacells'`, returns the metacells aggretated by pre-clusters 
#' and the corresponding cell indexes. If `returnType='graph'`, returns the
#' graph of (meta-)cells and the corresponding cell indexes.
#' 
#' @importFrom igraph cluster_louvain membership
#' @importFrom intrinsicDimension maxLikGlobalDimEst
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
fastcluster <- function( x, k=NULL, rdname="PCA", nstart=3, iter.max=20, 
                         ndims=NULL, nfeatures=1000, 
                         returnType=c("clusters","preclusters","metacells",
                                      "graph"),
                         BPPARAM=SerialParam() ){
  returnType <- match.arg(returnType)
  x <- .getDR(x, ndims=ndims, nfeatures=nfeatures, rdname=rdname)
  if(is.null(k)) k <- min(2500, floor(nrow(x)/10))
  if((returnType != "clusters" || nrow(x)>1000) && nrow(x)>k){
    k <- kmeans(x, k, iter.max=iter.max, nstart=nstart)$cluster
    if(returnType=="preclusters") return(k)
    x <- t(vapply(split(seq_along(k),k), FUN.VALUE=numeric(ncol(x)), 
               FUN=function(i) colMeans(x[i,,drop=FALSE])))
    if(returnType=="metacells") return(list(meta=x,idx=k))
  }else{
    k <- seq_len(nrow(x))
  }
  x <- makeKNNGraph(x, k=min(max(2,floor(sqrt(length(unique(k))))-1),10),
                    BPPARAM=BPPARAM)
  if(returnType=="graph") return(list(k=k, graph=x))
  cl <- membership(cluster_louvain(x))
  cl[k]
}

#' @importFrom scater runPCA 
#' @importFrom scuttle logNormCounts librarySizeFactors computeLibraryFactors 
#' @importFrom BiocSingular IrlbaParam
#' @import SingleCellExperiment
.prepSCE <- function(sce, ndims=30, nfeatures=1000){
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
                      ntop=min(nfeatures,nrow(sce)), BSPARAM=IrlbaParam())
    }
    sce
}

.getDR <- function(x, ndims=30, nfeatures=1000, rdname="PCA"){
  if(!(rdname %in% reducedDimNames(x)))
    x <- .prepSCE(x, ndims=ndims, nfeatures=nfeatures)
  x <- reducedDim(x, rdname)
  if(is.null(ndims)){
    ndims <- maxLikGlobalDimEst(x,k=20)$dim.est
    ndims <- min(c(50,ceiling(ndims),ncol(x)),na.rm=TRUE)
  }
  x[,seq_len(min(ncol(x),as.integer(ndims)))]
}

.getMetaGraph <- function(x, clusters, BPPARAM=SerialParam()){
  x <- t(vapply(split(seq_along(clusters),clusters), FUN.VALUE=numeric(ncol(x)), 
               FUN=function(i) colMeans(x[i,,drop=FALSE])))
  makeKNNGraph(x, k=min(max(2,floor(sqrt(length(unique(clusters))))-1),10),
               BPPARAM=BPPARAM)
}