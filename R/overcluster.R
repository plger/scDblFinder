#' overcluster
#'
#' This function deliberately overclusters based on the desired range of cluster
#' size. It uses iterative calls to `igraph::cluster_fast_greedy` until no 
#' cluster is above the size limits, and merges clusters that are too small.
#'
#' @param x A numeric matrix, with entities (e.g. cells) as columns and features
#'  (e.g. genes) as rows. Alternatively, an object of class `igraph`.
#' @param min.size The minimum cluster size (applies after splitting, and hence 
#' overrides `max.size`)
#' @param max.size The maximum cluster size. If omitted, will be calculated on 
#' the basis of the population size and initial number of clusters.
#' @param rdname Optional name of the reduced dimension to use.
#'
#' @return A vector of cluster labels.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' cc <- suppressWarnings(overcluster(m,min.size=5))
#' table(cc)
#' 
#' @importFrom scran buildSNNGraph
#' @importFrom igraph V cluster_fast_greedy membership
#' @export
overcluster <- function( x, min.size=50, max.size=NULL, rdname="PCA", ...){
  if(is.igraph(x)){
      N <- length(V(x))
      g <- x
  }else{
    x <- .prepSCE(x, ...)
    g <- buildSNNGraph(x, use.dimred=rdname, BNPARAM=AnnoyParam())
  }
  cl <- membership(cluster_fast_greedy(g, merges=FALSE, modularity=FALSE))
  resplitClusters(g, cl=cl, min.size=min.size, max.size=max.size)
}

.getMaxSize <- function(cl, max.size=NULL, min.size){
  if(!is.null(max.size)) return(max.size)
  ceiling(max(length(cl)/(2*length(unique(cl))),min.size+1))
}

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
#' @param nfeatures Number of features to use if doing PCA
#'
#' @return A vector of cluster labels
#' @export
fastcluster <- function( x, k=NULL, rdname="PCA", nstart=2, iter.max=20, 
                         ndims=30, nfeatures=1000){
  if(!(rdname %in% reducedDimNames(x)))
    x <- .prepSCE(x, ndims=ndims, nfeatures=nfeatures)
  x <- reducedDim(x, rdname)
  if(is.null(k)) k <- min(3000, floor(nrow(x)/10))
  k <- kmeans(x, k, iter.max=iter.max, nstart=nstart)$cluster
  ag <- sapply(split(names(k),k), FUN=function(i) colMeans(x[i,,drop=FALSE]))
  cl <- membership(cluster_louvain(buildKNNGraph(ag)))
  cl <- cl[k]
  cl
}


#' resplitClusters
#' 
#' Split (re-cluster) clusters of an existing graph-based clustering that are 
#' above a certain size
#'
#' @param g An object of class `igraph`
#' @param cl A vector of cluster labels corresponding to the nodes of `g`. If 
#' ommited, a new clustering will be run using `igraph::cluster_fast_greedy`.
#' @param max.size The maximum cluster size
#' @param min.size The minimum cluster size (default none). If given, this 
#' overrides `max.size`.
#' @param renameClusters Logical; whether to rename clusters
#' @param iterative Logical; whether to resplit until no cluster is above the 
#' size limit or no improvement is made (default TRUE). If FALSE, splits each 
#' cluster once.
#' @param nodesizes Optional size of each node of the graph
#'
#' @return A vector of cluster assignments.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' g <- scran::buildSNNGraph(rankTrans(m))
#' table(resplitClusters(g, min.size=2, max.size=20))
#' 
#' @importFrom igraph V E cluster_fast_greedy membership subgraph modularity
#' @export
resplitClusters <- function( g, cl=NULL, max.size=500, min.size=50, 
                             renameClusters=TRUE, iterative=TRUE, 
                             nodesizes=NULL ){
    if(!is.null(max.size) && !is.null(min.size) && max.size<min.size) 
        stop("max.size and min.size are incompatible")
    if(is.null(cl)){
        ## no initial clustering provided - run one
        cl <- membership(cluster_fast_greedy(g))
    }
    ll1 <- split(seq_len(length(V(g))), cl) # split nodes by cluster
    # restrict to clusters >limit
    ll1 <- ll1[which(vapply(ll1,FUN.VALUE=integer(1),FUN=length)>max.size)]
    # run clustering of clusters above size limit:
    ll2 <- lapply(ll1, FUN=function(x){
        membership(cluster_fast_greedy(
            suppressWarnings(subgraph(g, x))
        ))
    })
    ## update global cluster labels
    for(i in names(ll2)){
        if(length(unique(ll2[[i]]))>1){
            cl[ll1[[i]]] <- max(cl)+ll2[[i]]
        }
    }
    ## repeat until no more improvement or no cluster is above limit
    while(iterative && max(.getClusterSizes(cl, nodesizes))>max.size){
        newcl <- resplitClusters( g, cl, max.size=max.size, min.size=NULL, 
                                  renameClusters=FALSE, iterative=FALSE,
                                  nodesizes=nodesizes )
        if(identical(cl, newcl)) iterative <- FALSE
        cl <- newcl
    }
    if(!is.null(min.size)){
        ## merge clusters below minimum
        ## adapted from scran:::.merge_closest_graph()
        oldcl <- NULL
        while( min(.getClusterSizes(cl, nodesizes))<min.size && 
               !identical(oldcl, cl) ){
            cs <- .getClusterSizes(cl, nodesizes)
            oldcl <- cl
            min.cs <- names(cs)[which.min(cs)]
            other.cl <- setdiff(names(cs), min.cs)
            mod <- vapply(other.cl, FUN.VALUE=double(1), FUN=function(x){
                cl2 <- cl
                cl2[which(cl2==min.cs)] <- x
                suppressWarnings(modularity(g, cl2, weights = E(g)$weight))
            })
            cl[which(cl==min.cs)] <- other.cl[which.max(mod)]
            if(identical(cl,oldcl)){
                warning("Clusters could not meet criteria.")
                break
            }
        }
    }
    if(renameClusters){
        cl <- as.integer(as.factor(cl))
    }
    cl
}

.getClusterSizes <- function(cl, nodesizes=NULL){
  if(is.null(nodesizes)) return(table(cl))
  if(length(cl)!=length(nodesizes) && 
     (is.null(names(nodesizes)) || is.null(names(cl)) ||
     !all(names(cl) %in% names(nodesizes)) ) )
     stop("Cannot match nodes sizes and clusters.")
  if(!is.null(names(cl)) && is.null(names(cl)))
    nodesizes <- nodesizes[names(cl)]
  sapply(split(nodesizes, cl), FUN=sum)
}


#' @importFrom scater runPCA logNormCounts librarySizeFactors computeLibraryFactors 
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
        sce <- runPCA(sce, ncomponents=ndims, ntop=min(nfeatures,nrow(sce)),
                      BSPARAM=IrlbaParam())
    }
    sce
}
