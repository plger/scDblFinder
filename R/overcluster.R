#' overcluster
#'
#' This function deliberately overclusters based on the desired range of cluster
#' size. It first calculates a SNN network viar `scran::buildSNNGraph`, then 
#' runs `igraph::cluster_fast_greedy` until no cluster is above the size limits,
#'  and merges clusters that are too small. By default, `rankTrans` is used on 
#'  the counts before, because it tends to produce over-clustering influenced by 
#' library size, which is desirable for producing artificial doublets.
#'
#' @param x A numeric matrix, with entities (e.g. cells) as columns and features
#'  (e.g. genes) as rows. Alternatively, an object of class `igraph`.
#' @param rtrans Transformation to apply, either 'rankTrans' (default, dense 
#' step-preserving rank transformation, see `rankTrans`), 'scran' (default; see 
#' `scran::scaledColRanks`), or 'none' (data taken as-is). Ignored if `x` is an
#' `igraph`.
#' @param min.size The minimum cluster size (applies after splitting, and hence 
#' overrides `max.size`)
#' @param max.size The maximum cluster size. If omitted, will be calculated on 
#' the basis of the population size and initial number of clusters.
#'
#' @return A vector of cluster labels.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' cc <- suppressWarnings(overcluster(m,min.size=5))
#' table(cc)
#' 
#' @importFrom scran buildSNNGraph scaledColRanks
#' @import igraph
#' @export
overcluster <- function( x, rtrans=c("rankTrans","scran","none"), min.size=50, 
                         max.size=NULL){
  if(is.igraph(x)){
      N <- length(V(x))
      g <- x
  }else{
      x <- switch( match.arg(rtrans),
                   scran=scaledColRanks(x),
                   rankTrans=rankTrans(x),
                   x)
      g <- buildSNNGraph(x)
      N <- ncol(e)
  }
  cl <- membership(cluster_fast_greedy(g))
  if(is.null(max.size))
      max.size <- ceiling(max(N/(2*length(unique(cl))),min.size+1))
  resplitClusters(g, cl=cl, min.size=min.size, max.size=max.size)
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
#'
#' @return A vector of cluster assignments.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' g <- scran::buildSNNGraph(rankTrans(m))
#' table(resplitClusters(g, min.size=2, max.size=20))
#' 
#' @import igraph
#' @export
resplitClusters <- function( g, cl=NULL, max.size=500, min.size=50, 
                             renameClusters=TRUE, iterative=TRUE ){
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
    while(iterative && max(table(cl))>max.size){
        newcl <- resplitClusters( g, cl, max.size=max.size, min.size=NULL, 
                                  renameClusters=FALSE, iterative=FALSE )
        if(identical(cl, newcl)) iterative <- FALSE
        cl <- newcl
    }
    if(!is.null(min.size)){
        ## merge clusters below minimum
        ## adapted from scran:::.merge_closest_graph()
        oldcl <- NULL
        while(min(table(cl))<min.size && !identical(oldcl, cl)){
            cs <- table(cl)
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


#' fastClust
#'
#' @param sce An object of class `SingleCellExperiment`
#' @param nfeatures For the PCA
#' @param k number of nearest neighbors
#' @param dims number of PCA dimensions
#' @param graph.type either snn or knn
#' @param method Either 'louvain' or 'fast_greedy'
#' @param BPPARAM Passed to scran for KNN/SNN graph generation
#' @param ... passed to `overcluster`
#'
#' @return The `SingleCellExperiment` object with an additional colData column
#' scDblFinder.clusters
#' 
#' @importFrom scran buildSNNGraph buildKNNGraph
#' @importFrom igraph cluster_fast_greedy cluster_louvain
#' @export
fastClust <- function( sce, nfeatures=1000, k=10, dims=20, 
                       graph.type=c("snn","knn"),
                       method=c("louvain","fast_greedy","overcluster"),
                       BPPARAM=BiocParallel::SerialParam(), ...){
    method <- match.arg(method)
    if(graph.type=="snn"){ 
        gfn <- buildSNNGraph
    }else{
        gfn <- buildKNNGraph
    }
    sce <- .prepSCE(sce, nfeatures=nfeatures, ndims=dims) 
    g <- gfn(sce, BPPARAM=BPPARAM, use.dimred="PCA", k=k)
    sce$scDblFinder.clusters <- switch(method,
        louvain=igraph::cluster_louvain(g)$membership,
        fast_greedy=igraph::cluster_fast_greedy(g)$membership,
        overcluster=overcluster(g, ...)
    )
    sce
}

.prepSCE <- function(sce, ndims=30, nfeatures=1000){
    if(!("logcounts" %in% assayNames(sce))) sce <- scater::normalize(sce)
    if(!("PCA" %in% reducedDimNames(sce))){
        sce <- runPCA(sce, ncomponents=ndims, ntop=max(nfeatures,nrow(sce)))
    }
    sce
}
