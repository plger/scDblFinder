#' overcluster
#'
#' This function deliberately overclusters based on the desired range of cluster size.
#' It first calculates a SNN network viar `scran::buildSNNGraph`, then runs 
#' `igraph::cluster_fast_greedy` until no cluster is above the size limits, and merges 
#' clusters that are too small.
#'
#' @param e A numeric matrix, with entities (e.g. cells) as columns and features (e.g. 
#' genes) as rows.
#' @param rtrans Transformation to apply, either 'rankTrans' (default, dense step-
#' preserving rank transformation, see `rankTrans`), 'scran' (default; see 
#' `scran::scaledColRanks`), or 'none' (data taken as-is)
#' @param min.size The minimum cluster size (applies after splitting, and hence overrides
#' `max.size`)
#' @param max.size The maximum cluster size
#'
#' @return A vector of cluster labels.
#' @export
overcluster <- function(e, rtrans=c("scran","rankTrans","none"), min.size=20, max.size=300){
  e <- switch( match.arg(rtrans),
               scran=scran::scaledColRanks(e),
               rankTrans=rankTrans(e),
               e)
  g <- scran::buildSNNGraph(e)
  resplitClusters(g, min.size=min.size, max.size=max.size)
}


#' resplitClusters
#' 
#' Split (re-cluster) clusters of an existing graph-based clustering that are above a certain size
#'
#' @param g An object of class `igraph`
#' @param cl A vector of cluster labels corresponding to the nodes of `g`. If ommited,
#' a new clustering will be run using `igraph::cluster_fast_greedy`.
#' @param max.size The maximum cluster size
#' @param min.size The minimum cluster size (default none). If given, this overrides 
#' `max.size`.
#' @param renameClusters Logical; whether to rename clusters
#' @param iterative Logical; whether to resplit until no cluster is above the size limit 
#' or no improvement is made (default TRUE). If FALSE, splits each cluster once.
#'
#' @return A vector of cluster assignments.
#' @export
resplitClusters <- function(g, cl=NULL, max.size=500, min.size=50, renameClusters=TRUE, iterative=TRUE){
    library(igraph)
    if(is.null(cl)){
      # no initial clustering provided - run one
    	cl <- membership(cluster_fast_greedy(g))
    }
    ll1 <- split(1:length(V(g)), cl) # split nodes by cluster
    ll1 <- ll1[which(sapply(ll1,length)>max.size)] # restrict to cluster above size limit
    # run clustering of clusters above size limit:
    ll2 <- lapply(ll1, FUN=function(x){
    	membership(cluster_fast_greedy(suppressWarnings(subgraph(g, x))))
    })
    # update global cluster labels
    for(i in names(ll2)){
    	if(length(unique(ll2[[i]]))>1){
    		cl[ll1[[i]]] <- max(cl)+ll2[[i]]
    	}
    }
    # repeat until no more improvement or no cluster is above limit
    while(iterative && max(table(cl))>max.size){
    	newcl <- resplitClusters(g, cl, max.size=max.size, min.size=NULL, renameClusters=FALSE, iterative=FALSE)
    	if(identical(cl, newcl)) iterative <- FALSE
    	cl <- newcl
    }
    if(!is.null(min.size)){
      # merge clusters below minimum size
    	cl <- scran:::.merge_closest_graph(g, cl, min.size=min.size)
    }	
    if(renameClusters){
    	cl <- as.integer(as.factor(cl))
    }
    cl
}

