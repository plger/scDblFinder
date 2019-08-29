#' getArtificialDoublets
#'
#' Create expression profiles of random artificial doublets.
#'
#' @param x A count matrix, with features as rows and cells as columns.
#' @param n The approximate number of doublet to generate (default 3000).
#' @param prop.fullyRandom The proportion of the created doublets that are fully
#'  random (default 0); the rest will be doublets created across clusters. 
#'  Ignored if `clusters` is NULL.
#' @param clusters The optional clusters labels to use to build cross-cluster 
#' doublets.
#' @param n.meta.cells The number of meta-cell per cluster to create. If given,
#' additional doublets will be created from cluster meta-cells.
#' @param meta.triplets Logical; whether to create triplets from meta cells. 
#' Ignored if `clusters` is missing.
#'
#' @return A count matrix for artificial doublets.
#' 
#' @import Matrix
#' @export
getArtificialDoublets <- function( x, n=3000, prop.fullyRandom=0, clusters=NULL,
                                   n.meta.cells=1, meta.triplets=TRUE ){
  if(prop.fullyRandom>0 | is.null(clusters)){
    # create random combinations
    nr <- ifelse(is.null(clusters), n, ceiling(prop.fullyRandom*n))
    if(ncol(x)^2 <= nr){
      # few combinations, get them all
      ad <- expand.grid(1:ncol(x),1:ncol(x))
    }else{
      ad <- matrix(sample.int(ncol(x), 2*nr, replace=(2*nr >= ncol(x))),ncol=2)
    }
    # remove doublets of the same cell
    ad <- ad[which(ad[,1]!=ad[,2]),]
    # create doublets
    ad.m <- x[,ad[,1]]+x[,ad[,2]]
    colnames(ad.m) <- paste0("arificialDoublet.", seq_len(ncol(ad.m)))
    
    if(is.null(clusters)) return(ad.m)
  }else{
    ad.m <- NULL
  }
  
  # create doublets across clusters:
  n <- ceiling(n*(1-prop.fullyRandom))
  ca <- .getCellPairs(clusters, n=ifelse(n.meta.cells>0,ceiling(n*0.8),n))
  m2 <- x[,ca[,1]]+x[,ca[,2]]
  colnames(m2) <- paste0( "arificialDoublet.",
                          ifelse(is.null(ad.m),0,ncol(ad.m))+1:ncol(m2) )
  if(is.null(ad.m)){
    ad.m <- m2
  }else{
    ad.m <- cbind(ad.m, m2)
  }
  
  if(n.meta.cells>0){
    # create doublets from meta cells:
    meta <- .getMetaCells(x, clusters, n.meta.cells=n.meta.cells, 
                          meta.cell.size=30)
    clusters <- rep(unique(clusters),each=n.meta.cells)
    ca <- .getCellPairs(clusters, n=ceiling(n*0.2))
    m2 <- meta[,ca[,1]]+meta[,ca[,2]]
    colnames(m2) <- paste0("arificialMetaDoublet.",ncol(ad.m)+1:ncol(m2))
    ad.m <- cbind(ad.m, m2)
  }
  
  if(meta.triplets){
    # create triplets from meta cells:
    if(length(unique(clusters))^3 > n){
      tt <- table(clusters)
      topc <- min(10,length(tt))
      warning("Too many clusters - will create triplets only for the ", topc,
              " largest clusters.")
      tt <- names(tt)[order(tt,decreasing=TRUE)[seq_len(topc)]]
      w <- which(clusters %in% tt)
      x <- x[,w]
      clusters <- clusters[w]
    }
    meta <- .getMetaCells(x, clusters, n.meta.cells=1, meta.cell.size=100)
    i <- seq_len(ncol(meta))
    ca <- expand.grid(i, i, i)
    ca <- ca[apply(ca,1,FUN=function(x){ length(unique(x)) })==3,]
    m2 <- meta[,ca[,1]]+meta[,ca[,2]]+meta[,ca[,3]]
    colnames(m2) <- paste0("arificialTriplet.", ncol(ad.m)+1:ncol(m2))
    ad.m <- cbind(ad.m, m2)
  }
  ad.m
}
