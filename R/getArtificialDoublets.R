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
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(30,x) ) )
#' doublets <- getArtificialDoublets(m, 30)
#' 
#' @export
getArtificialDoublets <- function( x, n=3000, clusters=NULL,
                                   n.meta.cells=1, meta.triplets=TRUE ){
  if(is.null(clusters)){
    # create random combinations
    nr <- ifelse(is.null(clusters), n, ceiling(prop.fullyRandom*n))
    if(ncol(x)^2 <= nr){
      # few combinations, get them all
      ad <- expand.grid(seq_len(ncol(x)),seq_len(ncol(x)))
    }else{
      ad <- matrix(sample.int(ncol(x), 2*nr, replace=(2*nr >= ncol(x))),ncol=2)
    }
    # remove doublets of the same cell
    ad <- ad[which(ad[,1]!=ad[,2]),]
    # create doublets
    ad.m <- x[,ad[,1]]+x[,ad[,2]]
    colnames(ad.m) <- paste0("artDbl.", seq_len(ncol(ad.m)))
    return(ad.m)
  }
  
  if(length(unique(clusters))<3) n.meta.cells <- 0
    
  # create doublets across clusters:
  n <- ceiling(n)
  ca <- .getCellPairs(clusters, n=ifelse(n.meta.cells>0,ceiling(n*0.8),n))
  m2 <- x[,ca[,1]]+x[,ca[,2]]
  oc <- as.character(ca$orig.clusters)
  names(oc) <- colnames(m2) <- paste0( "artDbl.", seq_len(ncol(m2)) )
  ad.m <- m2
  rm(m2)
  gc(verbose=FALSE)
  
  if(n.meta.cells>0){
    # create doublets from meta cells:
    meta <- .getMetaCells(x, clusters, n.meta.cells=n.meta.cells, 
                          meta.cell.size=30)
    ca <- .getCellPairs(rep(unique(clusters),each=n.meta.cells), 
                        n=ceiling(n*0.2))
    m2 <- meta[,ca[,1]]+meta[,ca[,2]]
    if(is(ad.m,"dgCMatrix")) m2 <- as(round(m2),"dgCMatrix")
    oc2 <- ca$orig.clusters
    names(oc2) <- colnames(m2) <- 
      paste0("artMetaDbl.",ncol(ad.m)+seq_len(ncol(m2)))
    ad.m <- cbind(ad.m, m2)
    oc <- c(oc,as.character(oc2))
  }
  
  pc10 <- length(clusters)/10
  tt <- table(clusters)
  if(meta.triplets && length(tt)>2){
    # get clusters that have more than 10% of the cells
    cl2 <- names(tt)[tt>=pc10]
    # otherwise get the 3 largest clusters
    if(length(cl2)<3) cl2 <- names(sort(tt, decreasing=TRUE))[1:3]
    w <- which(clusters %in% cl2)
    # create triplets from meta cells:
    meta <- .getMetaCells(x[,w], clusters[w], n.meta.cells=1, meta.cell.size=100)
    i <- seq_len(ncol(meta))
    ca <- expand.grid(i, i, i)
    ca <- ca[ca[,1]<ca[,2] & ca[,2]<ca[,3],,drop=FALSE]
    m2 <- meta[,ca[,1],drop=FALSE]+meta[,ca[,2],drop=FALSE]+meta[,ca[,3],drop=FALSE]
    if(is(ad.m,"dgCMatrix")) m2 <- as(round(m2),"dgCMatrix")
    oc2 <- rep(NA_character_, ncol(m2))
    names(oc2) <- colnames(m2) <- 
      paste0("artTriplet.", ncol(ad.m)+seq_len(ncol(m2)))
    ad.m <- cbind(ad.m, m2)
    oc <- c(oc, oc2)
  }
  list( counts=ad.m, origins=as.factor(oc) )
}
