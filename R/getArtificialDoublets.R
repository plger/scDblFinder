#' getArtificialDoublets
#'
#' Create expression profiles of random artificial doublets.
#'
#' @param x A count matrix, with features as rows and cells as columns.
#' @param n The approximate number of doublet to generate (default 3000).
#' @param clusters The optional clusters labels to use to build cross-cluster 
#' doublets.
#' @param prop.fullyRandom The proportion of the created doublets that are fully
#'  random (default 0.1); the rest will be doublets created across clusters. 
#'  Ignored if `clusters` is NULL.
#' @param n.meta.cells The number of meta-cell per cluster to create. If given,
#' additional doublets will be created from cluster meta-cells.
#' @param meta.triplets Logical; whether to create triplets from meta cells. 
#' Ignored if `clusters` is missing.
#'
#' @return If `clusters` is Null, returns a count matrix of artificial doublets.
#' Otherwise, returns a list with two elements: `counts` (the count matrix of
#' the artificial doublets) and `origins` the clusters from which each 
#' artificial doublets originated.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(30,x) ) )
#' doublets <- getArtificialDoublets(m, 30)
#' 
#' @export
getArtificialDoublets <- function( x, n=3000, clusters=NULL, 
                                   prop.fullyRandom=0.1, n.meta.cells=1,
                                   meta.triplets=TRUE ){
  if(is.null(clusters)){
    # create random combinations
    if(ncol(x)^2 <= n){
      # few combinations, get them all
      ad <- expand.grid(seq_len(ncol(x)),seq_len(ncol(x)))
    }else{
      ad <- matrix(sample.int(ncol(x), 2*n, replace=(2*n >= ncol(x))),ncol=2)
    }
    # remove doublets of the same cell
    ad <- ad[which(ad[,1]!=ad[,2]),]
    # create doublets
    ad.m <- x[,ad[,1]]+x[,ad[,2]]
    colnames(ad.m) <- paste0("artDbl.", seq_len(ncol(ad.m)))
    return(ad.m)
  }
    
  if((nr <- ceiling(n*prop.fullyRandom))>0){
    ad.m <- getArtificialDoublets(x, n=nr, n.meta.cells=0)
    colnames(ad.m) <- paste0("R",colnames(ad.m))
    oc <- rep(NA_character_, ncol(ad.m))
    n <- ceiling(n*(1-prop.fullyRandom))
  }else{
    ad.m <- x[,c(),drop=FALSE]
    oc <- character()
  }
  
  if(length(unique(clusters))<3) n.meta.cells <- 0
  
  # create doublets across clusters:
  n <- ceiling(n)
  ca <- .getCellPairs(clusters, n=ifelse(n.meta.cells>0,ceiling(n*0.8),n))
  m2 <- x[,ca[,1]]+x[,ca[,2]]
  oc <- c(oc, as.character(ca$orig.clusters))
  names(oc) <- colnames(m2) <- paste0( "artDbl.", seq_len(ncol(m2)) )
  ad.m <- cbind(ad.m, m2)
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

# get random cross-cluster pairs of cells from a cluster assignment vector
.getCellPairs <- function(clusters, n=1000){
  cli <- split(seq_along(clusters), clusters)
  ca <- expand.grid(seq_along(cli), seq_along(cli))
  ca <- ca[ca[,1]<ca[,2],]
  n <- ceiling(n/nrow(ca))
  oc <- paste( names(cli)[ca[,1]], names(cli)[ca[,2]], sep="+")
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca <- data.frame(ca, orig.clusters=rep(as.factor(oc), each=n))
  ca[!duplicated(ca),]
}  

# creates within-cluster meta-cells from a count matrix
.getMetaCells <- function(x, clusters, n.meta.cells=20, meta.cell.size=20){
  if(is.factor(clusters)) clusters <- droplevels(clusters)
  cli <- split(seq_along(clusters), clusters)
  meta <- unlist(lapply(cli, FUN=function(x){
    lapply(seq_len(n.meta.cells), FUN=function(y){
      sample(x,min(ceiling(0.6*length(x)),meta.cell.size),replace=FALSE)
    })
  }), recursive=FALSE)
  meta <- vapply(meta, FUN.VALUE=double(nrow(x)), 
                 FUN=function(y){ Matrix::rowMeans(x[,y,drop=FALSE]) })
  colnames(meta) <- paste0("metacell.",seq_len(ncol(meta)))
  meta
}


