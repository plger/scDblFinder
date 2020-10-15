#' getArtificialDoublets
#'
#' Create expression profiles of random artificial doublets.
#'
#' @param x A count matrix, with features as rows and cells as columns.
#' @param n The approximate number of doublet to generate (default 3000).
#' @param clusters The optional clusters labels to use to build cross-cluster 
#' doublets.
#' @param adjustSize Logical; whether to adjust the size of the doublets using
#' the ratio between each cluster's median library size. Alternatively, a number
#' between 0 and 1 can be given, determining the proportion of the doublets for
#' which to perform the size adjustment.
#' @param propRandom The proportion of the created doublets that are fully
#'  random (default 0.1); the rest will be doublets created across clusters. 
#'  Ignored if `clusters` is NULL.
#' @param n.meta.cells The number of meta-cell per cluster to create. If given,
#' additional doublets will be created from cluster meta-cells.
#' @param meta.triplets Logical; whether to create triplets from meta cells. 
#' Ignored if `clusters` is missing.
#' @param traj.weightFun The function by which to weigh distances (steps) in 
#' trajectory-based doublet generation.
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
getArtificialDoublets <- function( x, n=3000, clusters=NULL, adjustSize=FALSE,
                                   propRandom=0.1, n.meta.cells=2,
                                   meta.triplets=TRUE, traj.weightFun=NULL ){
  ls <- Matrix::colSums(x)
  w <- which(ls>0 & ls>=quantile(ls,0.01) & ls<=quantile(ls,0.99))
  x <- x[,w,drop=FALSE]
  if(!is.null(clusters)){
    if(is.list(clusters)){
      clusters$k <- clusters$k[w]
      clo <- clusters
      clusters <- clusters$k
    }else{
      clo <- clusters <- clusters[w] 
    }
    clusters <- as.factor(clusters)
  }
  
  if(is.null(clusters) || propRandom==1){
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
    if(is.null(clusters)){
      ad.m <- createDoublets(x, ad, adjustSize=FALSE, prefix="rDbl.")
      oc <- NULL
    }else{
      ad.m <- createDoublets(x, ad, clusters=clusters, adjustSize=adjustSize, 
                             prefix="rDbl.")
      oc <- matrix(clusters[as.numeric(ad)],ncol=2)
      w <- which(oc[,1]>oc[,2])
      oc[w,1:2] <- oc[w,2:1]
      oc <- paste(oc[,1],oc[,2],sep="+")
    }
    row.names(ad.m) <- row.names(x)
    return(list( counts=ad.m, origins=as.factor(oc) ))
  }
  
  if((nr <- ceiling(n*propRandom))>0){
    ad.m <- getArtificialDoublets(x, n=nr, clusters=clusters, 
                                  propRandom=1, adjustSize=adjustSize)
    oc <- ad.m$origins
    ad.m <- ad.m$counts
    n <- ceiling(n*(1-propRandom))
  }else{
    ad.m <- x[,c(),drop=FALSE]
    oc <- character()
  }
  
  if(length(unique(clusters))<3) n.meta.cells <- 0
  
  # create doublets across clusters:
  n <- ceiling(n)
  ca <- getCellPairs(clo, n=ifelse(n.meta.cells>0,ceiling(n*0.9),n), 
                     weightFun=traj.weightFun)
  m2 <- createDoublets(x, ca, clusters=clusters, adjustSize=adjustSize)
  oc <- c(oc, as.character(ca$orig.clusters))
  names(oc) <- names(m2)
  ad.m <- cbind(ad.m, m2)
  rm(m2)
  gc(verbose=FALSE)
  
  if(n.meta.cells>0){
    # create doublets from meta cells:
    meta <- .getMetaCells(x, clusters, n.meta.cells=n.meta.cells, 
                          meta.cell.size=30)
    cl2 <- rep(unique(clusters),each=n.meta.cells)
    ca <- getCellPairs(cl2, n=ceiling(n*0.1))
    m2 <- createDoublets(meta, ca, clusters=cl2, adjustSize=adjustSize, 
                         prefix="artMetaDbl.")
    oc2 <- ca$orig.clusters
    names(oc2) <- colnames(m2)
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
    meta <- .getMetaCells(x[,w], clusters[w], n.meta.cells=1, 
                          meta.cell.size=100)
    i <- seq_len(ncol(meta))
    ca <- expand.grid(i, i, i)
    ca <- ca[ca[,1]<ca[,2] & ca[,2]<ca[,3],,drop=FALSE]
    m2 <- round((meta[,ca[,1],drop=FALSE]+meta[,ca[,2]]+meta[,ca[,3]])/2)
    if(is(ad.m,"dgCMatrix")) m2 <- as(m2,"dgCMatrix")
    oc2 <- rep(NA_character_, ncol(m2))
    names(oc2) <- colnames(m2) <- 
      paste0("artTriplet.", ncol(ad.m)+seq_len(ncol(m2)))
    ad.m <- cbind(ad.m, m2)
    oc <- c(oc, oc2)
  }
  
  list( counts=ad.m, origins=as.factor(oc) )
}

#' getCellPairs
#'
#' Given a vector of cluster labels, returns pairs of cross-cluster cells
#' 
#' @param x A vector of cluster labels for each cell, or a graph
#' @param n The number of cell pairs to obtain
#' @param ... Further arguments, for instance the `k` vector of precluster 
#' labels if `x` is a metacell graph.
#'
#' @return A data.frame with the columns
#' @export
#'
#' @examples
#' # create random labels
#' x <- sample(head(LETTERS), 100, replace=TRUE)
#' getCellPairs(x, n=6)
getCellPairs <- function(x, n=1000, ...){
  if(is(x,"igraph")) return(.getCellPairsFromGraph(x, n=n, ...))
  if(is.list(x) && all(c("graph","k") %in% names(x)))
     return(.getCellPairsFromGraph(x$graph, k=x$k, n=n, ...))
  .getCellPairsFromClusters(x, n=n, ...)
}

# get cross-cluster pairs of cells
.getCellPairsFromClusters <- function(clusters, n=1000, ls=NULL, q=c(0.2,0.9),
                                      ...){
  cli <- split(seq_along(clusters), clusters)
  if(!is.null(ls)){
    ls <- split(ls, clusters)
    for(f in names(cli)){
        qr <- as.numeric(quantile(ls[[f]], prob=q))
        cli[[f]] <- cli[[f]][which(ls[[f]]>=qr[1] & ls[[f]]<=qr[2])]
    }
  }
  ca <- expand.grid(seq_along(cli), seq_along(cli))
  ca <- ca[ca[,1]<ca[,2],]
  n <- ceiling(n/nrow(ca))
  oc <- paste( names(cli)[ca[,1]], names(cli)[ca[,2]], sep="+")
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca <- data.frame(ca, orig.clusters=rep(as.factor(oc), each=n))
  colnames(ca) <- c("cell1","cell2","orig.clusters")
  ca[!duplicated(ca),]
}

# get pairs of cells based on distances on the meta-cell graph
#' @importFrom igraph distances
.getCellPairsFromGraph <- function(g, k=seq_len(length(V(g))), n=1000, 
                                   weightFun=NULL, ...){
  cli <- split(seq_along(k), k)
  ca <- expand.grid(seq_along(cli), seq_along(cli))
  ca <- as.data.frame(ca[ca[,1]<ca[,2],])
  d <- distances(g)
  d[is.infinite(d)] <- max(d[!is.infinite(d)])
  if(is.null(weightFun)) weightFun <- function(d) sqrt(d)-0.5
  d <- weightFun(d)
  d[is.nan(d) | d<0] <- 0
  ca$dist <- as.integer(vapply(seq_len(nrow(ca)), FUN.VALUE=numeric(1), 
                    FUN=function(i) d[ca[i,1],ca[i,2]]))
  ca$n <- rpois(nrow(ca), n/sum(ca$dist)*ca$dist)
  ca <- ca[ca$n>0,]
  oc <- rep(factor(paste(names(cli)[ca[,1]], names(cli)[ca[,2]], sep="+")),ca$n)
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i){ 
    cbind( sample(cli[[ca[i,1]]],size=ca$n[i],replace=TRUE),
           sample(cli[[ca[i,2]]],size=ca$n[i],replace=TRUE) )
  }))
  ca <- data.frame(ca, orig.clusters=oc)
  colnames(ca) <- c("cell1","cell2","orig.clusters")
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


#' addDoublets
#'
#' Adds artificial doublets to an existing dataset
#'
#' @param x A count matrix of singlets, or a `SingleCellExperiment`
#' @param clusters A vector of cluster labels for each column of `x`
#' @param dbr The doublet rate
#' @param only.heterotypic Whether to add only heterotypic doublets.
#' @param adjustSize Whether to adjust the library sizes of the doublets.
#' @param prefix Prefix for the colnames generated.
#'
#' @return A `SingleCellExperiment`
#' @export
addDoublets <- function(x, clusters, dbr=(0.01*ncol(x)/1000), 
                        only.heterotypic=TRUE, adjustSize=FALSE, 
                        prefix="doublet."){
  if(is(x, "SingleCellExperiment")) x <- counts(x)
  ed <- round(getExpectedDoublets(clusters, dbr=dbr, 
                                  only.heterotypic=only.heterotypic))
  ed <- ed[ed>0]
  ca <- getCellPairs(clusters, n=length(ed)*(30+max(ed)), ls=Matrix::colSums(x))
  ca <- split(ca, ca$orig.clusters)
  ca <- do.call(rbind, lapply(names(ed), FUN=function(i){
    ca2 <- ca[[i]][,1:2]
    ca2[sample.int(nrow(ca2), ed[[i]]),,drop=FALSE]
  }))
  m2 <- createDoublets(x, ca, clusters, adjustSize=adjustSize, prefix=prefix)
  cd <- data.frame(type=rep(c("singlet","doublet"),c(ncol(x),ncol(m2))),
                   cluster=c(clusters,rep(NA, ncol(m2))))
  SingleCellExperiment(list(counts=cbind(x, m2)), colData=cd)
}

#' createDoublets
#' 
#' Creates artificial doublet cells by combining given pairs of cells
#' 
#' @param x A count matrix of real cells
#' @param dbl.idx A matrix or data.frame with pairs of cell indexes stored in 
#' the first two columns.
#' @param clusters An optional vector of cluster labels (for each column of `x`)
#' @param resamp Logical; whether to resample the doublets using the poisson
#' distribution.
#' @param adjustSize Logical; whether to adjust the size of the doublets using
#' the median sizes per cluster of the originating cells. Requires `clusters` to
#' be given. Alternatively to a logical value, a number between 0 and 1 can be 
#' given, determining the proportion of the doublets for which to perform the 
#' size adjustment.
#' @param prefix Prefix for the colnames generated.
#'
#' @return A matrix of artificial doublets.
#' @export
createDoublets <- function(x, dbl.idx, clusters=NULL, resamp=TRUE, 
                           adjustSize=FALSE, prefix="dbl."){
  dgc <- is(x,"dgCMatrix")
  adjustSize <- as.numeric(adjustSize)
  if(adjustSize>1 || adjustSize<0)
      stop("`adjustSize` should be a logical or a number between 0 and 1.")
  wAd <- sample.int(nrow(dbl.idx), size=round(adjustSize*nrow(dbl.idx)))
  wNad <- setdiff(seq_len(nrow(dbl.idx)),wAd)
  #message("wAd=",length(wAd)," wNad=",length(wNad))
  x1 <- x[,dbl.idx[wNad,1]]+x[,dbl.idx[wNad,2]]
  if(length(wAd)>0){
    if(is.null(clusters)) stop("If `adjustSize=TRUE`, clusters must be given.")
    dbl.idx <- as.data.frame(dbl.idx[wAd,,drop=FALSE])
    ls <- Matrix::colSums(x)
    csz <- vapply(split(ls,clusters), FUN=median, FUN.VALUE=numeric(1))
    dbl.idx$ls.ratio <- ls[dbl.idx[,1]]/(ls[dbl.idx[,1]]+ls[dbl.idx[,2]])
    ls1 <- csz[as.character(clusters[dbl.idx[,1]])]
    ls2 <- csz[as.character(clusters[dbl.idx[,2]])]
    dbl.idx$factor <- (dbl.idx$ls.ratio+ls1/(ls1+ls2))/2
    dbl.idx$factor[dbl.idx$factor>0.8] <- 0.8
    dbl.idx$factor[dbl.idx$factor<0.2] <- 0.2
    dbl.idx$ls <- (ls[dbl.idx[,1]]+ls[dbl.idx[,2]])/2
    x2 <- x[,dbl.idx[,1]]*dbl.idx$factor+x[,dbl.idx[,2]]*(1-dbl.idx$factor)
    x2 <- tryCatch(x2 %*% diag(dbl.idx$ls/Matrix::colSums(x2)),
                   error=function(e) t(t(x2)/Matrix::colSums(x2)))
    x1 <- cbind(x1,x2)
    rm(x2)
  }
  x <- x1
  rm(x1)
  if(resamp){
    x <- matrix(as.integer(rpois(length(x), as.numeric(as.matrix(x)))), 
                nrow=nrow(x))
  }else{
    x <- round(x)
  }
  if(dgc) x <- as(x,"dgCMatrix")
  colnames(x) <- paste0( prefix, seq_len(ncol(x)) )
  x
}