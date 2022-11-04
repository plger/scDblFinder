#' getArtificialDoublets
#'
#' Create expression profiles of random artificial doublets.
#'
#' @param x A count matrix, with features as rows and cells as columns.
#' @param n The approximate number of doublet to generate (default 3000).
#' @param clusters The optional clusters labels to use to build cross-cluster
#' doublets.
#' @param resamp Logical; whether to resample the doublets using the poisson
#' distribution. Alternatively, if a proportion between 0 and 1, the proportion
#' of doublets to resample.
#' @param halfSize Logical; whether to half the library size of doublets
#' (instead of just summing up the cells). Alternatively, a number between 0
#' and 1 can be given, determining the proportion of the doublets for which
#' to perform the size adjustment.
#' @param adjustSize Logical; whether to adjust the size of the doublets using
#' the ratio between each cluster's median library size. Alternatively, a number
#' between 0 and 1 can be given, determining the proportion of the doublets for
#' which to perform the size adjustment.
#' @param propRandom The proportion of the created doublets that are fully
#'  random (default 0.1); the rest will be doublets created across clusters.
#'  Ignored if `clusters` is NULL.
#' @param selMode The cell pair selection mode for inter-cluster doublet
#' generation, either 'uniform' (same number of doublets for each combination),
#' 'proportional' (proportion expected from the clusters' prevalences), or
#' 'sqrt' (roughly the square root of the expected proportion).
#' @param n.meta.cells The number of meta-cell per cluster to create. If given,
#' additional doublets will be created from cluster meta-cells. Ignored if
#' `clusters` is missing.
#' @param meta.triplets Logical; whether to create triplets from meta cells.
#' Ignored if `clusters` is missing.
#' @param trim.q A vector of two values between 0 and 1
#'
#' @return A list with two elements: `counts` (the count matrix of
#' the artificial doublets) and `origins` the clusters from which each
#' artificial doublets originated (NULL if `clusters` is not given).
#'
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=50),
#'                FUN=function(x) rpois(30,x) ) )
#' doublets <- getArtificialDoublets(m, 30)
#'
#' @export
getArtificialDoublets <- function( x, n=3000, clusters=NULL, resamp=0.25,
                                   halfSize=0.25, adjustSize=0.25, propRandom=0.1,
                                   selMode=c("proportional","uniform","sqrt"),
                                   n.meta.cells=2, meta.triplets=TRUE,
                                   trim.q=c(0.05,0.95) ){
  selMode <- match.arg(selMode)
  ls <- Matrix::colSums(x)
  stopifnot(is.numeric(trim.q), length(trim.q)==2)
  if(is.null(clusters)){
    w <- which(ls>0 & ls>=quantile(ls,min(trim.q)) & ls<=quantile(ls,max(trim.q)))
  }else{
    clo <- clusters
    if(is.list(clusters)) clusters <- clusters$k
    q <- sapply(split(ls, clusters), FUN=function(x){
      if(length(x)<10) return(c(0,max(x)))
      quantile(x, prob=sort(trim.q), na.rm=TRUE)
    })
    w <- which(ls>0 & ls>=q[1,clusters] & ls<=q[2,clusters] & !is.na(clusters))
    clusters <- droplevels(as.factor(clusters[w]))
    if(is.list(clo)){
      clo$k <- clo$k[w]
    }else{
      clo <- clusters
    }
  }
  x <- x[,w,drop=FALSE]

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
      ad.m <- createDoublets(x, ad, adjustSize=FALSE, resamp=resamp,
                             halfSize=halfSize, prefix="rDbl.")
      oc <- rep(NA,ncol(ad.m))
    }else{
      oc1 <- matrix(clusters[as.numeric(ad)],ncol=2)
      ad.m <- createDoublets(x, ad, clusters=clusters, adjustSize=adjustSize,
                             halfSize=halfSize, resamp=resamp, prefix="rDbl.")
      oc <- paste(oc1[,2],oc1[,1],sep="+")
      oc[which(oc1[,2]==oc1[,1])] <- NA
    }
    row.names(ad.m) <- row.names(x)
    return(list( counts=ad.m, origins=as.factor(oc) ))
  }

  if((nr <- ceiling(n*propRandom))>0){
    ad.m <- getArtificialDoublets(x, n=nr, clusters=clusters,
                                  halfSize=halfSize, resamp=resamp,
                                  propRandom=1, adjustSize=adjustSize)
    oc <- ad.m$origins
    ad.m <- ad.m$counts
    n <- ceiling(n*(1-propRandom))
  }else{
    ad.m <- as(as.matrix(x[,1:2])[,c(),drop=FALSE], "CsparseMatrix")
    oc <- character()
  }

  if(length(unique(clusters))<3) n.meta.cells <- 0

  # create doublets across clusters:
  n <- ceiling(n)
  ca <- getCellPairs(clo, n=ifelse(n.meta.cells>0,ceiling(n*0.9),n),
                     selMode=selMode)
  m2 <- createDoublets(x, ca, clusters=clusters, adjustSize=adjustSize,
                       halfSize=halfSize, resamp=resamp)
  oc <- c(oc, as.character(ca$orig.clusters))
  names(oc) <- names(m2)
  ad.m <- cbind(ad.m, m2)
  rm(m2)
  gc(verbose=FALSE)

  if(n.meta.cells>0){
    # create doublets from meta cells:
    meta <- .getMetaCells(x, clusters, n.meta.cells=n.meta.cells,
                          meta.cell.size=30)
    cl2 <- meta$clusters
    meta <- meta$mu
    ca <- getCellPairs(cl2, n=ceiling(n*0.1))
    m2 <- createDoublets(meta, ca, clusters=cl2, adjustSize=FALSE,
                         halfSize=FALSE, resamp=TRUE, prefix="artMetaDbl.")
    oc2 <- ca$orig.clusters
    names(oc2) <- colnames(m2)
    ad.m <- cbind(ad.m, m2)
    oc <- c(oc,as.character(oc2))
  }
  pc10 <- length(clusters)/10
  tt <- table(clusters)
  if(meta.triplets && sum(tt>=pc10)>2){
    # get clusters that have more than 10% of the cells
    cl2 <- names(tt)[tt>=pc10]
    # otherwise get the 3 largest clusters
    if(length(cl2)<3) cl2 <- names(sort(tt, decreasing=TRUE))[1:3]
    w <- which(clusters %in% cl2)
    # create triplets from meta cells:
    meta <- .getMetaCells(x[,w], clusters[w], n.meta.cells=1,
                          meta.cell.size=100)$mu
    i <- seq_len(ncol(meta))
    ca <- expand.grid(i, i, i)
    ca <- ca[ca[,1]<ca[,2] & ca[,2]<ca[,3],,drop=FALSE]
    m2 <- round((meta[,ca[,1],drop=FALSE]+meta[,ca[,2]]+meta[,ca[,3]])/2)
    if(is(ad.m,"sparseMatrix")) m2 <- as(m2,"CsparseMatrix")
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
#' @param clusters A vector of cluster labels for each cell, or a list containing
#' metacells and graph
#' @param n The number of cell pairs to obtain
#' @param ls Optional library sizes
#' @param q Library size quantiles between which to include cells (ignored if 
#' `ls` is NULL)
#' @param selMode How to decide the number of pairs of each kind to produce.
#' Either 'proportional' (default, proportional to the abundance of the 
#' underlying clusters), 'uniform' or 'sqrt'.
#' @param soft.min Minimum number of pairs of a given type.
#'
#' @return A data.frame with the columns
#' @export
#'
#' @examples
#' # create random labels
#' x <- sample(head(LETTERS), 100, replace=TRUE)
#' getCellPairs(x, n=6)
getCellPairs <- function(clusters, n=1000, ls=NULL, q=c(0.1,0.9),
                         selMode="proportional", soft.min=5){
  if(is.factor(clusters)) clusters <- droplevels(clusters)
  cli <- split(seq_along(clusters), clusters)
  if(!is.null(ls)){
    ls <- split(ls, clusters)
    for(f in names(cli)){
        qr <- as.numeric(quantile(ls[[f]], prob=q))
        cli[[f]] <- cli[[f]][which(ls[[f]]>=qr[1] & ls[[f]]<=qr[2])]
    }
  }
  ca <- expand.grid(seq_along(cli), seq_along(cli))
  ca <- as.data.frame(ca[ca[,1]<ca[,2],])
  ca$orig <- factor(paste( names(cli)[ca[,1]], names(cli)[ca[,2]], sep="+"))
  if(selMode=="uniform"){
    ca$n <- ceiling(n/nrow(ca))
  }else{
    ed <- getExpectedDoublets(clusters, only.heterotypic=TRUE)
    ed <- ed*n/sum(ed)
    if(selMode=="sqrt") ed <- sqrt(ed)
    ed <- soft.min + ed
    ed <- ceiling(ed*n/sum(ed))
    ca$n <- ed[as.character(ca$orig)]
  }
  ca <- ca[ca$n>0,]
  lvls <- levels(ca$orig)
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i){
    cbind( cell1=sample(cli[[ca[i,1]]], size=ca[i,4],
                        replace=ca[i,4]>length(cli[[ca[i,1]]])),
           cell2=sample(cli[[ca[i,2]]], size=ca[i,4],
                        replace=ca[i,4]>length(cli[[ca[i,2]]])),
           orig.clusters=rep(ca[i,3],ca[i,4]) )
  }))
  ca <- as.data.frame(ca)
  ca$orig.clusters <- factor(ca$orig.clusters, levels=seq_along(lvls), labels=lvls)
  ca[!duplicated(ca),]
}


# creates within-cluster meta-cells from a count matrix
.getMetaCells <- function(x, clusters, n.meta.cells=20, meta.cell.size=20){
  if(is.factor(clusters)) clusters <- droplevels(clusters)
  cli <- split(seq_along(clusters), clusters)
  meta <- lapply(cli, FUN=function(x){
    lapply(seq_len(n.meta.cells), FUN=function(y){
      sample(x,min(ceiling(0.6*length(x)),meta.cell.size),replace=FALSE)
    })
  })
  cl2 <- rep(names(cli), lengths(meta))
  meta <- vapply(unlist(meta, recursive=FALSE), FUN.VALUE=double(nrow(x)),
                 FUN=function(y){ Matrix::rowMeans(x[,y,drop=FALSE]) })
  colnames(meta) <- paste0("metacell.",seq_len(ncol(meta)))
  list(mu=meta, clusters=cl2)
}


#' addDoublets
#'
#' Adds artificial doublets to an existing dataset
#'
#' @param x A count matrix of singlets, or a
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param clusters A vector of cluster labels for each column of `x`
#' @param dbr The doublet rate
#' @param only.heterotypic Whether to add only heterotypic doublets.
#' @param adjustSize Whether to adjust the library sizes of the doublets.
#' @param prefix Prefix for the colnames generated.
#' @param ... Any further arguments to \code{\link{createDoublets}}.
#'
#' @return A `SingleCellExperiment` with the colData columns `cluster` and
#' `type` (indicating whether the cell is a singlet or doublet).
#'
#' @export
#' @examples
#' sce <- mockDoubletSCE(dbl.rate=0)
#' sce <- addDoublets(sce, clusters=sce$cluster)
addDoublets <- function(x, clusters, dbr=(0.01*ncol(x)/1000),
                        only.heterotypic=TRUE, adjustSize=FALSE,
                        prefix="doublet.", ...){
  ed <- round(getExpectedDoublets(clusters, dbr=dbr,
                                  only.heterotypic=only.heterotypic))
  ed <- ed[ed>0]
  if(length(ed)==0) return(x)
  if(is(x, "SingleCellExperiment")) x <- counts(x)
  ca <- getCellPairs(clusters, n=length(ed)*(30+max(ed)), ls=Matrix::colSums(x))
  ca <- split(ca, ca$orig.clusters)
  ca <- do.call(rbind, lapply(names(ed), FUN=function(i){
    ca2 <- ca[[i]]
    ca2[sample.int(nrow(ca2), ed[[i]]),,drop=FALSE]
  }))
  cl <- as.factor(c(as.character(clusters), as.character(ca[,3])))
  m2 <- createDoublets(x, ca, clusters, adjustSize=adjustSize, prefix=prefix)
  cd <- data.frame(type=rep(c("singlet","doublet"),c(ncol(x),ncol(m2))),
                   cluster=cl)
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
#' distribution. Alternatively, if a proportion between 0 and 1, the proportion
#' of doublets to resample.
#' @param halfSize Logical; whether to half the library size of doublets
#' (instead of just summing up the cells). Alternatively, a number between 0
#' and 1 can be given, determining the proportion of the doublets for which
#' to perform the size adjustment. Ignored if not resampling.
#' @param adjustSize Logical; whether to adjust the size of the doublets using
#' the median sizes per cluster of the originating cells. Requires `clusters` to
#' be given. Alternatively to a logical value, a number between 0 and 1 can be
#' given, determining the proportion of the doublets for which to perform the
#' size adjustment.
#' @param prefix Prefix for the colnames generated.
#'
#' @return A matrix of artificial doublets.
#' @export
#' @examples
#' sce <- mockDoubletSCE()
#' idx <- getCellPairs(sce$cluster, n=200)
#' art.dbls <- createDoublets(sce, idx)
createDoublets <- function(x, dbl.idx, clusters=NULL, resamp=0.5,
                           halfSize=0.5, adjustSize=FALSE, prefix="dbl."){
  adjustSize <- .checkPropArg(as.numeric(adjustSize),FALSE)
  halfSize <- .checkPropArg(as.numeric(halfSize),FALSE)
  resamp <- .checkPropArg(as.numeric(resamp),FALSE)
  if(is(x,"SingleCellExperiment")) x <- counts(x)
  if(adjustSize>1 || adjustSize<0)
      stop("`adjustSize` should be a logical or a number between 0 and 1.")
  if(halfSize>1 || halfSize<0)
      stop("`adjustSize` should be a logical or a number between 0 and 1.")
  wAd <- sample.int(nrow(dbl.idx), size=round(adjustSize*nrow(dbl.idx)))
  wNad <- setdiff(seq_len(nrow(dbl.idx)),wAd)
  x1 <- x[,dbl.idx[wNad,1],drop=FALSE]+x[,dbl.idx[wNad,2],drop=FALSE]
  if(length(wAd)>1){
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
    dbl.idx$ls <- (ls[dbl.idx[,1]]+ls[dbl.idx[,2]])
    x2 <- x[,dbl.idx[,1]]*dbl.idx$factor+x[,dbl.idx[,2]]*(1-dbl.idx$factor)
    x2 <- tryCatch(x2 %*% diag(dbl.idx$ls/Matrix::colSums(x2)),
                   error=function(e) t(t(x2)/Matrix::colSums(x2)))
    x1 <- cbind(x1,x2)
    rm(x2)
  }
  x <- x1
  rm(x1)
  if(halfSize>0){
    wAd <- sample.int(nrow(dbl.idx), size=ceiling(halfSize*nrow(dbl.idx)))
    if(length(wAd)>0)    x[,wAd] <- x[,wAd]/2
  }
  if(resamp>0){
    if(resamp!=halfSize) wAd <- sample.int(ncol(x), ceiling(resamp*ncol(x)))
    if(length(wAd)>0)
      x[,wAd] <- matrix(as.integer(rpois(nrow(x)*length(wAd),
                                       as.numeric(as.matrix(x[,wAd])))),
                         nrow=nrow(x))
  }else{
    x <- round(x)
  }
  x <- as(x,"CsparseMatrix")
  colnames(x) <- paste0( prefix, seq_len(ncol(x)) )
  x
}