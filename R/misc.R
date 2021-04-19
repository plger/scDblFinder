.tablevec <- function(x){
  if(!is(x, "table")) x <- table(x)
  y <- as.numeric(x)
  names(y) <- names(x)
  y
}

.getMostLikelyOrigins <- function(knn, known.origins=NULL){
  origins <- t(vapply( seq_len(nrow(knn$orig)), FUN.VALUE=character(2),
                       FUN=function(i){
           if(all(is.na(knn$orig[i,]))) return(c(NA_character_,NA_character_))
           n <- .tablevec(knn$orig[i,])
           if(length(n)==1) return(c(names(n),FALSE))
           if(sum(n==max(n))==1) return(c(names(n)[which.max(n)],FALSE))
           x <- split(knn$distance[i,], knn$orig[i,])
           x <- sort(vapply(x,min,numeric(1)))
           if( (x[2]-x[1])/x[1] > 0.2 ) return(c(names(x)[1], FALSE))
           return(c(names(x)[1], TRUE))
                       }))
  origins <- as.data.frame(origins)
  colnames(origins) <- c("mostLikelyOrigin", "originAmbiguous")
  if(!is.null(known.origins)){
    w <- which(!is.na(known.origins))
    origins[w,1] <- as.character(known.origins[w])
    origins[w,2] <- "FALSE"
  }
  origins[,1] <- factor(origins[,1])
  origins[,2] <- as.logical(origins[,2])
  origins
}

#' getExpectedDoublets
#'
#' @param x A vector of cluster labels for each cell
#' @param dbr The expected doublet rate.
#' @param only.heterotypic Logical; whether to return expectations only for
#' heterotypic doublets
#'
#' @return The expected number of doublets of each combination of clusters
#'
#' @examples
#' # random cluster labels
#' cl <- sample(head(LETTERS,4), size=2000, prob=c(.4,.2,.2,.2), replace=TRUE)
#' getExpectedDoublets(cl)
#' @export
getExpectedDoublets <- function(x, dbr=NULL, only.heterotypic=TRUE){
  if(is(x,"SingleCellExperiment")){
    clusters <- x$scDblFinder.clusters
  }else if(is(x,"list")){
    clusters <- x$k
  }else{
    clusters <- x
  }
  clusters <- droplevels(as.factor(clusters))
  ncells <- length(clusters)
  if(is.null(dbr)) dbr <- (0.01*ncells/1000)

  cs <- table(clusters)/ncells
  eg <- expand.grid(seq_along(cs), seq_along(cs))
  eg <- eg[eg[,1]<=eg[,2],]
  expected <- apply(eg,1, FUN=function(x){
    cs[[x[[1]]]]*cs[[x[[2]]]]
  })
  expected <- dbr*ncells*expected/sum(expected)
  names(expected) <- apply(eg,1,FUN=function(x){
    paste(names(cs)[x], collapse="+")
  })
  if(only.heterotypic) expected <- expected[which(eg[,1]!=eg[,2])]
  expected
}

.castorigins <- function(e, val=NULL){
  if(is.table(e) || is.null(dim(e))){
    e <- cbind(do.call(rbind, strsplit(names(e),"+",fixed=TRUE)),
               data.frame(val=as.numeric(e)))
  }
  if(is.null(val)) val <- rev(colnames(e))[1]
  names(n) <- n <- unique(as.character(as.matrix(e[,1:2])))
  vapply(n, FUN.VALUE=numeric(length(n)), FUN=function(x){
    vapply(n, FUN.VALUE=numeric(1), FUN=function(y){
      if(x==y) return(NA)
      w <- which(e[,1] %in% c(x,y) & e[,2] %in% c(x,y))
      if(length(w)==0) return(0)
      sum(e[w,val])
    })
  })
}

#' selFeatures
#'
#' Selects features based on cluster-wise expression or marker detection, or a
#' combination.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}},
#' \code{\link[SingleCellExperiment]{SingleCellExperiment-class}} with a
#' 'counts' assay.
#' @param clusters Optional cluster assignments. Should either be a vector of
#' labels for each cell.
#' @param nfeatures The number of features to select.
#' @param propMarkers The proportion of features to select from markers (rather
#' than on the basis of high expression). Ignored if `clusters` isn't given.
#' @param FDR.max The maximum marker binom FDR to be included in the selection.
#' (see \code{\link[scran]{findMarkers}}).
#'
#' @return A vector of feature (i.e. row) names.
#' @export
#'
#' @importFrom scuttle sumCountsAcrossCells
#' @importFrom scran findMarkers
#' @examples
#' sce <- mockDoubletSCE()
#' selFeatures(sce, clusters=sce$cluster, nfeatures=5)
selFeatures <- function(sce, clusters=NULL, nfeatures=1000, propMarkers=0, FDR.max=0.05){
  if(nrow(sce)<=nfeatures) return(row.names(sce))
  if(is.null(clusters)) propMarkers <- 0
  g <- c()
  if((ng <- ceiling((1-propMarkers)*nfeatures))>0){
    if(is.null(clusters)){
      g <- row.names(sce)[order(Matrix::rowMeans(counts(sce),na.rm=TRUE),
                                decreasing=TRUE)[seq_len(ng)]]
    }else{
      g <- tryCatch({
        cl.means <- as.matrix(assay(scuttle::sumCountsAcrossCells(counts(sce), clusters)))
        g <- unique(as.numeric(t(apply(cl.means, 2, FUN=function(x){
          order(x, decreasing=TRUE)[seq_len(nfeatures)]
        }))))[seq_len(ng)]
        row.names(sce)[g]
      }, error=function(e){
        g <- row.names(sce)[order(Matrix::rowMeans(counts(sce),na.rm=TRUE),
                                  decreasing=TRUE)[seq_len(ng)]]
      })
    }
  }
  if(ng==nfeatures) return(g)
  mm <- scran::findMarkers(sce, groups=clusters, test.type="binom", assay.type="counts")
	mm <- dplyr::bind_rows(lapply(mm, FUN=function(x){
	  x <- x[x$FDR<FDR.max,]
	  data.frame(gene=row.names(x), Top=x$Top, FDR=x$FDR, stringsAsFactors=FALSE)
	}), .id = "cluster")
	g2 <- unique(c(g,mm$gene))
	if(length(g2)<nfeatures) return(g2)
	i <- nfeatures/(2*length(unique(clusters)))
	while(length(g <- unique(c(g,mm$gene[mm$Top<=i])))<nfeatures) i<-i+1
	return(head(g,nfeatures))
}


#' mockDoubletSCE
#'
#' Creates a mock random single-cell experiment object with doublets
#'
#' @param ncells A positive integer vector indicating the number of cells per
#' cluster (min 2 clusters)
#' @param ngenes The number of genes to simulate. Ignored if `mus` is given.
#' @param mus A list of cluster averages.
#' @param dbl.rate The doublet rate
#' @param only.heterotypic Whether to create only heterotypic doublets
#'
#' @return A SingleCellExperiment object, with the colData columns `type`
#' indicating whether the cell is a singlet or doublet, and `cluster`
#' indicating from which cluster (or cluster combination) it was simulated.
#'
#' @export
#' @import SingleCellExperiment
#' @importFrom stats rnorm rpois
#' @examples
#' sce <- mockDoubletSCE()
mockDoubletSCE <- function(ncells=c(200,300), ngenes=200, mus=NULL,
                           dbl.rate=0.1, only.heterotypic=TRUE){
  if(length(ncells)<2)
    stop("ncells should be a positive integer vector of length >=2")
  if(is.null(names(ncells))) names(ncells)<-paste0("cluster",seq_along(ncells))
  if(is.null(mus)){
    mus <- lapply(ncells, FUN=function(x) 2^rnorm(ngenes))
  }else{
    if(!is.list(mus) || length(mus)!=length(ncells) ||
       length(unique(lengths(mus)))!=1)
      stop("If provided, `mus` should be a list of length equal to that of ",
           "`ncells`, with each slot containing a numeric vector of averages ",
           "for each gene.")
    names(mus) <- names(ncells)
  }

  # non-doublets
  counts <- do.call(cbind, lapply(seq_along(mus), FUN=function(i){
    matrix(rpois(ncells[[i]]*ngenes, mus[[i]]), nrow=ngenes)
  }))
  sce <- SingleCellExperiment(list(counts=counts),
                              colData=data.frame(type="singlet",
                                           origin=rep(names(ncells), ncells)))

  # doublets
  expected <- getExpectedDoublets(sce$origin, dbl.rate,
                                  only.heterotypic=only.heterotypic)
  simdbl <- rpois(length(expected), expected)
  if(sum(simdbl)>0){
    mus <- lapply(strsplit(names(expected),"+",fixed=TRUE), FUN=function(x){
      mus[[x[[1]]]]+mus[[x[[2]]]]
    })
    doublets <- do.call(cbind, lapply(seq_along(mus), FUN=function(i){
      matrix(rpois(simdbl[[i]]*ngenes, mus[[i]]), nrow=ngenes)
    }))
    sce2 <- SingleCellExperiment( list(counts=doublets),
                                  colData=data.frame(type="doublet",
                                          origin=rep(names(expected), simdbl)))
    sce <- cbind(sce, sce2)
  }

  # set original cluster for homotypic doublets
  sce$cluster <- factor(sce$origin, c(names(ncells),names(expected)))
  names(n) <- n <- levels(sce$cluster)
  n[paste(names(ncells),names(ncells),sep="+")] <- names(ncells)
  levels(sce$cluster) <- as.character(n)
  sce$cluster <- droplevels(sce$cluster)

  colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  row.names(sce) <- paste0("gene",seq_len(ngenes))
  sce$type <- factor(sce$type, c("singlet","doublet"))
  sce
}


.checkColArg <- function(sce, x, acceptNull=TRUE){
  arg <- deparse(substitute(x))
  if(is.null(x)){
    if(!acceptNull) stop("Missing argument `",arg,"`!")
    return(NULL)
  }
  if(is.character(x) && length(x)==1){
    if(!(x %in% colnames(colData(sce))))
      stop("Could not find `", arg, "` column in colData!")
    x <- colData(sce)[[x]]
  }else if(length(x)!=ncol(sce)){
    stop("`",arg,"` should have a length equal to the number of columns in ",
         "`sce`.")
  }
  x
}
.checkPropArg <- function(x, acceptNull=TRUE){
  arg <- deparse(substitute(x))
  if( is.null(x) && acceptNull ) return(NULL)
  if(is.null(x) || length(x)!=1 || !is.numeric(x) || x>1 || x<0)
    stop("`",arg,"` should be a positive value between 0 and 1.")
  x
}

#' cxds2
#'
#' Calculates a coexpression-based doublet score using the method developed by
#' \href{https://doi.org/10.1093/bioinformatics/btz698}{Bais and Kostka 2020}.
#' This is the original implementation from the
#' \href{https://www.bioconductor.org/packages/release/bioc/html/scds.html}{scds}
#' package, but enabling scores to be calculated for all cells while the gene
#' coexpression is based only on a subset (i.e. excluding known/artificial doublets).
#'
#' @param x A matrix of counts, or a `SingleCellExperiment` containing a
#' 'counts'
#' @param whichDbls The columns of `x` which are known doublets.
#' @param ntop The number of top features to keep.
#' @param binThresh The count threshold to be considered expressed.
#' @references \url{https://doi.org/10.1093/bioinformatics/btz698}
#' @return A cxds score or, if `x` is a `SingleCellExperiment`, `x` with an
#' added `cxds_score` colData column.
#' @export
#' @importFrom stats pbinom
#' @examples
#' sce <- mockDoubletSCE()
#' sce <- cxds2(sce)
cxds2 <- function(x, whichDbls=c(), ntop=500, binThresh=0){
  if(is(x,"SingleCellExperiment")){
    x$cxds_score <- cxds2(counts(x), whichDbls=whichDbls, ntop=ntop,
                          binThresh=binThresh)
    return(x)
  }
  Bp <- x <- x > binThresh
  ps <- Matrix::rowMeans(x)
  if(nrow(x)>ntop){
    hvg <- order(ps * (1 - ps), decreasing=TRUE)[seq_len(ntop)]
    Bp <- x <- x[hvg, ]
    ps <- ps[hvg]
  }
  if(length(whichDbls)>0) Bp <- Bp[,-whichDbls]
  prb <- outer(ps, 1 - ps)
  prb <- prb + t(prb)
  obs <- Bp %*% (1 - Matrix::t(Bp))
  obs <- obs + Matrix::t(obs)
  S <- suppressWarnings({
    stats::pbinom(as.matrix(obs) - 1, prob = prb, size=ncol(Bp),
        lower.tail=FALSE, log.p=TRUE)
  })
  if(any(w <- is.infinite(S))){
    smin <- min(S[!is.infinite(S)])
    S[S<smin] <- smin
  }
  s <- -Matrix::colSums(x * (S %*% x))
  s <- s - min(s)
  s/max(s)
}

#' @importFrom stats cor
.clustSpearman <- function(e, clusters, nMarkers=30){
  if(is.null(dim(clusters))){
    e2 <- e[,seq_along(clusters)]
    g <- seq_len(nrow(e2))
    if(nMarkers>0 & nMarkers<nrow(e)){
      suppressWarnings(mm <- scran::findMarkers(e2, groups=clusters, test.type="binom"))
      g <- unique(unlist(lapply(mm, FUN=function(x) row.names(x)[seq_len(nMarkers)])))
    }
    e2 <- scuttle::sumCountsAcrossCells(e2[g,], ids=clusters)
    clusters <- as.matrix(assay(e2))
  }else{
    if(ncol(clusters)>500)
      warning("You're using a very large `clustCor` matrix, are you sure that the",
              "columns are cell types, and not individual cells?")
    g <- intersect(row.names(e), row.names(clusters))
    if(length(g)<5){
      warning("Too few marker genes in common between `clustCor` and the data.",
        "No correlation will be used. Consider checking that the gene identifiers match.")
      return(NULL)
    }
    clusters <- clusters[g,]
  }
  cor(as.matrix(e[g,]), clusters, method="spearman")
}

# gets a global doublet rate from samples' doublet rates
.gdbr <- function(d, dbr=NULL){
  if(!is.null(dbr)){
    if(length(dbr)==1) return(dbr)
    stopifnot(!is.null(d$sample))
    sl <- as.numeric(table(d$sample, d$src=="real")[,2])
    return(sum(dbr*sl)/sum(sl))
  }
  if(is.null(d$sample)){
    sl <- sum(d$src=="real")
  }else{
    ## estimate a global doublet rate
    sl <- as.numeric(table(d$sample, d$src=="real")[,2])
  }
  dbr <- (0.01*sl/1000)
  sum(dbr*sl)/sum(sl)
}


#' propHomotypic
#'
#' Computes the proportion of pairs expected to be made of elements from the same cluster.
#'
#' @param clusters A vector of cluster labels
#'
#' @return A numeric value between 0 and 1.
#' @export
#'
#' @examples
#' clusters <- sample(LETTERS[1:5], 100, replace=TRUE)
#' propHomotypic(clusters)
propHomotypic <- function(clusters){
  p <- as.numeric(table(clusters)/length(clusters))
  p <- p %*% t(p)
  sum(diag(p))/sum(p)
}

.totExpectedHeteroDoublets <- function(clusters, dbr=NULL){
  if(is.factor(clusters)) clusters <- droplevels(clusters)
  ncells <- length(clusters)
  if(is.null(dbr)) dbr <- (0.01*ncells/1000)
  p <- as.numeric(table(clusters)/ncells)
  p <- (p %*% t(p)) * ncells*dbr
  diag(p) <- 0
  sum(p)
}


# procedure to 0-1 rescale scores across samples so that the mins, maxs, and trimmed mean
# of real and artificial cells is the same across samples; then doing quantile
# normalization between these fixed points
.rescaleSampleScores <- function(d, byType=is.null(q), q=NULL, what="score", newName=NULL,
                                 mode=c("quantile","linear")){
  mode <- match.arg(mode)
  if(what=="score" && is.null(newName)) d$perSampleScore <- d$score
  if(is.null(newName)){
    newName <- what
  }else{
    d[[newName]] <- NA_real_
  }
  if(is.null(q)){
    if(byType){
      wD <- which(d$type=="doublet")
      q <- t(vapply(split(d[[what]], d$sample), FUN.VALUE=numeric(2L), FUN=function(x){
        sort(c(mean(x[-wD], na.rm=TRUE, trim=0.5), mean(x[wD], na.rm=TRUE, trim=0.5))) }))
      colnames(q) <- c("singlet","doublet")
    }else{
      q <- as.matrix(vapply(split(d[[what]], d$sample), FUN.VALUE=numeric(1L),
                            FUN=mean, trim=0.5, na.rm=TRUE))
    }
  }else{
    if(byType) stop("Cannot scale by type with given quantiles.")
    if(is.null(dim(q))) q <- as.matrix(q)
  }
  if(mode=="quantile"){
    .qscale <- function(x, by){
      si <- split(seq_along(x),droplevels(by))
      sf <- lapply(si, ecdf)
      for(s in names(si)) x[si[[s]]] <- sf[[s]](x[si[[s]]])
      x
    }
  }else{
    .qscale <- function(x, by){
      mins <- vapply(split(x,by),FUN.VALUE=numeric(1),FUN=function(x){
        if(length(x)==0) return(0); min(x,na.rm=TRUE) })
      maxs <- vapply(split(x,by),FUN.VALUE=numeric(1),FUN=function(x){
        if(length(x)==0) return(1); max(x,na.rm=TRUE) })
      (x-mins[by])/(maxs[by]-mins[by])
    }
  }
  or.sc <- d[[what]]
  if(ncol(q)>1){
    qmeds <- apply(q,2,median)
    w <- which(or.sc <= q[d$sample,1])
    d[[newName]][w] <- .qscale(or.sc[w],by=d$sample[w])*qmeds[1]
    w <- which(or.sc > q[d$sample,1] & or.sc <=q[d$sample,2])
    d[[newName]][w] <- qmeds[1]+.qscale(or.sc[w],by=d$sample[w])*(qmeds[2]-qmeds[1])
    w <- which(or.sc > q[d$sample,2])
    d[[newName]][w] <- qmeds[2]+.qscale(or.sc[w],by=d$sample[w])*(1-qmeds[2])
    maxs <- vapply(split(d[[newName]],d$sample),FUN.VALUE=numeric(1),na.rm=TRUE,FUN=max)
    d[[newName]] <- d[[newName]]/maxs[d$sample]
  }else{
    thm <- median(q)
    ll <- split(or.sc, d$sample)
    ths <- q[d$sample,]
    w <- which(or.sc<=ths)
    d[[newName]][w] <- thm*d[[what]][w]/ths[w]
    w <- which(or.sc>ths)
    d[[newName]][w] <- thm + (d[[what]][w]-ths[w])*(1-thm)/(1-ths[w])
  }
  return(d)
}

#' aggregateFeatures
#'
#' Aggregates similar features (rows).
#'
#' @param x A integer/numeric (sparse) matrix, or a `SingleCellExperiment`
#' including a `counts` assay.
#' @param dims.use The PCA dimensions to use for clustering rows.
#' @param k The approximate number of meta-features desired
#' @param num_init The number of initializations used for k-means clustering.
#' @param use.mbk Logical; whether to use minibatch k-means (see
#' \code{\link[mbkmeans]{mbkmeans}}). If NULL, the minibatch approach will be
#' used if there are more than 30000 features.
#' @param use.subset How many cells (columns) to use to cluster the features.
#' @param use.TFIDF Logical; whether to use \link{TFIDF} normalization (instead
#' of standard normalization) to assess the similarity between features.
#' similarity
#' @param ... Passed to \code{\link[mbkmeans]{mbkmeans}}. Can for instance be
#' used to pass the `BPPARAM` argument for multithreading.
#'
#' @return An aggregated version of `x` (either an array or a
#' `SingleCellExperiment`, depending on the input).
#'
#' @importFrom scuttle logNormCounts
#' @importFrom BiocSingular runPCA IrlbaParam
aggregateFeatures <- function(x, dims.use=seq(2L,12L), k=1000, num_init=2,
                               use.mbk=NULL, use.subset=5000, use.TFIDF=TRUE,
                              ...){
  xo <- x
  if(ncol(x)>use.subset){
    if(is(x,"SingleCellExperiment")){
      cs <- Matrix::colSums(counts(x))
    }else{
      cs <- Matrix::colSums(x)
    }
    x <- x[,order(cs,decreasing=TRUE)[seq_len(use.subset)]]
  }
  if(use.TFIDF){
    if(is(x,"SingleCellExperiment")) x <- counts(x)
    x <- TFIDF(x)
  }else{
    if(is(x,"SingleCellExperiment")){
      if("logcounts" %in% assayNames(x)) x <- scuttle::logNormCounts(x)
      x <- logcounts(x)
    }
    x <- t(normalizeCounts(t(x)))
  }
  pca <- runPCA(x, BSPARAM=IrlbaParam(), center=FALSE,
                rank=max(dims.use))$x[,dims.use]
  if(is.null(use.mbk)) use.mbk <- nrow(x) > 30000
  if(use.mbk && suppressWarnings(require("mbkmeans", quietly=TRUE))){
    fc <- mbkmeans:::mbkmeans(t(pca), k, num_init=num_init, ...)$Clusters
  }else{
    fc <- kmeans(pca, k, nstart=num_init, iter.max=100)$cluster
  }
  x <- scuttle::sumCountsAcrossFeatures(xo, fc)
  row.names(x) <- paste0("feat",seq_len(nrow(x)))
  if(is(xo,"SingleCellExperiment")){
    x <- SingleCellExperiment(list(counts=x), colData=colData(xo),
                              reducedDims=reducedDims(xo))
  }
  x
}

#
#' TFIDF
#'
#' The Term Frequency - Inverse Document Frequency (TF-IDF) normalization, as
#' implemented in Stuart & Butler et al. 2019.
#'
#' @param x The matrix of occurrences
#' @param sf Scaling factor
#'
#' @return An array of same dimensions as `x`
#' @export
#' @importFrom Matrix tcrossprod Diagonal rowSums colSums
#'
#' @examples
#' m <- matrix(rpois(500,1),nrow=50)
#' m <- TFIDF(m)
TFIDF <- function(x, sf=10000){
  if(!is(x,"sparseMatrix")) x <- as(x, "sparseMatrix")
  tf <- Matrix::tcrossprod(x, Diagonal(x=1L/Matrix::colSums(x)))
  idf <- ncol(x)/Matrix::rowSums(x)
  x <- log1p(sf*(Diagonal(length(idf), x=idf) %*% tf))
  x[is.na(x)] <- 0
  x
}
