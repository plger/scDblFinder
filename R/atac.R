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
aggregateFeatures <- function(x, dims.use=seq(2L,12L), k=1000, num_init=3,
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
  if(use.mbk && suppressWarnings(requireNamespace("mbkmeans", quietly=TRUE))){
    fc <- mbkmeans::mbkmeans(t(pca), k, num_init=num_init, ...)$Clusters
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

#' getFragmentOverlaps
#' 
#' Count the number of overlapping fragments.
#' 
#' @param x The path to a fragments file, or a GRanges object containing the
#' fragments (with the `name` column containing the barcode, and optionally 
#' the `score` column containing the count).
#' @param barcodes Optional character vector of cell barcodes to consider
#' @param minFrags Minimum number of fragments for a barcode to be 
#' considered. If `uniqueFrags=TRUE`, this is the minimum number of unique 
#' fragments. Ignored if `barcodes` is given.
#' @param regionsToExclude A GRanges of regions to exclude. As per the original
#'   Amulet method, we recommend excluding repeats, as well as sex and 
#'   mitochondrial chromosomes.
#' @param uniqueFrags Logical; whether to use only unique fragments.
#' @param maxFragSize Integer indicating the maximum fragment size to consider
#' @param removeHighOverlapSites Logical; whether to remove sites that have
#'   more than two reads in unexpectedly many cells.
#' @param fullInMemory Logical; whether to process all chromosomes together.
#'   This will speed up the process but at the cost of a very high memory 
#'   consumption (as all fragments will be loaded in memory). This is anyway the
#'   default mode when `x` is not Tabix-indexed.
#' @param verbose Logical; whether to print progress messages.
#'
#' @details When used on normal (or compressed) fragment files, this 
#' implementation is relatively fast (except for reading in the data) but it 
#' has a large memory footprint since the overlaps are performed in memory. It 
#' is therefore recommended to compress the fragment files using bgzip and index
#' them with Tabix; in this case each chromosome will be read and processed 
#' separately, leading to a considerably lower memory footprint.
#'
#' @return A data.frame with counts and overlap statistics for each barcode.
#' 
#' @importFrom GenomicRanges reduce GRanges
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges overlapsAny IRanges
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom rtracklayer import
#' @export
getFragmentOverlaps <- function(x, barcodes=NULL, regionsToExclude=GRanges(
                                  c("M","chrM","MT","X","Y","chrX","chrY"),
                                  IRanges(1L,width=10^8)), 
                                minFrags=500L, uniqueFrags=TRUE, 
                                maxFragSize=1000L, removeHighOverlapSites=TRUE,
                                fullInMemory=FALSE, verbose=TRUE){
  # prepare inputs
  stopifnot(is.null(barcodes) || is.character(barcodes))
  if(!is.null(barcodes) && length(barcodes)==1 && file.exists(barcodes))
    barcodes <- readLines(barcodes) # assume barcodes to be a text file
  if(!is.null(regionsToExclude)){
    if(length(regionsToExclude)==1 && file.exists(regionsToExclude)){
      if(verbose) message("Reading regions to exclude...")
      regionsToExclude <- rtracklayer::import(regionsToExclude)
    }else{
      stopifnot(is.null(regionsToExclude) || is(regionsToExclude, "GRanges"))
    }
    regionsToExclude <- sort(regionsToExclude)
  }
  # prepare empty output for returns
  emptyOutput <- data.frame(row.names=character(0), nFrags=integer(0), 
                            uniqFrags=integer(0), nAbove2=integer(0))
  
  if(is.character(x) && length(x)==1){
    if(!file.exists(x)) stop("x should be a fragment file!")
    if(!fullInMemory &&
       is(tf <- tryCatch(TabixFile(x), error=function(e) NULL), "TabixFile")){
      if(verbose) message("Reading Tabix-indexed fragment file and ",
                          "computing overlaps...")
      x <- lapply(seqnamesTabix(tf), FUN=function(x){ # count for each chr
        if(verbose) cat(paste0(x,", "))
        getFragmentOverlaps(
          rtracklayer::import(tf, format="bed", 
                              which=GRanges(x, IRanges(1,10^8))), 
          barcodes, regionsToExclude=regionsToExclude, 
          minFrags=0.00001, uniqueFrags=uniqueFrags, 
          verbose=FALSE, maxFragSize=maxFragSize
        )
      })
      if(verbose){
        cat("\n")
        message("Merging...")
      }
      x <- x[unlist(lapply(x,nrow))>0]
      if(length(x)==0) return(emptyOutput)

      if(is.null(barcodes)){
        barcodes <- rowsum(
          unlist(lapply(x, FUN=function(x) x$nFrags)), 
                        unlist(lapply(x, row.names)))[,1]
        barcodes <- names(barcodes[barcodes>=minFrags])
      }
      return(as.data.frame(Reduce("+", lapply(x, FUN=function(x){
        x <- as.matrix(x[barcodes,])
        x[is.na(x)] <- 0
        x
      })), row.names=barcodes))
    }else{
      if(!fullInMemory){
        warning("Fragment file is not tabix-indexed, requiring the",
                "whole file to be imported in memory.")
      }else if(verbose){
        message("Reading full fragments...")
      }
      gr <- rtracklayer::import(x, format="bed")
    }
  }else{
    if(!is(x, "GRanges") || !("name" %in% colnames(mcols(x))))
      stop("`x` should either be a path to a fragments file, or a GRanges ",
           "with the 'name' column containing the cell barcode (and optionally
           the 'score' column containing the counts).")
    gr <- x
  }
  gr <- gr[(width(gr)<=maxFragSize),]
  gr$name <- as.factor(gr$name)
  if(!is.null(regionsToExclude))
    gr <- gr[!overlapsAny(gr, regionsToExclude)]
  
  if(verbose) message("Splitting and subsetting barcodes...")
  uniqFrags <- table(gr$name)
  if(minFrags<1L & minFrags>0L) minFrags <- round(minFrags*length(gr))
  if(is.null(barcodes)){
    uniqFrags <- uniqFrags[uniqFrags>=minFrags]
    gr <- gr[gr$name %in% names(uniqFrags)]
  }else{
    if((mis <- length(setdiff(barcodes, names(uniqFrags))))>0)
      if(verbose)
        warning("Some barcodes (", mis, " or ",round(100*mis/length(gr),1),"%)",
                " are missing from the fragments file!")
    uniqFrags <- uniqFrags[intersect(names(uniqFrags), barcodes)]
    gr <- gr[gr$name %in% names(uniqFrags)]
  }
  if(length(gr)==0) return(emptyOutput)
  if(!uniqueFrags && !is.null(score(gr))){
    i <- rep(seq_along(gr), as.integer(score(gr)))
    gr <- GenomicRanges::split(granges(gr)[i], gr$name[i])
    rm(i)
  }else{
    gr <- GenomicRanges::split(granges(gr), gr$name)
  }
  if(verbose) message("Obtaining overlaps...")
  nFrags <- lengths(gr)
  gr <- unlist(reduce(gr, with.revmap=TRUE, min.gapwidth=0L))
  gr <- gr[lengths(gr$revmap)>2L]
  if(removeHighOverlapSites){
    hor <- reduce(gr, with.revmap=TRUE, min.gapwidth=0L)
    hlw <- width(hor)
    hor <- hor$revmap
    hl <- lengths(hor)
    if(FALSE){ # correct for width
      p <- tryCatch({
        mod <- MASS::glm.nb(hl~log(hlw), maxit=1000)
        pnbinom(hl, mu=pmax(predict(mod),0), lower.tail=FALSE,
                size=mod$theta)
      }, error=function(e){
        mod <- glm(hl~log(hlw), family="poisson")
        ppois(hl, lambda=pmax(predict(mod),0), lower.tail=FALSE)
      })
    }else{
      p <- tryCatch({
        pnbinom(hl, mu=mean(hl), lower.tail=FALSE,
                size=MASS::theta.ml(hl, mean(hl)))
      }, error=function(e){
        ppois(hl, mean(hl), lower.tail=FALSE)
      })
    }
    if(length(indices2remove <- unlist(hor[which(p<0.01)]))>0)
      gr <- gr[-indices2remove]
  }
  if(length(gr)==0) return(emptyOutput)
  tt <- table(names(gr))
  nAbove2 <- setNames(rep(0L,length(uniqFrags)), names(uniqFrags))
  nAbove2[names(tt)] <- as.integer(tt)
  rm(gr,tt,hor)
  gc(verbose=FALSE)
  return(data.frame(row.names=names(uniqFrags), 
                    nFrags=as.integer(nFrags[names(uniqFrags)]),
                    uniqFrags=as.integer(uniqFrags), nAbove2=nAbove2))
}





.overlapsFromCounts <- function(x, maxWidth=500L, ...){
  if(is(x, "SingleCellExperiment")){
    if(!is.null(ranges(x)) && 
       !all(sum(IRanges::width(ranges(x)))==0L) ){
      y <- counts(x)[which(IRanges::width(ranges(x))<=maxWidth),]
    }else{
      y <- counts(x)
    }
  }else if(is.matrix(x)){
    y <- x
  }else{
    stop("Object should be a count matrix or a SingleCellExperiment")
  }
  libsizes <- colSums(y)
  y <- y>2L
  rs <- rowSums(y)
  rp <- p.adjust(ppois(rs, lambda=mean(rs), lower.tail=FALSE))
  y <- y[rp>0.01,]
  data.frame(row.names=colnames(y), nFrags=libsizes, nAbove2=colSums(y))
}

#' amuletFromCounts
#' 
#' A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021), based on tile/peak 
#' counts. Note that this is only a fast approximation to the original Amulet 
#' method, and performs considerably worse; for an equivalent implementation, 
#' see \code{\link{amuletFromFragments}}.
#'
#' @param x A `SingleCellExperiment` object, or a matrix of counts with cells
#' as columns. If the rows represent peaks, it is recommended to limite their
#' width (see details).
#' @param maxWidth the maximum width for a feature to be included. This is 
#' ignored unless `x` is a `SingleCellExperiment` with `rowRanges`.
#'
#' @return If `x` is a `SingleCellExperiment`, returns the object with an 
#' additional `amuletFromCounts.q` colData column. Otherwise returns a vector of
#'  the amulet doublet q-values for each cell.
#' 
#' @details 
#' The rationale for the amulet method is that a single diploid cell should not 
#' have more than two reads covering a single genomic location, and the method 
#' looks for cells enriched with sites covered by more than two reads.
#' If the method is applied on a peak-level count matrix, however, larger peaks
#' can however contain multiple reads even though no single nucleotide is 
#' covered more than once. Therefore, in such case we recommend to limit the 
#' width of the peaks used for this analysis, ideally to maximum twice the upper
#' bound of the fragment size. For example, with a mean fragment size of 250bp 
#' and standard deviation of 125bp, peaks larger than 500bp are very likely to 
#' contain non-overlapping fragments, and should therefore be excluded using the
#' `maxWidth` argument.
#' 
#' @seealso \code{\link{amuletFromFragments}}
#' 
#' @importFrom IRanges width
#' @importFrom SummarizedExperiment ranges
#' @export
#' @examples
#' x <- mockDoubletSCE()
#' x <- amuletFromCounts(x)
#' table(call=x$amulet.q<0.05, truth=x$type)
amuletFromCounts <- function(x, maxWidth=500L){
  d <- .overlapsFromCounts(x, maxWidth)
  q <- p.adjust(ppois(d$nAbove2, lambda=mean(d$nAbove2), lower.tail=FALSE), 
                method="fdr")
  if(is(x,"SummarizedExperiment")){
    x$amuletFromCounts.q <- q
    return(x)
  }
  q
}



#' amulet
#' 
#- A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021)
#'
#' @param x The path to a fragments file, or a GRanges object containing the
#' fragments (with the `name` column containing the barcode, and the `score`
#' column containing the count).
#' @param ... Any argument to \code{\link{getFragmentOverlaps}}.
#'
#' @details When used on normal (or compressed) fragment files, this 
#' implementation is relatively fast (except for reading in the data) but it 
#' has a large memory footprint since the overlaps are performed in memory. It 
#' is therefore recommended to compress the fragment files using bgzip and index
#' them with Tabix; in this case each chromosome will be read and processed 
#' separately, leading to a considerably lower memory footprint.
#'
#' @return A vector of the number of overlaps above `above`, named with cell 
#' barcodes
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges overlapsAny
#' @export
amulet <- function(x, ...){
  d <- getFragmentOverlaps(x, ...)
  d$p.value <- ppois(d$nAbove2, mean(d$nAbove2), lower.tail=FALSE)
  d$q.value <- p.adjust(d$p.value, method="BH")
  d
}


#' amulet2
#' 
#- A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021)
#'
#' @param x The path to a fragments file, or a GRanges object containing the
#' fragments (with the `name` column containing the barcode, and the `score`
#' column containing the count).
#' @param ... Any argument to \code{\link{getFragmentOverlaps}}.
#'
#' @details When used on normal (or compressed) fragment files, this 
#' implementation is relatively fast (except for reading in the data) but it 
#' has a large memory footprint since the overlaps are performed in memory. It 
#' is therefore recommended to compress the fragment files using bgzip and index
#' them with Tabix; in this case each chromosome will be read and processed 
#' separately, leading to a considerably lower memory footprint.
#'
#' @return A vector of the number of overlaps above `above`, named with cell 
#' barcodes
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges overlapsAny
#' @export
amulet2 <- function(x, ...){
  d <- getFragmentOverlaps(x, ...)
  
  d$p.value <- ppois(d$nAbove2, mean(d$nAbove2), lower.tail=FALSE)
  d$q.value <- p.adjust(d$p.value, method="BH")
  d
}


#' conformalAmulet
#'
#' @param x A data.frame with the columns `nFrags` (the library size),  
#' `nAbove2` (the number of loci with more than 2 reads) and optionally 
#' `sample`; or paths to fragment file(s); or (not recommended) a count matrix 
#' or \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param conformal Logical, whether to use conformal inference
#' @param verbose Logical; whether to print progress messages.
#' @param conf.alpha Alpha level for the conformal inference
#' @param n.cyc Number of \link[gamlss]{gamlss} cycles
#' @param fitSeparately Logical; whether to fit multiple samples separately
#' @param local Logical; whether to use local weighting
#' @param dbl.alpha The q-value threshold for doublets (during iterations)
#' @param maxIt The maximum number of iterations to run
#' @param correction.method The multiple testing correction method (see 
#' \code{\link{p.adjust}})
#' @param BPPARAM Multithreading parameters for calculating overlaps when 
#' reading from fragment file(s), and for fitting when samples are fit 
#' separately
#'
#' @return A data.frame of cell features, including doublet p-values
#' @export
#'
#' @examples
conformalAmulet <- function(x, conformal=TRUE, local=TRUE, fitSeparately=TRUE,
                            conf.alpha=0.005, n.cyc=500, dbl.alpha=0.05, 
                            dbr=NULL, maxIt=5, verbose=TRUE, qth.start=0.9,
                            correction.method="BH", BPPARAM=SerialParam(), ...){
  if(is(x, "GRanges") || is(x, "SingleCellExperiment") || is.matrix(x))
    x <- list(x)
  if(is.character(x) || (is.list(x) && 
                         all(unlist(lapply(x,class2="GRanges",FUN=is))))){
    d <- lapply(x, FUN=function(x){
      if(verbose & is.character(x)) message(x)
      getFragmentOverlaps(x, barcodes=barcodes, BPPARAM=BPPARAM,
                 minFrags=minFrags, uniqueFrags=uniqueFrags, 
                 regionsToExclude=regionsToExclude, verbose=verbose)
    })
    d <- dplyr::bind_rows(d, .id="sample")
  }else if(is(x[[1]], "SingleCellExperiment") || is.matrix(x[[1]])){
    d <- lapply(x, FUN=.overlapsFromCounts,
                correction.method=correction.method)
    d <- dplyr::bind_rows(d, .id="sample")
  }else{
    d <- as.data.frame(x)
  }
  stopifnot(is.data.frame(d) && all(c("nFrags","nAbove2") %in% colnames(d)))
  if(uniqueFrags && "uniqFrags" %in% colnames(d)){
    d$x <- d$uniqFrags
  }else{
    d$x <- d$nFrags
  }
  d$y <- d$nAbove2
  d <- d[d$x>=minFrags,]
  
  if(conformal){
    ff <- function(dat, newdat=NULL)
      .getConf(dat, newdat, conf.alpha=conf.alpha, n.cyc=n.cyc, local=local,
               adjmethod=correction.method)
  }else{
    ff <- function(dat, newdat=NULL) .simpleFragNB(dat, newdat, maxit=n.cyc,
                                                   adjmethod=correction.method)
  }
  if(!fitSeparately || is.null(d$sample)){
    sspl <- rep(1L,nrow(d))
  }else{
    sspl <- d$sample
  }
  if(verbose & maxIt>1) cat("Iterative fit ")
  d <- lapply(split(d, sspl), FUN=function(x){
    if(is.null(dbr)){
      nc <- nrow(x)
      if(!is.null(x$sample) && length(unique(x$sample))>1)
        nc <- table(x$sample)
      dbr <- weighted.mean(0.01*nc/1000, nc/sum(nc))
    }
    maxNDbl <- (1.1*dbr)*nrow(x)
    x[,c("expected","p","adj.p")] <- ff(x[x$y<=quantile(x$y,qth.start),],x)
    wKeepOld <- seq_len(nrow(x))
    # wKeep <- which(x$adj.p>dbl.alpha | rank(x$p, ties.method="min")>maxNDbl)
    wKeep <- which(x$p>dbl.alpha | rank(x$p, ties.method="min")>maxNDbl)
    i <- 1L
    if(verbose & maxIt>1) cat(".")
    while(i<maxIt && !identical(wKeepOld, wKeep)){
      wKeepOld <- wKeep
      x[,c("expected","p","adj.p")] <- ff(x[wKeep,], x)
      wKeep <- which(x$p>dbl.alpha | rank(x$p, ties.method="min")>maxNDbl)
      i <- i + 1L
      if(verbose) cat(".")
    }
    if(verbose & maxIt>1) cat("\n")
    x
  })
  d <- do.call(rbind, d)
  d$x <- d$y <- NULL
  d
}

.simpleFragNB <- function(dat, newdat=NULL, maxit=100, adjmethod="BH"){
  if(is.null(newdat)) newdat <- dat
  mod <- tryCatch(MASS::glm.nb(y~log(x), data=dat, maxit=maxit),
                  error=function(e)
                    glm(y~log(x), family="poisson", data=dat))
  expected <- predict(mod, newdata=newdat)
  expected[expected<0] <- 0
  p <- rep(1, length(expected))
  w0 <- which(newdat$y>0)
  if(is.null(mod$theta)){
    p[w0] <- ppois(newdat$y[w0], expected[w0], lower.tail=FALSE)
  }else{
    p[w0] <- pnbinom(newdat$y[w0], size=mod$theta, mu=expected[w0], 
                     lower.tail=FALSE)
  }
  data.frame(expected=expected, p=p, adj.p=p.adjust(p, method=adjmethod))
}
  
.getGamlss <- function(dat, xpred=NULL, n.cyc=500){
  dat <- dat[,c("x","y")]
  h <- gamlss::gamlss(y~ps(x), family=NBI, data=dat, 
                      sigma.formula=~cs(x, c.spar=c(-1.5, 2.5)),
                      control=gamlss.control(n.cyc=n.cyc, trace=FALSE))
  if(is.null(xpred)){
    xpred <- dat
  }else{
    if(!is.data.frame(xpred)) xpred <- data.frame(x=xpred)
    xpred <- xpred[,intersect(c("x","y"),colnames(xpred)),drop=FALSE]
  }
  mu <- predict(h, newdata=xpred, type="response", data=dat)
  sigma <- predict(h, newdata=xpred, type="response", what="sigma", data=dat)
  list(mu=mu, sigma=sigma, mod=h)
}

.conformalSemiwidth <- function(d, alpha, local=TRUE, ...){
  i <- sample.int(nrow(d), floor(nrow(d)/2))
  gfit <- .getGamlss(d[i,], d[-i,], ...)
  resi <- abs(d$y[-i]-gfit[[1]])
  kOut <- ceiling((length(i)+1)*(1-alpha))
  if(local){
    stdevgamlss <- sqrt(gfit[[1]]+gfit[[2]]*gfit[[1]]^2)
    resi <- resi/stdevgamlss
  }
  resi[order(resi)][kOut]
}

#' @importFrom gamlss gamlss pNBI NBI cs gamlss.control
.getConf <- function(d, d2=NULL, conf.alpha=0.005, n.cyc=500, local=TRUE, 
                     adjmethod="BH"){
  capture.output({
    gfit <- .getGamlss(d, d2, n.cyc=n.cyc)
    semi <- .conformalSemiwidth(d, alpha=conf.alpha, local=local)
  })
  conformal <- gfit$mu+gfit$sigma*semi
  if(is.null(d2)) d2 <- d
  sigint <- max(10^-3,as.numeric(gfit$mod$sigma.coefficients))
  p <- pNBI(d2$y, mu=conformal, sigma=sigint, lower.tail=FALSE)
  #print(list(semi=semi, sigint=sigint))
  data.frame(expected=conformal, p=p, adj.p=p.adjust(p,method=adjmethod))
}
