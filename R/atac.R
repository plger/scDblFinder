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



#' amuletFromCounts
#' 
#' A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021), based on tile/peak 
#' counts. Note that this is only a fast approximation to the original Amulet 
#' method; for an equivalent implementation, see 
#' \code{\link{amuletFromFragments}}.
#'
#' @param x A `SingleCellExperiment` object, or a matrix of counts with cells
#' as columns. If the rows represent peaks, it is recommended to limite their
#' width (see details).
#' @param feature.na2.min the minimum number of cells in which a feature has 
#' more than 2 reads in order for the feature to be kept
#' @param maxWidth the maximum width for a feature to be included. This is 
#' ignored unless `x` is a `SingleCellExperiment` with `rowRanges`.
#' @param feature.q.threshold significance threshold for features to be excluded
#' @param correction.method multiple testing corrected method (passed to 
#' `p.adjust`)
#'
#' @return If `x` is a `SingleCellExperiment`, returns the object with an 
#' additional `amulet.q` colData column. Otherwise returns a vector of the 
#' amulet doublet q-values for each cell.
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
#' x <- amulet(x)
#' table(call=x$amulet.q<0.05, truth=x$type)
amuletFromCounts <- function(x, feature.na2.min=max(2,0.01*ncol(x)), 
                             maxWidth=500L, feature.q.threshold=0.01, 
                             correction.method="BH"){
  if(is(x, "SingleCellExperiment")){
    if(!is.null(ranges(x)) && 
       !all(sum(IRanges::width(ranges(x)))==0L) ){
      y <- counts(x)[which(IRanges::width(ranges(x))<=maxWidth),]
    }else{
      y <- counts(x)
    }
    x$amulet.q <- amuletFromCounts( y, feature.na2.min=feature.na2.min, 
                                    feature.q.threshold=feature.q.threshold, 
                                    correction.method=correction.method )
    return(x)
  }
  x <- x>2L
  rs <- rowSums(x)
  x <- x[rs>=feature.na2.min,]
  rs <- rs[rs>=feature.na2.min]
  if(feature.q.threshold<1){
    rp <- p.adjust(ppois(rs, lambda = mean(rs), lower.tail=FALSE))
    x <- x[rp>feature.q.threshold,]
  }
  
  cs <- colSums(x)
  p.adjust(ppois(cs, lambda = mean(cs), lower.tail=FALSE), 
           method=correction.method)
}


#' amuletFromFragments
#' 
#- A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021)
#'
#' @param x The path to a fragments file
#' @param barcodes Optional character vector of cell barcodes to consider
#' @param minFrags Minimum number of fragments for a barcode to be 
#' considered. If `uniqueFrags=TRUE`, this is the minimum number of unique 
#' fragments. Ignored if `barcodes` is given.
#' @param above The number of overlaps above which to count
#' @param regionsToExclude A GRanges of regions to exclude (e.g. repeats)
#' @param uniqueFrags Logical; whether to use only unique fragments.
#' @param verbose Logical; whether to print progress messages.
#'
#' @details This implementation is relatively fast (except for reading in the 
#' data) but it has a large memory footprint (easily goes above 20GB) since 
#' the overlaps are performed in memory.
#'
#' @return A vector of the number of overlaps above `above`, named with cell 
#' barcodes
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors splitAsList
#' @importFrom IRanges overlapsAny
#' @export
amuletFromFragments <- function(x, barcodes=NULL, minFrags=500L, above=2L,
                                regionsToExclude=NULL, uniqueFrags=FALSE,
                                verbose=TRUE){
  stopifnot(is.null(barcodes) || is.character(barcodes))
  if(is.character(x) && length(x)==1){
    if(!file.exists(x)) stop("x should be a fragment file!")
    if(verbose) message("Reading fragment file...")
    gr <- .import.bed(x)
  }else{
    if(!is(x, "GRanges") || !all(c("name","score") %in% colnames(mcols(x))))
      stop("`x` should either be a path to a fragments file, or a GRanges ",
           "the 'name' and 'score' columns (respectively the cell barcode and",
           " the counts.")
    gr <- x
  }
  gr$name <- as.factor(gr$name)
  if(!is.null(regionsToExclude)){
    if(length(regionsToExclude)==1 && file.exists(regionsToExclude)){
      if(verbose) message("Reading regions to exclude...")
      regionsToExclude <- .import.bed(regionsToExclude)
    }else{
      stopifnot(is.null(regionsToExclude) || is(regionsToExclude, "GRanges"))
    }
    gr <- gr[!overlapsAny(gr, regionsToExclude)]
  } 
  if(verbose) message("Splitting and subsetting barcodes...")
  if(!uniqueFrags) gr <- gr[rep(seq_along(gr),as.integer(gr$score))]
  gr <- GenomicRanges::split(gr, gr$name)
  if(!is.null(barcodes)){
    if(length(barcodes)==1 && file.exists(barcodes))
      barcodes <- readLines(barcodes)
    if((mis <- length(setdiff(barcodes, names(gr))))>0)
      warning("Some barcodes (", mis, " or ",round(100*mis/length(gr),1),"%)",
              " are missing from the fragments file!")
    gr <- gr[barcodes]
    if(length(gr)==0) stop("No barcode found!")
  }else{
    gr <- gr[lengths(gr)>=minFrags]
  }
  if(verbose) message("Obtaining overlaps...")
  nU <- lengths(gr)
  gr <- unlist(reduce(gr,with.revmap=TRUE))
  nA2 <- sum(splitAsList(lengths(gr$revmap)>above,names(gr)))
  if(verbose) message("Stats and wrap-up...")
  d <- data.frame( row.names=names(nA2),
                   nUnique=nU[names(nA2)],
                   nAbove=nA2,
                   p.value=ppois(nA2, mean(nA2), lower.tail=FALSE))
  d$q.value <- p.adjust(d$p.value, method="BH")
  d
}
