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
