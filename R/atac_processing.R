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
#' @param minCount The minimum number of counts for a region to be included.
#' @param use.mbk Logical; whether to use minibatch k-means (see
#' \code{\link[mbkmeans]{mbkmeans}}). If NULL, the minibatch approach will be
#' used if there are more than 30000 features.
#' @param use.subset How many cells (columns) to use to cluster the features.
#' @param norm.fn The normalization function to use on the un-clustered data (a
#'   function taking a count matrix as a single argument and returning a matrix
#'   of the same dimensions). \link{TFIDF} by default.
#' @param twoPass Logical; whether to perform the procedure twice, so in the 
#'   second round cells are aggregated based on the meta-features of the first 
#'   round, before re-clustering the features. Ignored if the dataset has fewer
#'   than `use.subset` cells.
#'   
#' @param ... Passed to \code{\link[mbkmeans]{mbkmeans}}. Can for instance be
#' used to pass the `BPPARAM` argument for multithreading.
#'
#' @return An aggregated version of `x` (either an array or a
#' `SingleCellExperiment`, depending on the input). If `x` is a 
#' `SingleCellExperiment`, the feature clusters will also be stored in 
#' `metadata(x)$featureGroups`
#'
#' @importFrom scuttle logNormCounts
#' @importFrom BiocSingular runPCA IrlbaParam
#' @export
aggregateFeatures <- function(x, dims.use=seq(2L,12L), k=1000, num_init=3,
                              use.mbk=NULL, use.subset=20000, minCount=1L,
                              norm.fn=TFIDF, twoPass=FALSE, ...){
  xo <- x
  
  if(ncol(x)>use.subset){
    if(is(x,"SingleCellExperiment")){
      cs <- Matrix::colSums(counts(x))
    }else{
      cs <- Matrix::colSums(x)
    }
    # get rid of the cells with low libsize
    x <- x[,head(order(cs,decreasing=TRUE),
                 min(mean(use.subset,ncol(x)),2L*use.subset))]
    # if needed, sample randomly the remaining
    if(ncol(x)>use.subset)
      x <- x[,sample.int(ncol(x), use.subset, replace=FALSE)]
  }
  if(is(x,"SingleCellExperiment")) x <- counts(x)
  
  rs <- Matrix::rowSums(x)
  xo <- xo[which(rs>=minCount),]
  x <- x[which(rs>=minCount),]
  
  x <- norm.fn(x)

  fc <- .clusterFeaturesStep(x, k=k, dims.use=dims.use, use.mbk=use.mbk,
                             num_init=num_init, ...)

  if(twoPass & use.subset<ncol(xo)){
    message("Second iteration...")
    x <- t(normalizeCounts(scuttle::sumCountsAcrossFeatures(xo, fc)))
    cellclust <- kmeans(scale(x), min(1000,ceiling(use.subset/2)), iter.max=100, 
                        nstart=num_init)$cluster
    x <- sumCountsAcrossCells(xo, cellclust)
    if(is(x,"SummarizedExperiment")) x <- assay(x)
    x <- norm.fn(x)
    fc <- .clusterFeaturesStep(x, k=k, dims.use=seq_len(max(dims.use)),
                               use.mbk=use.mbk, num_init=num_init, ...)
  }
  fg <- setNames(fc, row.names(xo))
  x <- scuttle::sumCountsAcrossFeatures(xo, fc)
  row.names(x) <- paste0("feat",seq_len(nrow(x)))
  if(is(xo,"SingleCellExperiment")){
    x <- SingleCellExperiment(list(counts=x), colData=colData(xo),
                              reducedDims=reducedDims(xo))
    metadata(x)$featureGroups <- fg
  }
  x
}

# used by aggregateFeatures
.clusterFeaturesStep <- function(x, k, dims.use=2L:12L, use.mbk=NULL, 
                                 num_init=3L, ...){
  pca <- runPCA(x, BSPARAM=IrlbaParam(), center=FALSE,
                rank=max(dims.use))$x[,dims.use]
  if(is.null(use.mbk)) use.mbk <- nrow(x) > 30000
  if(use.mbk && requireNamespace("mbkmeans", quietly=TRUE)){
    fc <- mbkmeans::mbkmeans(t(pca), k, num_init=num_init, ...)$Clusters
  }else{
    fc <- kmeans(pca, k, nstart=num_init, iter.max=100)$cluster
  }
  fc
}
