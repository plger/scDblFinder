#' Smooth Doublet Scores
#'
#' This function applies a distance-weighted kNN smoothing to single-cell 
#' doublet probabilities, amplifying high-probability doublet clusters while 
#' suppressing isolated noisy predictions.
#'
#' @param x A \code{SingleCellExperiment} object (containing the 
#'   'scDblFinder.score' colData column), or a numeric vector of doublet 
#'   probabilities.
#' @param knn Optional k nearest neighbors. A list containing \code{index} and 
#'   \code{distance} matrices, typically the output of 
#'   \code{\link[BiocNeighbors]{findKNN}}. If \code{NULL}, the kNN graph will 
#'   be computed.
#' @param coords Either a character scalar indicating the name of the reduced 
#'   dimension to use for kNN computation, or a matrix of such reduced 
#'   dimensions, with cells as rows and dimensions as columns. Ignored if 
#'   \code{knn} is given. Default "PCA".
#' @param k Integer. The number of nearest neighbors to compute if \code{knn} 
#'   is \code{NULL}. Defaults to 30.
#' @param alpha Numeric [0, 1]. The blending factor between the cell's own 
#'   score (0) and the neighborhood consensus (1). Defaults to 0.5.
#' @param gamma Numeric >= 1.0. The non-linear amplification exponent. Values 
#'   above 1 exponentially amplify clusters of high-probability cells. Defaults
#'    to 2.0.
#' @param weights Logical. If \code{TRUE} (default), uses adaptive Gaussian 
#'   kernel weighting based on continuous distance. If \code{FALSE}, 
#'   applies uniform weights across all neighbors.
#' @param decayByDistance Whether (and how) to reduce the importance of the
#'   neighborhood for isolated cells (based on the first neighbor distance),
#'   forcing them to rely primarily on their own raw score. Default 'linear'.
#' @param scoreColumn Character. The column name in \code{colData(x)} 
#'   containing the raw doublet scores. Defaults to \code{"scDblFinder.score"}.
#' @param outColumn Character. The column name to store the smoothed scores if 
#'   \code{x} is a \code{SingleCellExperiment}. Defaults to 
#'   \code{"smoothedDoubletScore"}.
#' @param ... Passed to \code{\link[BiocNeighbors]{findKNN}}.
#'
#' @return If \code{x} is a \code{SingleCellExperiment}, returns the object with 
#'   an added \code{colData} column containing the smoothed scores. If \code{x} 
#'   is a numeric vector, returns a numeric vector of smoothed scores.
#'
#' @importFrom BiocNeighbors findKNN
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods is
#' @export
smoothDoubletScores <- function(x,
                                knn = NULL,
                                coords = "PCA",
                                k = 20, 
                                alpha = 0.5,
                                gamma = 1,
                                weights = TRUE,
                                decayByDistance = c("linear", "exp", "none"),
                                scoreColumn = "scDblFinder.score",
                                outColumn = "smoothedDoubletScore",
                                ... ){
  
  decayByDistance <- match.arg(decayByDistance)
  if(is(x, "SingleCellExperiment")){
    stopifnot(scoreColumn %in% colnames(colData(x)))
    probs <- colData(x)[[scoreColumn]]
    if(is.null(knn)){
      if(is.null(dim(coords))){
        stopifnot(length(coords) == 1, is.character(coords),
                  coords %in% reducedDimNames(x))
        coords <- reducedDim(x, coords)
      }
    }
  } else if(is.numeric(x)){
    probs <- x
    if(is.null(knn) && (is.null(coords) || is.null(dim(coords))))
      stop("If 'x' is a numeric vector, provide either 'knn' or a reducedDim ",
           "matrix as 'coords'.")
  } else {
    stop("'x' must be a SingleCellExperiment or a numeric vector.")
  }
  stopifnot(all(probs>=0 & probs<=1))
  
  if(is.null(knn)){
    knn <- BiocNeighbors::findKNN(coords, k = k, ...)
  }
  
  probs <- pmax(1e-6, pmin(1 - 1e-6, probs))
  logits <- log(probs / (1 - probs))
  neighborLogits <- array(logits[knn$index], dim = dim(knn$index))
  
  # Compute neighborhood consensus weights
  if (weights) {
    # Adaptive continuous Gaussian distance weights
    minDist <- min(knn$distance[which(knn$distance > 0)]) / 10
    sigma <- apply(knn$distance, 1, median) + minDist
    weights <- exp(- (knn$distance^2) / (2 * (sigma^2)))
    weights <- sweep(weights, 1, rowSums(weights), "/")
  } else {
    # Flat uniform weights (1 / k for all neighbors)
    weights <- matrix(1/ncol(knn$index), nrow=nrow(knn$index),
                      ncol=ncol(knn$index))
  }
  
  # Neighborhood consensus
  if(gamma==1){
    weightedConsensus <- rowSums(weights * neighborLogits)
  }else{
    signedPowerLogits <- sign(neighborLogits) * (abs(neighborLogits)^gamma)
    weightedConsensus <- rowSums(weights * signedPowerLogits)
    weightedConsensus <- sign(weightedConsensus) * 
                                (abs(weightedConsensus)^(1 / gamma))
  }

  # Determine per-cell alpha (fixed vs distance-decayed)
  if(decayByDistance!="none"){
    # Alpha drops when the first neighbor is farther from the median first 
    # neighbor distance
    firstDist <- knn$distance[,1]
    medFirstDist <- median(firstDist) + 1e-8
    if(decayByDistance=="linear"){
      cellAlpha <- alpha*pmin(firstDist/medFirstDist, 1)
    }else{
      cellAlpha <- alpha * exp(- (firstDist / medFirstDist))
    }
  }else{
    cellAlpha <- alpha
  }
  
  # Blend self with neighborhood and apply sigmoid to return to prob space
  finalLogits <- (1 - cellAlpha) * logits + cellAlpha * weightedConsensus
  smoothedProbs <- 1 / (1 + exp(-finalLogits))
  
  if(!is(x, "SingleCellExperiment")) return(smoothedProbs)
  
  colData(x)[[outColumn]] <- smoothedProbs
  x
}
