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
#'   \code{\link[BiocNeighbors](findKNN)}. If \code{NULL}, the kNN graph will 
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
#' @param scoreColumn Character. The column name in \code{colData(x)} 
#'   containing the raw doublet scores. Defaults to \code{"scDblFinder.score"}.
#' @param outColumn Character. The column name to store the smoothed scores if 
#'   \code{x} is a \code{SingleCellExperiment}. Defaults to 
#'   \code{"smoothedDoubletScore"}.
#'
#' @return If \code{x} is a \code{SingleCellExperiment}, returns the object with 
#'   an added \code{colData} column containing the smoothed scores. If \code{x} 
#'   is a numeric vector, returns a numeric vector of smoothed scores.
#'
#' @importFrom BiocNeighbors findKNN
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods is
#' @export
smoothDoubletScores <- function(x, knn = NULL, coords = "PCA", k=30, 
                                alpha = 0.5, gamma = 2.0, 
                                scoreColumn = "scDblFinder.score",
                                outColumn = "smoothedDoubletScore"){
    
    if(is(x, "SingleCellExperiment")){
      stopifnot(scoreColumn %in% colnames(colData(x)))
      probs <- colData(x)[[scoreColumn]]
      if(is.null(knn)){
        if(is.null(dim(coords))){
          stopifnot(length(coords)==1 && is.character(coords) &&
                      coords %in% names(reducedDim(x)))
          coords <- reducedDim(x, coords)
        }
      }
    }else if(is.numeric(x)){
      probs <- x
      if(is.null(knn) && (is.null(coords) || is.null(dim(coords))))
        stop("If 'x' is a numeric vector, provide either 'knn' or a reducedDim',
           'matrix as 'coords'.")
    }else{
      stop("'x' must be a SingleCellExperiment or a numeric vector.")
    }
    
    if(is.null(knn)){
        knn <- BiocNeighbors::findKNN(coords, k = k)
    }
    
    probs <- pmax(1e-6, pmin(1 - 1e-6, probs))
    logits <- log(probs / (1 - probs))
    neighborLogits <- array(logits[knn$index], dim = dim(knn$index))

    # compute adaptive continuous Gaussian distance weights
    minDist <- min(knn$distance[which(knn$distance>0)])/10
    sigma <- apply(knn$distance, 1, median) + minDist
    weights <- exp(- (knn$distance^2) / (2 * (sigma^2)))
    weights <- sweep(weights, 1, rowSums(weights), "/")
    
    # non-linear power scaling
    signedPowerLogits <- sign(neighborLogits) * (abs(neighborLogits)^gamma)
    
    # weighted neighborhood consensus
    weightedConsensus <- rowSums(weights * signedPowerLogits)
    
    # bring back to standard logit space
    nl <- sign(weightedConsensus) * (abs(weightedConsensus)^(1 / gamma))
    
    # blend self with neighborhood and apply sigmoid to return to prob space
    finalLogits <- (1 - alpha) * logits + alpha * nl
    smoothedProbs <- 1 / (1 + exp(-finalLogits))
    
    if(!is(x, "SingleCellExperiment")) return(smoothedProbs)
    
    colData(x)[[outColumn]] <- smoothedProbs
    x
}
