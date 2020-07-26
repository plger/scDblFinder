#' Recover intra-sample doublets
#'
#' Recover intra-sample doublets that are neighbors to known inter-sample doublets in a multiplexed experiment.
#'
#' @param x A log-expression matrix for all cells (including doublets) in columns and genes in rows.
#' If \code{transposed=TRUE}, this should be a matrix of low-dimensional coordinates where each row corresponds to a cell.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing 
#' (i) a log-expression matrix in the \code{\link{assays}} as specified by \code{assay.type},
#' or (ii) a matrix of reduced dimensions in the \code{\link{reducedDims}} as specified by \code{use.dimred}.
#' @param doublets A logical, integer or character vector specifying which cells in \code{x} are known (inter-sample) doublets.
#' @param samples A numeric vector containing the relative proportions of cells from each sample,
#' used to determine how many cells are to be considered as intra-sample doublets.
#' @param k Integer scalar specifying the number of nearest neighbors to use for computing the local doublet proportions.
#' @param transposed Logical scalar indicating whether \code{x} is transposed, i.e., cells in the rows.
#' @param subset.row A logical, integer or character vector specifying the genes to use for the neighbor search. 
#' Only used when \code{transposed=FALSE}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the nearest neighbor search.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization to use for the nearest neighbor search.
#' @param ... For the generic, additional arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment method, additional arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, additional arguments to pass to the SummarizedExperiment method.
#' @param assay.type A string specifying which assay values contain the log-expression matrix. 
#' @param use.dimred A string specifying whether existing values in \code{\link{reducedDims}(x)} should be used.
#'
#' @return
#' A \linkS4class{DataFrame} containing one row per cell and the following fields:
#' \itemize{
#' \item \code{proportion}, a numeric field containing the proportion of neighbors that are doublets.
#' \item \code{known}, a logical field indicating whether this cell is a known inter-sample doublet.
#' \item \code{predicted}, a logical field indicating whether this cell is a predicted intra-sample doublet.
#' }
#' The \code{\link{metadata}} contains \code{intra}, a numeric scalar containing the expected number of intra-sample doublets. 
#'
#' @details
#' In multiplexed single-cell experiments, we can detect doublets as libraries with labels for multiple samples.
#' However, this approach fails to identify doublets consisting of two cells with the same label.
#' Such cells may be problematic if they are still sufficiently abundant to drive formation of spurious clusters.
#'
#' This function identifies intra-sample doublets based on the similarity in expression profiles to known inter-sample doublets.
#' For each cell, we compute the proportion of the \code{k} neighbors that are known doublets.
#' Of the \dQuote{unmarked} cells that are not known doublets, those with top \eqn{X} largest proportions are considered to be intra-sample doublets.
#' We use \code{samples} to obtain a reasonable estimate for \eqn{X}, see the vignette for details.
#'
#' A larger value of \code{k} provides more stable estimates of the doublet proportion in each cell.
#' However, this comes at the cost of assuming that each cell actually has \code{k} neighboring cells of the same state.
#' For example, if a doublet cluster has fewer than \code{k} members,
#' its doublet proportions will be \dQuote{diluted} by inclusion of unmarked cells in the next-closest cluster.
#'
#' @author Aaron Lun
#' 
#' @seealso
#' \code{\link{doubletCells}} and \code{\link{doubletCluster}},
#' for alternative methods of doublet detection when no prior doublet information is available.
#'
#' \code{hashedDrops} from the \pkg{DropletUtils} package,
#' to identify doublets from cell hashing experiments.
#'
#' More detail on the mathematical background of this function is provided in the corresponding vignette at
#' \code{vignette("recoverDoublets", package="scDblFinder")}.
#'
#' @examples
#' # Mocking up an example.
#' set.seed(100)
#' ngenes <- 1000
#' mu1 <- 2^rnorm(ngenes, sd=2)
#' mu2 <- 2^rnorm(ngenes, sd=2)
#'
#' counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes) # Pure type 1
#' counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes) # Pure type 2
#' counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes) # Doublets (1 & 2)
#' all.counts <- cbind(counts.1, counts.2, counts.m)
#' lcounts <- scuttle::normalizeCounts(all.counts)
#' 
#' # Pretending that half of the doublets are known. Also pretending that 
#' # the experiment involved two samples of equal size.
#' known <- 200 + seq_len(10) 
#' out <- recoverDoublets(lcounts, doublets=known, k=10, samples=c(1, 1))
#' out
#'
#' @name recoverDoublets 
NULL

#' @importFrom Matrix t
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom utils head
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom scuttle .subset2index
#' @importFrom BiocParallel SerialParam
.doublet_recovery <- function(x, doublets, samples,
    k=50, transposed=FALSE, subset.row=NULL, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }

    is.doublet <- logical(nrow(x))
    is.doublet[.subset2index(doublets, x, byrow=TRUE)] <- TRUE

    fout <- findKNN(x, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    neighbors <- fout$index
    neighbors[] <- is.doublet[neighbors]
    P <- rowMeans(neighbors)

    expected.intra <- sum(samples^2)/sum(samples)^2
    intra.doublets <- sum(is.doublet) * expected.intra/(1 - expected.intra)

    predicted <- logical(nrow(x))
    o <- order(P[!is.doublet], decreasing=TRUE)
    predicted[!is.doublet][head(o, intra.doublets)] <- TRUE

    output <- DataFrame(proportion=P, known=is.doublet, predicted=predicted)
    metadata(output)$intra <- intra.doublets
    output
}

#' @export
#' @rdname recoverDoublets
#' @import methods
setGeneric("recoverDoublets", function(x, ...) standardGeneric("recoverDoublets"))

#' @export
#' @rdname recoverDoublets
setMethod("recoverDoublets", "ANY", .doublet_recovery)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname recoverDoublets
setMethod("recoverDoublets", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .doublet_recovery(assay(x, assay.type), ...)
})

#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @rdname recoverDoublets
setMethod("recoverDoublets", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .doublet_recovery(reducedDim(x, use.dimred), transposed=TRUE, ...)
    } else {
        callNextMethod(x=x, ...)
    }
})

