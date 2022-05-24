#' Compute the density of simulated doublets
#' 
#' Identify potential doublet cells based on the local density of simulated doublet expression profiles.
#' This replaces the older \code{doubletCells} function from the \pkg{scran} package. 
#' 
#' @param x A numeric matrix-like object of count values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param size.factors.norm A numeric vector of size factors for normalization of \code{x} prior to PCA and distance calculations.
#' If \code{NULL}, defaults to size factors derived from the library sizes of \code{x}.
#' 
#' For the SingleCellExperiment method, the default values are taken from \code{\link{sizeFactors}(x)}, if they are available.
#' @param size.factors.content A numeric vector of size factors for RNA content normalization of \code{x} prior to simulating doublets.
#' This is orthogonal to the values in \code{size.factors.norm}, see Details.
#' @param k An integer scalar specifying the number of nearest neighbours to use to determine the bandwidth for density calculations.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param niters An integer scalar specifying how many simulated doublets should be generated.
#' @param block An integer scalar controlling the rate of doublet generation, to keep memory usage low.
#' @param dims An integer scalar specifying the number of components to retain after the PCA.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link{findNeighbors}}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the neighbour searches should be parallelized.
#' @param ... For the generic, additional arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values contain the count matrix. 
#' 
#' @return 
#' A numeric vector of doublet scores for each cell in \code{x}.
#' 
#' @details
#' This function simulates doublets by adding the count vectors for two randomly chosen cells in \code{x}.
#' For each original cell, we compute the density of neighboring simulated doublets and compare it to the density of neighboring original cells.
#' Genuine doublets should have a high density of simulated doublets relative to the density of its neighbourhood.
#' Thus, the doublet score for each cell is defined as the ratio of densities of simulated doublets to the density of the original cells.
#' 
#' Densities are calculated in low-dimensional space after a PCA on the log-normalized expression matrix of \code{x}.
#' Simulated doublets are projected into the low-dimensional space using the rotation vectors computed from the original cells.
#' For each cell, the density of simulated doublets is computed for a hypersphere with radius set to the median distance to the \code{k} nearest neighbour.
#' This is normalized by \code{niters}, \code{k} and the total number of cells in \code{x} to yield the final score.
#' 
#' The two size factor arguments have different roles:
#' \itemize{
#' \item \code{size.factors.norm} contains the size factors to be used for normalization prior to PCA and distance calculations.
#' This defaults to the values returned by \code{\link{librarySizeFactors}} but can be explicitly set to ensure that the low-dimensional space is consistent with that in the rest of the analysis.
#' \item \code{size.factors.content} is much more important, and represents the size factors that preserve RNA content differences.
#' This is usually computed from spike-in RNA and ensures that the simulated doublets have the correct ratio of contributions from the original cells.
#' }
#' It is possible to set both of these arguments as they are orthogonal to each other.
#' Setting \code{size.factors.content} will not affect the calculation of log-normalized expression values from \code{x}.
#' Conversely, setting \code{size.factors.norm} will not affect the ratio in which cells are added together when simulating doublets.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' # Mocking up an example.
#' set.seed(100)
#' ngenes <- 1000
#' mu1 <- 2^rnorm(ngenes)
#' mu2 <- 2^rnorm(ngenes)
#' mu3 <- 2^rnorm(ngenes)
#' mu4 <- 2^rnorm(ngenes)
#' 
#' counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes) # Pure type 1
#' counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes) # Pure type 2
#' counts.3 <- matrix(rpois(ngenes*100, mu3), nrow=ngenes) # Pure type 3
#' counts.4 <- matrix(rpois(ngenes*100, mu4), nrow=ngenes) # Pure type 4
#' counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes) # Doublets (1 & 2)
#' 
#' counts <- cbind(counts.1, counts.2, counts.3, counts.4, counts.m)
#' clusters <- rep(1:5, c(rep(100, 4), ncol(counts.m)))
#' 
#' # Find potential doublets.
#' scores <- computeDoubletDensity(counts)
#' boxplot(split(log10(scores), clusters))
#' 
#' @references
#' Lun ATL (2018).
#' Detecting doublet cells with \emph{scran}.
#' \url{https://ltla.github.io/SingleCellThoughts/software/doublet_detection/bycell.html}
#'
#' @seealso
#' \code{\link{findDoubletClusters}}, to detect doublet clusters.
#'
#' \code{\link{scDblFinder}}, which uses a hybrid approach involving simulation and overclustering.
#'
#' More detail on the mathematical background of this function is provided in the corresponding vignette at
#' \code{vignette("computeDoubletDensity", package="scDblFinder")}.
#'
#' @name computeDoubletDensity
NULL

#' @importFrom scuttle librarySizeFactors normalizeCounts .bpNotSharedOrUp
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom BiocParallel SerialParam bpmapply bpstart bpstop
#' @importFrom Matrix rowMeans
#' @importFrom stats median
#' @importFrom BiocNeighbors findKNN findNeighbors queryNeighbors queryKNN buildIndex
#' @importFrom BiocSingular runPCA bsparam
#' @importFrom methods is
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.doublet_cells <- function(x, size.factors.norm=NULL, size.factors.content=NULL,
    k=50, subset.row=NULL, niters=max(10000, ncol(x)), block=10000, dims=25, 
    BNPARAM=KmknnParam(), BSPARAM=bsparam(), BPPARAM=SerialParam())
{
    # Setting up the parallelization.
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    if (.bpNotSharedOrUp(BPPARAM)){ 
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (is.null(size.factors.norm)) {
        size.factors.norm <- librarySizeFactors(x, BPPARAM=BPPARAM)
    }
    if(!all(size.factors.norm>0))
        stop("Some size.factors are not positive. This typically happens ",
             "because some cells have no reads in the features specified by ",
             "`subset.row` -- these should be filtered out.")

    # Manually controlling the size factor centering here to ensure the final counts are on the same scale.
    size.factors.norm <- size.factors.norm/mean(size.factors.norm)
    if (!is.null(size.factors.content)) {
        x <- normalizeCounts(x, size.factors.content, log=FALSE, center_size_factors=FALSE)
        size.factors.norm <- size.factors.norm/size.factors.content
    }
    y <- normalizeCounts(x, size.factors.norm, center_size_factors=FALSE)

    # Running the PCA.
    pc.out <- runPCA(t(y), center=TRUE, BSPARAM=BSPARAM, rank=dims, BPPARAM=BPPARAM)
    pcs <- pc.out$x
    sim.pcs <- .spawn_doublet_pcs(x, size.factors.norm, V=pc.out$rotation, centers=rowMeans(y), niters=niters, block=block)

    # Computing densities, using a distance computed from the kth nearest neighbor.
    pre.pcs <- buildIndex(pcs, BNPARAM=BNPARAM)
    self.dist <- findKNN(BNINDEX=pre.pcs, k=k, BPPARAM=BPPARAM, last=1, get.index=FALSE, warn.ties=FALSE)$distance

    sim.n <- queryNeighbors(sim.pcs, query=pcs, 
        threshold=self.dist * 1.00000001, # bump it up to avoid issues with numerical precision during tests.
        BNPARAM=BNPARAM, BPPARAM=BPPARAM, 
        get.distance=FALSE, get.index=FALSE)

    sim.prop <- sim.n/niters
    sim.prop/(k/ncol(x))
}

#' @importFrom Matrix crossprod
#' @importFrom scuttle normalizeCounts
#' @importFrom DelayedArray sweep
.spawn_doublet_pcs <- function(x, size.factors, V, centers, niters=10000L, block=10000L) {
    collected <- list()
    counter <- 1L
    current <- 0L
    mean.correction <- colSums(centers * V)

    while (current < niters) {
        to.make <- min(block, niters - current)
        left <- sample(ncol(x), to.make, replace=TRUE)
        right <- sample(ncol(x), to.make, replace=TRUE)
        sim.x <- x[,left,drop=FALSE] + x[,right,drop=FALSE]

        # Do not center, otherwise the simulated doublets will always have higher normalized counts
        # than actual doublets (as the latter will have been normalized to the level of singlets).
        sim.sf <- size.factors[left] + size.factors[right]
        sim.y <- normalizeCounts(sim.x, sim.sf, center_size_factors=FALSE)

        # Projecting onto the PC space of the original data.
        sim.pcs <- crossprod(sim.y, V)
        sim.pcs <- sweep(sim.pcs, 2L, mean.correction, FUN="-", check.margin=FALSE)
        collected[[counter]] <- sim.pcs
        counter <- counter + 1L
        current <- current + block
    }

    do.call(rbind, collected)
}

##############################
# S4 method definitions here #
##############################

#' @export
#' @rdname computeDoubletDensity
setGeneric("computeDoubletDensity", function(x, ...) standardGeneric("computeDoubletDensity"))

#' @export
#' @rdname computeDoubletDensity
setMethod("computeDoubletDensity", "ANY", .doublet_cells)

#' @export
#' @rdname computeDoubletDensity
#' @importFrom SummarizedExperiment assay
setMethod("computeDoubletDensity", "SummarizedExperiment", function(x, ..., assay.type="counts")
{
    .doublet_cells(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname computeDoubletDensity
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
setMethod("computeDoubletDensity", "SingleCellExperiment", function(x, size.factors.norm=sizeFactors(x), ...) {
    callNextMethod(x=x, size.factors.norm=size.factors.norm, ...)
})
