#' Detect doublet clusters
#'
#' Identify potential clusters of doublet cells based on whether they have intermediate expression profiles,
#' i.e., their profiles lie between two other \dQuote{source} clusters.
#'
#' @param x A numeric matrix-like object of count values,
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param clusters A vector of length equal to \code{ncol(x)}, containing cluster identities for all cells.
#' If \code{x} is a SingleCellExperiment, this is taken from \code{\link{colLabels}(x)} by default.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param threshold A numeric scalar specifying the FDR threshold with which to identify significant genes.
#' @param ... For the generic, additional arguments to pass to specific methods.
#'
#' For the ANY method, additional arguments to pass to \code{\link{findMarkers}}.
#'
#' For the SummarizedExperiment method, additional arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, additional arguments to pass to the SummarizedExperiment method.
#' @param assay.type A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.
#' @param get.all.pairs Logical scalar indicating whether statistics for all possible source pairings should be returned.
#'
#' @return
#' A \linkS4class{DataFrame} containing one row per query cluster with the following fields:
#' \describe{
#' \item{\code{source1}:}{String specifying the identity of the first source cluster.}
#' \item{\code{source2}:}{String specifying the identity of the second source cluster.}
#' \item{\code{num.de}:}{Integer, number of genes that are significantly non-intermediate
#' in the query cluster compared to the two putative source clusters.}
#' \item{\code{median.de}:}{Integer, median number of genes that are significantly non-intermediate
#' in the query cluster across all possible source cluster pairings.}
#' \item{\code{best}:}{String specifying the identify of the top gene with the lowest p-value
#' against the doublet hypothesis for this combination of query and source clusters.}
#' \item{\code{p.value}:}{Numeric, containing the adjusted p-value for the \code{best} gene.}
#' \item{\code{lib.size1}:}{Numeric, ratio of the median library sizes for the first source cluster to the query cluster.}
#' \item{\code{lib.size2}:}{Numeric, ratio of the median library sizes for the second source cluster to the query cluster.}
#' \item{\code{prop}:}{Numeric, proportion of cells in the query cluster.}
#' \item{\code{all.pairs}:}{A \linkS4class{SimpleList} object containing the above statistics
#' for every pair of potential source clusters, if \code{get.all.pairs=TRUE}.}
#' }
#' Each row is named according to its query cluster.
#'
#' @details
#' This function detects clusters of doublet cells in a manner similar to the method used by Bach et al. (2017).
#' For each \dQuote{query} cluster, we examine all possible pairs of \dQuote{source} clusters,
#' hypothesizing that the query consists of doublets formed from the two sources.
#' If so, gene expression in the query cluster should be strictly intermediate
#' between the two sources after library size normalization.
#'
#' We apply pairwise t-tests to the normalized log-expression profiles to reject this null hypothesis.
#' This is done by identifying genes that are consistently up- or down-regulated in the query compared to \emph{both} sources.
#' We count the number of genes that reject the null hypothesis at the specified FDR \code{threshold}.
#' For each query cluster, the most likely pair of source clusters is that which minimizes the number of significant genes.
#'
#' Potential doublet clusters are identified using the following characteristics, in order of importance:
#' \itemize{
#' \item Low number of significant genes (i.e., \code{num.de}).
#' Ideally, \code{median.de} is also high to indicate that the absence of strong DE is not due to a lack of power.
#' \item A reasonable proportion of cells in the cluster, i.e., \code{prop}.
#' This requires some expectation of the doublet rate in the experimental protocol.
#' \item Library sizes of the source clusters that are below that of the query cluster, i.e., \code{lib.size*} values below unity.
#' This assumes that the doublet cluster will contain more RNA and have more counts than either of the two source clusters.
#' }
#'
#' For each query cluster, the function will only report the pair of source clusters with the lowest \code{num.de}.
#' Setting \code{get.all.pairs=TRUE} will retrieve statistics for all pairs of potential source clusters.
#' This can be helpful for diagnostics to identify relationships between specific clusters.
#'
#' The reported \code{p.value} is of little use in a statistical sense, and is only provided for inspection.
#' Technically, it could be treated as the Simes combined p-value against the doublet hypothesis for the query cluster.
#' However, this does not account for the multiple testing across all pairs of clusters for each chosen cluster,
#' especially as we are chosing the pair that is most concordant with the doublet null hypothesis.
#'
#' We use library size normalization (via \code{\link{librarySizeFactors}}) even if existing size factors are present.
#' This is because intermediate expression of the doublet cluster is not guaranteed for arbitrary size factors.
#' For example, expression in the doublet cluster will be higher than that in the source clusters if normalization was performed with spike-in size factors.
#'
#' @author
#' Aaron Lun
#'
#' @references
#' Bach K, Pensa S, Grzelak M, Hadfield J, Adams DJ, Marioni JC and Khaled WT (2017).
#' Differentiation dynamics of mammary epithelial cells revealed by single-cell RNA sequencing.
#' \emph{Nat Commun.} 8, 1:2128.
#'
#' @seealso
#' \code{\link{findMarkers}}, to detect DE genes between clusters.
#'
#' @examples
#' # Mocking up an example.
#' library(SingleCellExperiment)
#' sce <- mockDoubletSCE(c(200,300,200))
#'
#' # Compute doublet-ness of each cluster:
#' dbl <- findDoubletClusters(counts(sce), sce$cluster)
#' dbl
#'
#' # Narrow this down to clusters with very low 'N':
#' library(scuttle)
#' isOutlier(dbl$num.de, log=TRUE, type="lower")
#'
#' # Get help from "lib.size" below 1.
#' dbl$lib.size1 < 1 & dbl$lib.size2 < 1
#'
#' @name findDoubletClusters
NULL

#' @importFrom scuttle librarySizeFactors logNormCounts
#' @importFrom scran findMarkers .logBH
#' @importFrom BiocGenerics "sizeFactors<-" sizeFactors
#' @importFrom stats p.adjust median
#' @importFrom methods as
#' @importClassesFrom S4Vectors SimpleList
.doublet_cluster <- function(x, clusters, subset.row=NULL, threshold=0.05, get.all.pairs=FALSE, ...) {
    if (length(unique(clusters)) < 3L) {
        stop("need at least three clusters to detect doublet clusters")
    }

    # Computing normalized counts using the library size (looking for compositional differences!)
    sce <- SingleCellExperiment(list(counts=x))
    sizeFactors(sce) <- librarySizeFactors(x, subset_row=subset.row)
    sce <- logNormCounts(sce)

    degs <- findMarkers(sce, clusters, subset.row=subset.row, full.stats=TRUE, ...)
    med.lib.size <- vapply(split(sizeFactors(sce), clusters), FUN=median, FUN.VALUE=0)
    n.cluster <- table(clusters)/length(clusters)

    # Setting up the output.
    all.clusters <- names(degs)
    collected.top <- collected.all <- vector("list", length(all.clusters))
    names(collected.top) <- names(collected.all) <- all.clusters

    # Running through all pairs of clusters and testing against the third cluster.
    for (ref in all.clusters) {
        ref.stats <- degs[[ref]]
        remnants <- setdiff(all.clusters, ref)

        num <- length(remnants) * (length(remnants) - 1L)/2L
        all.N <- med.N <- all.gene <- all.parent1 <- all.parent2 <- integer(num)
        all.p <- numeric(num)
        idx <- 1L

        for (i1 in seq_along(remnants)) {
            stats1 <- ref.stats[[paste0("stats.", remnants[i1])]]
            for (i2 in seq_len(i1-1L)) {
                stats2 <- ref.stats[[paste0("stats.", remnants[i2])]]

                # Obtaining the IUT and setting opposing log-fold changes to 1.
                max.log.p <- pmax(stats1$log.p.value, stats2$log.p.value)
                max.log.p[sign(stats1$logFC) != sign(stats2$logFC)] <- 0

                # Correcting across genes. We use [1] to get NA when there are
                # no genes, which avoids an nrow() mismatch in DataFrame().
                log.adj.p <- .logBH(max.log.p)
                best.gene <- which.min(max.log.p)[1]

                all.N[idx] <- sum(log.adj.p <= log(threshold), na.rm=TRUE)
                all.gene[idx] <- best.gene
                all.p[idx] <- exp(log.adj.p[best.gene])
                all.parent1[idx] <- i1
                all.parent2[idx] <- i2
                idx <- idx + 1L
            }
        }

        # Formatting the output.
        parent1 <- remnants[all.parent1]
        parent2 <- remnants[all.parent2]

        stats <- DataFrame(source1=parent1, source2=parent2,
            num.de=all.N,
            median.de=rep(0, length(all.N)), # placeholder, see below.
            best=rownames(ref.stats)[all.gene],
            p.value=all.p,
            lib.size1=unname(med.lib.size[parent1]/med.lib.size[ref]),
            lib.size2=unname(med.lib.size[parent2]/med.lib.size[ref]))

        o <- order(all.N, -all.p)
        top <- cbind(stats[o[1],], prop=n.cluster[[ref]])
        med.de <- median(all.N)
        top$median.de <- med.de
        rownames(top) <- ref
        collected.top[[ref]] <- top

        if (get.all.pairs) {
            stats$median.de <- NULL
            collected.all[[ref]] <- stats[o,]
        }
    }

    # Returning the DataFrame of compiled results.
    out <- do.call(rbind, collected.top)
    if (get.all.pairs) {
        out$all.pairs <- as(collected.all, "SimpleList")
    }
    out[order(out$num.de),]
}

##############################
# S4 method definitions here #
##############################

#' @export
#' @rdname findDoubletClusters
setGeneric("findDoubletClusters", function(x, ...) standardGeneric("findDoubletClusters"))

#' @export
#' @rdname findDoubletClusters
setMethod("findDoubletClusters", "ANY", .doublet_cluster)

#' @export
#' @rdname findDoubletClusters
#' @importFrom SummarizedExperiment assay
setMethod("findDoubletClusters", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .doublet_cluster(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname findDoubletClusters
#' @importFrom SingleCellExperiment colLabels
setMethod("findDoubletClusters", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ...) {
    callNextMethod(x=x, clusters=clusters, ...)
})
