#' rankTrans
#'
#' A dense rank transformation that preserves 0 and rank step size in the presence of many
#' ties.
#'
#' @param x A matrix, with samples/cells as columns and genes/features as rows.
#'
#' @return rank-transformed x.
#' @export
rankTrans <- function(x){
  y <- apply(x,2,ties.method="dense",FUN=data.table::frank)-1
  y <- matrix(as.integer(t(max(y)*t(y)/colMaxs(y))),nrow=nrow(y))
  dimnames(y) <- dimnames(x)
  y
}

#' plotROCs
#' 
#' Plot ROC curves for given scores.
#'
#' @param scores A data.frame with the different types of scores as columns.
#' @param truth A vector of the true class corresponding to each row of `scores`
#'
#' @return a ggplot
#' @export
plotROCs <- function(scores, truth){
  truth <- as.integer(as.factor(truth))-1
  library(ggplot2)
  scores <- as.data.frame(scores)
  roclist <- lapply(scores, FUN=function(x){
    labels <- truth[order(x, decreasing=TRUE)]
    data.frame(FPR=cumsum(!labels)/sum(!labels),
               TPR=cumsum(labels)/sum(labels))
  })
  d <- data.frame( method=factor( rep(names(roclist), sapply(roclist, FUN=function(x) length(x$TPR))), levels=names(roclist)),
                   FPR=unlist(lapply(roclist,FUN=function(x) x$FPR)),
                   TPR=unlist(lapply(roclist,FUN=function(x) x$TPR)) )
  ggplot(d, aes(FPR, TPR, colour=method)) + geom_line(size=1.2)
}

# get random cross-cluster pairs of cells from a cluster assignment vector
.getCellPairs <- function(clusters, n=1000){
  cli <- split(1:length(clusters), clusters)
  ca <- expand.grid(unique(clusters), unique(clusters))
  ca <- ca[which(ca[,1]!=ca[,2]),]
  ca <- do.call(rbind, lapply(1:nrow(ca), n=ceiling(n/nrow(ca)), FUN=function(i,n){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca[!duplicated(ca),]
}  

# creates within-cluster meta-cells from a count matrix
.getMetaCells <- function(x, clusters, n.meta.cells=20, meta.cell.size=20){
  cli <- split(1:length(clusters), clusters)
  meta <- unlist(lapply(cli, FUN=function(x){
    lapply(1:n.meta.cells, FUN=function(y){
      sample(x,min(ceiling(0.6*length(x)),meta.cell.size),replace=FALSE)
    })
  }), recursive=FALSE)
  meta <- sapply(meta, FUN=function(y){ Matrix::rowMeans(x[,y,drop=FALSE]) })
  colnames(meta) <- paste0("metacell.",1:ncol(meta))
  meta
}

#' scdsWrapper
#' 
#' A wrapper around scds's hybrid score, used for comparison (requires the `scds` package 
#' to be installed)
#'
#' @param sce An object of class `SingleCellExperiment`
#'
#' @return The updated `sce` object, with the colData field `hybrid_score`
#' @export
scdsWrapper <- function(sce){
  require(scds)
  sce <- scds::cxds(sce)
  sce <- scds::bcds(sce)
  scds::cxds_bcds_hybrid(sce)
}


#' dblFinderWrapper
#'
#' A wrapper around \url{https://github.com/chris-mcginnis-ucsf/DoubletFinder}[DoubletFinder],
#' used for comparison. Requires the `DoubletFinder` and `Seurat` version 3 packages to be
#'  installed.
#' 
#' @param sce An object of class `SingleCellExperiment`
#' @param doublet.formation.rate The expected doublet rate, default 0.025
#' @param dims The dimensions ot use, default 1:10
#'
#' @return The updated `sce` object, with the colData field `DF.score`
#' @export
dblFinderWrapper <- function(sce, doublet.formation.rate=0.025, dims=1:10){
  se <- CreateSeuratObject(counts(sce))    
  se <- NormalizeData(se, verbose=FALSE)
  se <- ScaleData(se, verbose=FALSE)
  se <- FindVariableFeatures(se, verbose=FALSE)
  se <- RunPCA(se, pcs.print=0, verbose=FALSE)
  se <- RunTSNE(se, dims.use=dims, verbose=FALSE)
  
  ## pK Identification
  sweep.res <- paramSweep_v3(se, PCs=dims)
  
  bcmvn <- find.pK(summarizeSweep(sweep.res, GT=FALSE))
  pK <- as.numeric(as.character(bcmvn[order(bcmvn$BCmetric,decreasing=T)[1],2]))
  message("pK=",pK)
  
  ## Homotypic Doublet Proportion Estimate
  clusters <- scran::quickCluster(sce)
  homotypic.prop <- modelHomotypic(clusters)
  nExp_poi <- round(doublet.formation.rate*ncol(sce))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder
  se <- doubletFinder_v3(se, PCs=dims, pK=pK, nExp=nExp_poi.adj)
  n <- ncol(se@meta.data)
  df <- se@meta.data[,(n-1):n]
  sce$DF.score <- df[,1]
  sce$DF.class <- df[,2]
  sce
}
