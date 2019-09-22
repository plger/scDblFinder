#' rankTrans
#'
#' A dense rank transformation that preserves 0 and rank step size in the 
#' presence of many ties.
#'
#' @param x A matrix, with samples/cells as columns and genes/features as rows.
#'
#' @return rank-transformed x.
#' 
#' @examples
#' m <- t(sapply( seq(from=0, to=5, length.out=30), 
#'                FUN=function(x) rpois(30,x) ) )
#' m2 <- rankTrans(m)
#' 
#' @importFrom data.table frank
#' @importFrom matrixStats colMaxs
#' @export
rankTrans <- function(x){
  y <- apply(x,2,ties.method="dense",FUN=frank)-1
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
#' @import ggplot2
#' 
#' @examples
#' myscores <- list( test=1:10 )
#' truth <- sample(c(TRUE,FALSE), 10, TRUE)
#' plotROCs(  myscores, truth )
#' 
#' @export
plotROCs <- function(scores, truth){
  truth <- as.integer(as.factor(truth))-1
  scores <- as.data.frame(scores)
  roclist <- lapply(scores, FUN=function(x){
    labels <- truth[order(x, decreasing=TRUE)]
    data.frame(FPR=cumsum(!labels)/sum(!labels),
               TPR=cumsum(labels)/sum(labels))
  })
  methods <- factor( rep(names(roclist), 
                         sapply(roclist, FUN=function(x) length(x$TPR))),
                     levels=names(roclist) )
  d <- data.frame( method=methods,
                   FPR=unlist(lapply(roclist,FUN=function(x) x$FPR)),
                   TPR=unlist(lapply(roclist,FUN=function(x) x$TPR)) )
  ggplot(d, aes(FPR, TPR, colour=method)) + geom_line(size=1.2)
}

# get random cross-cluster pairs of cells from a cluster assignment vector
.getCellPairs <- function(clusters, n=1000){
  cli <- split(seq_along(clusters), clusters)
  ca <- expand.grid(unique(clusters), unique(clusters))
  ca <- ca[which(ca[,1]!=ca[,2]),]
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), 
                               n=ceiling(n/nrow(ca)), 
                               FUN=function(i,n){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca[!duplicated(ca),]
}  

# creates within-cluster meta-cells from a count matrix
.getMetaCells <- function(x, clusters, n.meta.cells=20, meta.cell.size=20){
  cli <- split(seq_along(clusters), clusters)
  meta <- unlist(lapply(cli, FUN=function(x){
    lapply(seq_len(n.meta.cells), FUN=function(y){
      sample(x,min(ceiling(0.6*length(x)),meta.cell.size),replace=FALSE)
    })
  }), recursive=FALSE)
  meta <- sapply(meta, FUN=function(y){ Matrix::rowMeans(x[,y,drop=FALSE]) })
  colnames(meta) <- paste0("metacell.",seq_len(ncol(meta)))
  meta
}

