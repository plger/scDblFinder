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
#' @param called.class A numeric vector (with names corresponding to the columns
#' of `scores`) indicating, for each method, the number of cases called as true
#' (i.e. the threshold decided). Those will be plotted as points.
#' @param nbT Logical; whether to add a dot for each score at the number of 
#' true positives (default TRUE).
#' @param dot.size The size of the dots.
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
plotROCs <- function(scores, truth, called.class=NULL, nbT=TRUE, dot.size=5){
  truth <- as.integer(as.factor(truth))-1
  scores <- as.data.frame(scores)
  roclist <- lapply(scores, FUN=function(x){
    labels <- truth[order(x, decreasing=TRUE)]
    data.frame(FPR=cumsum(!labels)/sum(!labels),
               TPR=cumsum(labels)/sum(labels))
  })
  methods <- factor( rep(names(roclist), 
                         vapply(roclist, FUN.VALUE=integer(1), 
                                FUN=function(x) length(x$TPR))
                         ),
                     levels=names(roclist) )
  d <- data.frame( method=methods,
                   FPR=unlist(lapply(roclist,FUN=function(x) x$FPR)),
                   TPR=unlist(lapply(roclist,FUN=function(x) x$TPR)) )
  p <- ggplot(d, aes(FPR, TPR, colour=method)) + geom_line(size=1.2)
  if(nbT){
      nbt <- sum(truth)
      d2 <- data.frame( method=names(roclist),
                        FPR=sapply(roclist,FUN=function(x) x$FPR[nbt]),
                        TPR=unlist(lapply(roclist,FUN=function(x) x$TPR[nbt])) )
      p <- p + geom_point(data=d2, alpha=0.5, size=dot.size)
  }
  if(!is.null(called.class) &&
     length(mm <- intersect(names(called.class), names(roclist)))>0){
      d2 <- data.frame( method=mm,
                        FPR=mapply(r=roclist[mm], cc=called.class[mm], FUN=function(r,cc) r$FPR[cc]),
                        TPR=mapply(r=roclist[mm], cc=called.class[mm], FUN=function(r,cc) r$TPR[cc]) )
      p <- p + geom_point(data=d2, size=dot.size, shape=8, show.legend = FALSE)
  }
  
  p
}

# get random cross-cluster pairs of cells from a cluster assignment vector
.getCellPairs <- function(clusters, n=1000){
  cli <- split(seq_along(clusters), clusters)
  ca <- expand.grid(unique(clusters), unique(clusters))
  ca <- t(apply(ca, 1, FUN=function(x) sort(x)))
  ca <- !duplicated(ca)
  ca <- ca[which(ca[,1]!=ca[,2]),]
  n <- ceiling(n/nrow(ca))
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i,n){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  oc <- paste(rep(ca[,1], each=n), rep(ca[,1], each=n), sep="+")
  ca <- data.frame(ca, orig.clusters=as.factor(oc))
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
  meta <- vapply(meta, FUN.VALUE=double(nrow(x)), 
                 FUN=function(y){ Matrix::rowMeans(x[,y,drop=FALSE]) })
  colnames(meta) <- paste0("metacell.",seq_len(ncol(meta)))
  meta
}

