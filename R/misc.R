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
  ca <- expand.grid(seq_along(cli), seq_along(cli))
  ca <- ca[ca[,1]<ca[,2],]
  n <- ceiling(n/nrow(ca))
  oc <- paste( names(cli)[ca[,1]], names(cli)[ca[,2]], sep="+")
  ca <- do.call(rbind, lapply( seq_len(nrow(ca)), FUN=function(i){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca <- data.frame(ca, orig.clusters=rep(as.factor(oc), each=n))
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



.tablevec <- function(x){
  if(!is(x, "table")) x <- table(x)
  y <- as.numeric(x)
  names(y) <- names(x)
  y
}

.getMostLikelyOrigins <- function(knn, known.origins=NULL){
  origins <- t(vapply( seq_len(nrow(knn$orig)), FUN.VALUE=character(2),
                       FUN=function(i){
           if(all(is.na(knn$orig[i,]))) return(c(NA_character_,NA_character_))
           n <- .tablevec(knn$orig[i,])
           if(length(n)==1) return(c(names(n),FALSE))
           if(sum(n==max(n))==1) return(c(names(n)[which.max(n)],FALSE))
           x <- split(knn$distance[i,], knn$orig[i,])
           x <- sort(vapply(x,min,numeric(1)))
           if( (x[2]-x[1])/x[1] > 0.2 ) return(c(names(x)[1], FALSE))
           return(c(names(x)[1], TRUE))
                       }))
  origins <- as.data.frame(origins)
  colnames(origins) <- c("mostLikelyOrigin", "originAmbiguous")
  if(!is.null(known.origins)){
    w <- which(!is.na(known.origins))
    origins[w,1] <- known.origins[w]
    origins[w,2] <- "FALSE"
  }
  origins[,1] <- factor(origins[,1])
  origins[,2] <- as.logical(origins[,2])
  origins
}

.adjust.dbr.homotypy <- function(clusters, dbr=NULL){
  if(is.null(dbr)) dbr <- (0.01*length(clusters)/1000)
  # adjust expected dbr for expected homotypic doublets
  homotypic.prop <- sum((table(clusters)/length(clusters))^2)
  dbr*(1-homotypic.prop)
}

getExpectedDoublets <- function(x, dbr=NULL){
  if(is(x,"SingleCellExperiment")){
    clusters <- x$scDblFinder.clusters
  }else{
    clusters <- x
  }
  ncells <- length(clusters)
  dbr <- .adjust.dbr.homotypy(clusters, dbr)
  cs <- table(clusters)
  eg <- expand.grid(seq_along(unique(clusters)), seq_along(unique(clusters)))
  eg <- eg[eg[,1]<eg[,2],]
  nd <- ceiling(dbr*ncells)
  expected <- apply(eg,1, FUN=function(x){
    nd*(cs[[x[[1]]]]/ncells)*(cs[[x[[2]]]]/ncells)
  })
  names(expected) <- apply(eg,1,FUN=function(x){
    paste(names(cs)[x], collapse="+")
  })
  expected
}

.castorigins <- function(e, val=NULL){
  if(is.table(e) || is.null(dim(e))){
    e <- cbind(do.call(rbind, strsplit(names(e),"+",fixed=TRUE)),
               data.frame(val=as.numeric(e)))
  }
  if(is.null(val)) val <- rev(colnames(e))[1]
  names(n) <- n <- unique(as.character(as.matrix(e[,1:2])))
  sapply(n, FUN=function(x){
    sapply(n, FUN=function(y){
      if(x==y) return(NA)
      w <- which(e[,1] %in% c(x,y) & e[,2] %in% c(x,y))
      if(length(w)==0) return(0)
      sum(e[w,val])
    })
  })
}

.clusterTopG <- function(sce, clusters=NULL, nfeatures=1000){
  if(is.null(clusters))
    return(row.names(sce)[order(Matrix::rowMeans(counts(sce)),decreasing=TRUE)])
  # get mean expression across clusters
  cli <- split(seq_len(ncol(sce)), clusters)
  cl.means <- vapply(cli, FUN.VALUE=double(nrow(sce)), FUN=function(x){
    Matrix::rowMeans(counts(sce)[,x,drop=FALSE])
  })
  # grab the top genes in each cluster
  g <- unique(as.numeric(apply(cl.means, 2, FUN=function(x){
    order(x, decreasing=TRUE)[seq_len(nfeatures)]
  })))
  g
}

#' @export
plotDoubletMap <- function(sce, colorBy="enrichment", labelBy="observed", 
                           addSizes=TRUE, col=NULL, column_title="Clusters", 
                           row_title="Clusters", column_title_side="bottom", 
                           na_col="white", ...){
  s <- metadata(sce)$scDblFinder.stats
  s$enrichment <- (s$observed+1)/(s$expected+1)
  colorBy <- match.arg(colorBy, colnames(s))
  labelBy <- match.arg(labelBy, colnames(s))
  comb <- do.call(rbind,strsplit(s$combination,"+",fixed=TRUE))
  colnames(comb) <- paste0("cluster",1:2)
  s <- cbind(comb, s)
  ob <- .castorigins(s, val="observed")
  en <- .castorigins(s, val=colorBy)
  if(colorBy=="enrichment"){
    en <- log2(en)
    colorBy <- "log2\nenrichment"
  }
  if(addSizes){
    sizes <- table(sce$scDblFinder.cluster)
    n <- paste0(colnames(ob), " (", as.numeric(sizes[colnames(ob)]),")")
    colnames(ob) <- row.names(ob) <- colnames(en) <- row.names(en) <- n
    if(is.null(col))
      col <- circlize::colorRamp2(c(min(en,na.rm=TRUE),0,max(en,na.rm=TRUE)),
                                 colors=c("blue","white","red"))
  }else{
    if(is.null(col)) col <- viridisLite::viridis(100)
  }
  
  Heatmap(en, name=colorBy, column_title=column_title, 
          row_title=row_title, column_title_side=column_title_side, 
          col=col, na_col=na_col,
          cell_fun = function(j, i, x, y, width, height, fill){
            if(is.na(ob[i, j])) return(NULL)
            grid.text(as.character(ob[i, j]), x, y, gp=gpar(fontsize=10))
          }, ...)
}