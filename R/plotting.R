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
#' @param fdr logical; whether to plot FDR instead of FPR on the x axis
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
plotROCs <- function(scores, truth, called.class=NULL, nbT=TRUE, dot.size=5, fdr=FALSE){
  truth <- as.integer(as.factor(truth))-1
  scores <- as.data.frame(scores)
  roclist <- lapply(scores, FUN=function(x){
    labels <- truth[order(x, decreasing=TRUE)]
    data.frame(FPR=cumsum(!labels)/sum(!labels),
               FDR=cumsum(!labels)/seq_along(labels),
               TPR=cumsum(labels)/sum(labels))
  })
  methods <- factor( rep(names(roclist), 
                         vapply(roclist, FUN.VALUE=integer(1), 
                                FUN=function(x) length(x$TPR))
                         ),
                     levels=names(roclist) )
  d <- data.frame( method=methods,
                   FPR=unlist(lapply(roclist,FUN=function(x) x$FPR)),
                   TPR=unlist(lapply(roclist,FUN=function(x) x$TPR)),
                   FDR=unlist(lapply(roclist,FUN=function(x) x$FDR)))
  if(fdr){
    p <- ggplot(d, aes(FDR, TPR, colour=method)) + geom_path(size=1.2)
  }else{
    p <- ggplot(d, aes(FPR, TPR, colour=method)) + geom_line(size=1.2)
  }
  if(nbT){
      nbt <- sum(truth)
      d2 <- data.frame( method=names(roclist),
                        FPR=sapply(roclist,FUN=function(x) x$FPR[nbt]),
                        FDR=sapply(roclist,FUN=function(x) x$FDR[nbt]),
                        TPR=unlist(lapply(roclist,FUN=function(x) x$TPR[nbt])) )
      p <- p + geom_point(data=d2, alpha=0.5, size=dot.size)
  }
  if(!is.null(called.class) &&
     length(mm <- intersect(names(called.class), names(roclist)))>0){
      d2 <- data.frame( method=mm,
                        FDR=mapply(r=roclist[mm], cc=called.class[mm], 
                                   FUN=function(r,cc) r$FDR[cc]),
                        FPR=mapply(r=roclist[mm], cc=called.class[mm], 
                                   FUN=function(r,cc) r$FPR[cc]),
                        TPR=mapply(r=roclist[mm], cc=called.class[mm], 
                                   FUN=function(r,cc) r$TPR[cc]) )
      p <- p + geom_point(data=d2, size=dot.size, shape=8, show.legend = FALSE)
  }
  p
}


#' plotDoubletMap
#' 
#' Plots a heatmap of observed versus expected doublets
#' 
#' @param sce A SingleCellExperiment object on which `scDblFinder` has been run.
#' @param colorBy Determines the color mapping. Either "enrichment" (for 
#' log2-enrichment over expectation) or any column of 
#' `metadata(sce)$scDblFinder.stats`
#' @param labelBy Determines the cell labels. Either "enrichment" (for 
#' log2-enrichment over expectation) or any column of
#'  `metadata(sce)$scDblFinder.stats`
#' @param addSizes Logical; whether to add the sizes of clusters to labels
#' @param col The colors scale to use (passed to `ComplexHeatmap::Heatmap`)
#' @param column_title passed to `ComplexHeatmap::Heatmap`
#' @param row_title passed to `ComplexHeatmap::Heatmap`
#' @param column_title_side passed to `ComplexHeatmap::Heatmap`
#' @param na_col color for NA cells
#' @param ... passed to `ComplexHeatmap::Heatmap`
#'
#' @return a Heatmap object
#' 
#' @importFrom ComplexHeatmap Heatmap
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
