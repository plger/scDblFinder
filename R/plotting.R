#' plotDoubletMap
#' 
#' Plots a heatmap of observed versus expected doublets. 
#' Requires the `ComplexHeatmap` package.
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
#' @export
#' @importFrom stats aggregate
plotDoubletMap <- function(sce, colorBy="enrichment", labelBy="observed", 
                           addSizes=TRUE, col=NULL, column_title="Clusters", 
                           row_title="Clusters", column_title_side="bottom", 
                           na_col="white", ...){
  if(is.data.frame(sce)){
    s <- sce
  }else{
    s <- metadata(sce)$scDblFinder.stats
  }
  if(is.null(s)) stop("Could not find doublet metadata. Was scDblFinder run?")
  if(isMultiSample <- is(s,"list")) s <- dplyr::bind_rows(s, .id="sample")
  s$enrichment <- log2((s$observed+1)/(s$expected+1))
  colorBy <- match.arg(colorBy, colnames(s))
  labelBy <- match.arg(labelBy, colnames(s))
  comb <- do.call(rbind,strsplit(s$combination,"+",fixed=TRUE))
  colnames(comb) <- paste0("cluster",1:2)
  s <- cbind(comb, s)
  if(isMultiSample)
    ag <- aggregate(s[,c(labelBy,colorBy)], by=s[,1:2], na.rm=TRUE, FUN=mean)
  doag <- function(x) isMultiSample && !(x %in% c("observed","expected"))
  ob <- .castorigins(switch(as.character(doag(labelBy)),
                           "TRUE"=ag, "FALSE"=s), val=labelBy)
  en <- .castorigins(switch(as.character(doag(colorBy)),
                           "TRUE"=ag, "FALSE"=s), val=colorBy)
  if(colorBy=="enrichment"){
    colorBy <- "log2\nenrichment"
    if(is.null(col))
        col <- circlize::colorRamp2(c(min(en,na.rm=TRUE),0,max(en,na.rm=TRUE)),
                                     colors=c("blue","white","red"))
  }else if(is.null(col)){
    col <- viridisLite::viridis(100)
  }
  if(doag(colorBy)) colorBy <- paste0("mean\n", colorBy)
  if(addSizes && !is.null(sce$scDblFinder.cluster)){
    sizes <- table(sce$scDblFinder.cluster)
    n <- paste0(colnames(ob), " (", as.numeric(sizes[colnames(ob)]),")")
    colnames(ob) <- row.names(ob) <- colnames(en) <- row.names(en) <- n
  }
  ComplexHeatmap::Heatmap(en, name=colorBy, column_title=column_title, 
          row_title=row_title, column_title_side=column_title_side, 
          col=col, na_col=na_col,
          cell_fun = function(j, i, x, y, width, height, fill){
            if(is.na(ob[i, j])) return(NULL)
            grid::grid.text(as.character(ob[i, j]), x, y, 
                            gp=grid::gpar(fontsize=10))
          }, ...)
}
