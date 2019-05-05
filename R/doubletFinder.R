#' plDoubletFinder
#' 
#' An approach similar to `DoubletFinder` for identifying doublets, but much 
#' simpler/faster and independent from normalization (executed directly on counts).
#'
#' @param sce An object of class `SingleCellExperiment`
#' @param artificialDoublets The approximate number of artificial doublets to create. If
#' NULL, will be the maximum of the number of cells or `3*nbClusters^2`.
#' @param clusters The optional cluster assignments (if ommitted, will run `quickCluster`).
#' This is used to make doublets more efficiently.
#' @param minClusSize The minimum cluster size for `quickCluster` (default 50)
#' @param d The number of dimensions used to build the KNN network (default 10)
#' @param dbr The expected doublet rate, defaults to 0.021 on the basis of the 
#' on the basis of the mixology 10x datasets results, where the proportion of
#' demuxlet doublet is estimated respectively at 0.012 and 0.029 (after accounting
#' for homotypic doublets).
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.01.
#' @param graph.type Either 'snn' or 'knn' (default).
#' @param fullTable Logical; whether to return the full table including artificial 
#' doublets (default FALSE), rather than the table for real cells only.
#' @param verbose Logical; whether to print messages and the thresholding plot (default
#' TRUE).
#' @param ... arguments passed to `scran::buildSNNGraph`; can for instance be used to
#' pass `BPPARAM` and speed up the graph building.
#'
#' @return A data.frame including, for each cells, the number of neighbors considered,
#' the number of these neighbors that are from artificial doublets, the ratio of 
#' neighbors that are from artificial doublets, the enrichment in artificial doublets,
#'  and the classification (after applying the threshold).
#' 
#' @export
plDoubletFinder <- function( sce, artificialDoublets=NULL, clusters=NULL, minClusSize=50,
                             d=10, dbr=0.021, dbr.sd=0.01, graph.type=c("knn","snn"), 
                             fullTable=FALSE, verbose=TRUE, ...){
  graph.type <- match.arg(graph.type)
  if(is.null(clusters)){
    # we first simplify the dataset and identify rough clusters
    # for reasons of speed, we'll use the fast greedif there are 
    clust.method <- ifelse(ncol(sce)>=5000,"igraph","hclust")
    if(verbose) message("Quick ", clust.method, " clustering:")
    if(nrow(sce)>3000){
      sce <- sce[order(Matrix::rowMeans(counts(sce)), decreasing=TRUE)[1:3000],]
    }
    clusters <- scran::quickCluster(sce, min.size=minClusSize, method=clust.method)
    if(verbose) print(table(clusters))
  }
  cli <- split(1:ncol(sce), clusters)
  
  if(verbose) message("Identifying top genes per cluster...")
  # get mean expression across clusters
  cl.means <- sapply(cli, FUN=function(x) Matrix::rowMeans(counts(sce[,x])) )
  # grab the top genes in each cluster
  g <- unique(as.numeric(apply(cl.means, 2, FUN=function(x){
    order(x, decreasing=TRUE)[1:1500]
  })))
  sce <- sce[g,]
  cells <- colnames(sce)
  
  if(is.null(artificialDoublets)){
    artificialDoublets <- max(ncol(sce), 3*length(unique(clusters))^2)
  }
  if(verbose) message("Creating ~", artificialDoublets, " artifical doublets...")
  ad <- getArtificialDoublets(counts(sce), n=artificialDoublets, clusters=clusters)
  
  if(verbose) message("Building ", toupper(graph.type), " graph...")
  if(graph.type=="knn"){
    qr <- scran::scaledColRanks(cbind(counts(sce), ad))
    sce <- sce[,intersect(colnames(sce),colnames(qr))]
    ad <- ad[,intersect(colnames(ad),colnames(qr))]
    graph <- scran::buildKNNGraph(qr, d=10, pc.approx=TRUE, ...)
  }else{
    graph <- suppressWarnings(scran::buildSNNGraph(cbind(counts(sce), ad), type="rank", ...))
  }
  ctype <- rep(c("real","artificial"), c(ncol(sce),ncol(ad)))
  
  if(verbose) message("Evaluating cell neighborhoods...")
  graph <- igraph::get.adjlist(graph)
  d <- t(sapply(graph,FUN=function(x){
    x <- as.numeric(x)
    c(length(x), sum(ctype[x]=="artificial"))
  }))
  rm(graph)
  
  d <- as.data.frame(d)
  colnames(d) <- c("nb.neighbors", "artificial.neighbors")
  d$ratio <- d[,2]/d[,1]
  d$enrichment <- d[,2]/(d[,1]*ncol(ad)/c(ncol(ad)+ncol(sce)))
  d$type <- ctype
  
  if(verbose) message("Finding threshold...")
  # adjust expected doublet proportion for expected homotypic doublets
  homotypic.prop <- sum((table(clusters)/length(clusters))^2)
  dbr <- dbr*(1-homotypic.prop)
  f1 <- ecdf(d$ratio[which(d$type!="artificial")])
  f2 <- ecdf(d$ratio[which(d$type=="artificial")])
  # accuracy
  accfn <- function(x){ (1-f1(x))+max(0,f2(x)-homotypic.prop) }
  accfn2 <- function(x){ (1-f1(x))+f2(x) }
  # deviation from expected doublet proportion
  dbr.dev <- function(x){
    p <- sum(d$ratio>=x & d$type=="real")/ncol(sce)
    abs(pcauchy(abs(dbr-p), scale=dbr.sd)*2-1)^2
  }
  totfn <- function(x){
    accfn(x)+dbr.dev(x)
  }
  th <- optimize(totfn, c(0,1), maximum=FALSE)$minimum
  if(verbose){
    # we plot the thresholding data
    x <- 1:99/100
    acc <- sapply(x, accfn)
    plot(x, acc, type="l", lty="dashed", col="blue", ylim=c(0,max(acc)), lwd=2,
         main="Thresholding", ylab="Error rate", 
         xlab="Ratio of artificial doublets in neighborhood")
    lines(x, sapply(x, accfn2), col="blue", lwd=2)
    lines(x, sapply(x, dbr.dev), col="red", lwd=2)
    lines(x, sapply(x, FUN=function(x) sum(d$ratio>=x & d$type=="real")/ncol(sce)), col="darkgrey")
    abline(v=th, lty="dashed")
    leg <- c( "Accuracy in classifying real cells vs artificial doublets",
              "Accuracy in classifying artificial doublets, adjusted for homotypy",
              "Deviation from expected doublet rate",
              "Proportion of real cells considered doublet" )
    legend("top", col=c("blue","blue", "red","darkgrey"), lty=c(1,2,1,1), lw=2, legend=leg)
  }
  if(verbose) message("Threshold found:", round(th,3))
  d$classification <- ifelse(d$ratio >= th, "doublet", "singlet")
  if(fullTable) return(d)
  d <- d[1:ncol(sce),]
  d$type <- NULL
  row.names(d) <- colnames(sce)
  d <- d[cells,]
  d
}
