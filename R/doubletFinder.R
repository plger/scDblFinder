#' scDblFinder
#' 
#' Identification of doublets in single-cell RNAseq directly from counts using 
#' overclustering-based generation of artifical doublets.
#'
#' @param sce An object of class `SingleCellExperiment`
#' @param artificialDoublets The approximate number of artificial doublets to create. If
#' NULL, will be the maximum of the number of cells or `3*nbClusters^2`.
#' @param clusters The optional cluster assignments (if ommitted, will run `quickCluster`).
#' This is used to make doublets more efficiently.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if `clusters` is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if `clusters` is
#' given. If NA, clustering will be performed using `quickCluster`, otherwise via
#' `overcluster`. If missing, the default value will be estimated by `overcluster`.
#' @param d The number of dimensions used to build the KNN network (default 10)
#' @param dbr The expected doublet rate. Defaults to 0.025 the basis of the mixology 10x 
#' datasets demuxlet results, where the proportion of demuxlet doublet is estimated 
#' respectively at 0.012 and 0.029 (after accounting adding expected homotypic doublets). 
#' However, note that this value might not be appropriate for all datasets, especially 
#' older 10x datasets, where the proportion of doublets can go up to 0.1.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.02.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param graph.type Either 'snn' or 'knn' (default).
#' @param fullTable Logical; whether to return the full table including artificial 
#' doublets (default FALSE), rather than the table for real cells only.
#' @param trans The transformation to use before computing the KNN network. The default,
#' `scran::scaledColRanks`, gave the best result in our hands.
#' @param verbose Logical; whether to print messages and the thresholding plot (default
#' TRUE).
#' @param BPPARAM Passed to scran; doesn't work so well...
#'
#' @return A data.frame including, for each cells (rows, in the order in which they appear
#' in the `sce` object), the number of neighbors considered, the number of these neighbors
#'  that are from artificial doublets, the ratio of neighbors that are from artificial 
#'  doublets (which represents a quantitative doublet the score), the enrichment in 
#'  artificial doublets, and the classification (after applying the threshold).
#' 
#' @export
scDblFinder <- function( sce, artificialDoublets=NULL, clusters=NULL, minClusSize=50,
                          maxClusSize=NULL, d=10, dbr=0.025, dbr.sd=0.015, k=5, 
                          graph.type=c("knn","snn"), fullTable=FALSE, 
                          trans=c("scran", "rankTrans", "none", "lognorm"), verbose=TRUE, 
                          BPPARAM=SerialParam()){
  library(BiocParallel)
  graph.type <- match.arg(graph.type)
  trans <- match.arg(trans)
  sce2 <- sce
  if(is.null(clusters)){
    # we first simplify the dataset and identify rough clusters
    if(nrow(sce)>3000){
      sce2 <- sce[order(Matrix::rowMeans(counts(sce)), decreasing=TRUE)[1:3000],]
    }
    if(is.null(maxClusSize) || !is.na(maxClusSize)){
      if(verbose) message("Overclustering...")
      clusters <- overcluster(counts(sce2), min.size=minClusSize, max.size=maxClusSize)
      tt <- table(clusters)
      if(verbose) message(length(tt), " clusters of sizes ranging from ", min(tt)," to ", max(tt))
    }else{
      # for reasons of speed, we'll use the fast greedy algorithm if there are many cells
      clust.method <- ifelse(ncol(sce2)>=5000,"igraph","hclust")
      if(verbose) message("Quick ", clust.method, " clustering:")
      clusters <- scran::quickCluster(sce2, min.size=minClusSize, method=clust.method, use.ranks=TRUE, BPPARAM=BPPARAM)
      if(verbose) print(table(clusters))
    }
  }
  maxSameDoublets <- prod(sort(table(clusters), decreasing=TRUE)[1:2]/length(clusters))*length(clusters)*(dbr+2*dbr.sd)
  if(min(table(clusters)) <= maxSameDoublets){
    warning("In light of the expected rate of doublets given, and of the size of the ",
            "clusters, it is possible that some of the smaller clusters are composed of ",
            "doublets of the same type.\n Consider increasing `min.sze`, or breaking down",
            " the larger clusters (e.g. see `?overcluster`).")
  }
  cli <- split(1:ncol(sce2), clusters)
  
  if(nrow(sce2)>2000){
    if(verbose) message("Identifying top genes per cluster...")
    # get mean expression across clusters
    cl.means <- sapply(cli, FUN=function(x) Matrix::rowMeans(counts(sce2[,x])) )
    # grab the top genes in each cluster
    g <- unique(as.numeric(apply(cl.means, 2, FUN=function(x){
      order(x, decreasing=TRUE)[1:1500]
    })))
    sce2 <- sce2[g,]
  }
  
  # get the artificial doublets
  if(is.null(artificialDoublets)){
    artificialDoublets <- max(ncol(sce2), min(40000, 5*length(unique(clusters))^2))
  }
  if(verbose){
    message("Creating ~", artificialDoublets, " artifical doublets...")
    ad <- getArtificialDoublets(counts(sce2), n=artificialDoublets, clusters=clusters)
  }else{
    ad <- suppressWarnings(getArtificialDoublets( counts(sce2), n=artificialDoublets, 
                                                  clusters=clusters))
  }
  
  # build graph and evaluate neigbhorhood
  if(verbose) message("Building ", toupper(graph.type), " graph...")
  qr <- switch(trans,
               scran=scran::scaledColRanks(cbind(counts(sce2), ad)),
               rankTrans=rankTrans(cbind(counts(sce2), ad)),
               none=cbind(counts(sce2), ad),
               norm=logcounts(scater::normalize(SingleCellExperiment(list(counts=cbind(counts(sce2),ad))))) )
  sce2 <- sce2[,intersect(colnames(sce2),colnames(qr))]
  ad <- ad[,intersect(colnames(ad),colnames(qr))]
  if(graph.type=="knn"){
    graph <- suppressWarnings(scran::buildKNNGraph(qr, d=10, k=k, pc.approx=TRUE, BPPARAM=BPPARAM))
  }else{
    graph <- suppressWarnings(scran::buildSNNGraph(qr, BPPARAM=BPPARAM, 
                                type=ifelse(trans %in% c("norm","none"),"number","rank")))
  }
  ctype <- rep(c("real","artificial"), c(ncol(sce2),ncol(ad)))
  
  if(verbose) message("Evaluating cell neighborhoods...")
  graph <- igraph::get.adjlist(graph)
  
  d <- as.data.frame(t(sapply(graph, FUN=function(x){
    x <- as.numeric(x)
    c(length(x), sum(ctype[x]=="artificial"))
  })))

  rm(graph)
  
  colnames(d) <- c("nb.neighbors", "artificial.neighbors")
  d$ratio <- d[,2]/d[,1]
  d$enrichment <- d[,2]/(d[,1]*ncol(ad)/c(ncol(ad)+ncol(sce2)))
  d$type <- ctype

  if(verbose) message("Finding threshold...")
  th <- doubletThresholding(d$ratio, d$type, clusters=clusters, dbr=dbr, dbr.sd=dbr.sd, do.plot=verbose)
  if(verbose) message("Threshold found:", round(th,3))
  d$classification <- ifelse(d$ratio >= th, "doublet", "singlet")
  if(fullTable) return(d)
  
  d <- d[1:ncol(sce2),]
  d$type <- NULL
  row.names(d) <- colnames(sce2)
  d <- d[colnames(sce),]
  sce$scDblFinder.neighbors <- d$nb.neighbors
  sce$scDblFinder.ratio <- d$ratio
  sce$scDblFinder.class <- d$classification
  sce
}
