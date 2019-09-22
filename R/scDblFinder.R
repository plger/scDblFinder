#' scDblFinder
#' 
#' Identification of doublets in single-cell RNAseq directly from counts using 
#' overclustering-based generation of artifical doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param artificialDoublets The approximate number of artificial doublets to 
#' create. If NULL, will be the maximum of the number of cells or 
#' `3*nbClusters^2`.
#' @param clusters The optional cluster assignments (if ommitted, will run 
#' `quickCluster`). This is used to make doublets more efficiently.
#' @param samples A vector of the same length as cells (or the name of a column of 
#' `colData(sce)`), indicating to which sample each cell belongs. Here, a sample
#' is understood as being processed independently. If ommitted, doublets will be
#' searched for with all cells together. If given, doublets will be searched for 
#' independently for each sample, which is preferable if they represent different
#' captures.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if `clusters` is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if 
#' `clusters` is given. If NA, clustering will be performed using `quickCluster`,
#' otherwise via `overcluster`. If missing, the default value will be estimated 
#' by `overcluster`.
#' @param d The number of dimensions used to build the KNN network (default 10)
#' @param dbr The expected doublet rate. By default this is assumed to be 
#' 0.1*nbCells percent, which is appropriate for 10x datasets.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.015.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param graph.type Either 'snn' or 'knn' (default).
#' @param fullTable Logical; whether to return the full table including artificial 
#' doublets (default FALSE), rather than the table for real cells only.
#' @param trans The transformation to use before computing the KNN network. The 
#' default, `scran::scaledColRanks`, gave the best result in our hands.
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when 
#' `samples!=NULL`); otherwise passed to scran (which doesn't work so well)...
#'
#' @return The `sce` object with the following additional colData columns: 
#' `scDblFinder.neighbors` (number of neighbors considered), `scDblFinder.ratio`
#' (ratio of aritifical doublets among neighbors), and `scDblFinder.class` 
#' (whether the cell is called as 'doublet' or 'singlet'). Alternatively, if 
#' `fullTable=TRUE`, a data.frame will be returned with information about real 
#' and artificial cells.
#' 
#' @import SingleCellExperiment scran Matrix BiocParallel
#' @importFrom dplyr bind_rows
#' @importFrom testthat test_that expect_equal expect_that is_a
#' 
#' @examples
#' library(SingleCellExperiment)
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' sce <- SingleCellExperiment( list(counts=m) )
#' sce <- scDblFinder(sce, minClusSize=2, maxClusSize=20, verbose=FALSE)
#' 
#' @export
scDblFinder <- function( sce, artificialDoublets=NULL, clusters=NULL, 
                         samples=NULL, minClusSize=50, maxClusSize=NULL, d=10, 
                         dbr=NULL, dbr.sd=0.015, k=5, graph.type=c("knn","snn"), 
                         fullTable=FALSE, 
                         trans=c("scran", "rankTrans", "none", "lognorm"), 
                         verbose=is.null(samples), BPPARAM=SerialParam()){
  graph.type <- match.arg(graph.type)
  trans <- match.arg(trans)
  expect_that(sce, is_a("SingleCellExperiment"))
  if( !("counts" %in% assayNames(sce)) ) 
      stop("`sce` should have an assay named 'counts'")
  if(!is.null(samples)){
    test_that( "Samples vector matches number of cells", 
               expect_equal(ncol(sce), length(samples)) )
    # splitting by samples
    if(length(samples)==1 && samples %in% colnames(colData(sce)))
        samples <- colData(sce)[[samples]]
    test_that( "Samples vector matches number of cells", 
               expect_equal(ncol(sce), length(samples)) )
    if(fullTable) warning("`fullTable` param ignored when splitting by samples")
    if(verbose) message("`verbose` param ignored when splitting by samples.")
    cs <- split(seq_along(samples), samples, drop=TRUE)
    CD <- bind_rows(bplapply(cs, BPPARAM=BPPARAM, FUN=function(x){ 
        if(!is.null(clusters) && length(clusters)>1) clusters <- clusters[x]
        x <- colData( scDblFinder(sce[,x], artificialDoublets=artificialDoublets, 
                                  clusters=clusters, minClusSize=minClusSize, 
                                  maxClusSize=maxClusSize, d=d, dbr=dbr, 
                                  dbr.sd=dbr.sd, k=k, graph.type=graph.type, 
                                  trans=trans, verbose=FALSE) )
        as.data.frame(x[,paste0("scDblFinder.",c("neighbors","ratio","class"))])
    }))
    for(f in colnames(CD)) colData(sce)[[f]] <- CD[unlist(cs),f]
    return(sce)
  }
  if(is.null(dbr)){
      ## dbr estimated as for chromium data, 1% per 1000 cells captured:
      dbr <- (0.01*ncol(sce)/1000)
  }
  if(is.null(clusters)) clusters <- .getClusters(sce, maxClusSize, minClusSize)
  maxSameDoublets <- length(clusters)*(dbr+2*dbr.sd) * 
      prod(sort(table(clusters), decreasing=TRUE)[1:2] / length(clusters))
  if(min(table(clusters)) <= maxSameDoublets){
    warning("In light of the expected rate of doublets given, and of the size ",
            "of the clusters, it is possible that some of the smaller clusters",
            " are composed of doublets of the same type.\n ",
            "Consider increasing `min.sze`, or breaking down the larger ",
            "clusters (e.g. see `?overcluster`).")
  }
  cli <- split(1:ncol(sce), clusters)
  
  sce2 <- sce
  if(is.null(colnames(sce2)))
      colnames(sce2) <- paste0("cell",seq_len(ncol(sce2)))
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
    artificialDoublets <- max( ncol(sce2), 
                               min(40000, 5*length(unique(clusters))^2))
  }
  if(verbose){
    message("Creating ~", artificialDoublets, " artifical doublets...")
    ad <- getArtificialDoublets( counts(sce2), n=artificialDoublets, 
                                 clusters=clusters )
  }else{
    ad <- suppressWarnings( getArtificialDoublets(counts(sce2), 
                                                  n=artificialDoublets, 
                                                  clusters=clusters ) )
  }
  
  # build graph and evaluate neigbhorhood
  if(verbose) message("Building ", toupper(graph.type), " graph...")
  qr <- switch(trans,
               scran=scran::scaledColRanks(cbind(counts(sce2), ad)),
               rankTrans=rankTrans(cbind(counts(sce2), ad)),
               none=cbind(counts(sce2), ad),
               norm=logcounts(scater::normalize(SingleCellExperiment(
                   list(counts=cbind(counts(sce2),ad)))))
               )
  sce2 <- sce2[,intersect(colnames(sce2),colnames(qr))]
  ad <- ad[,intersect(colnames(ad),colnames(qr))]
  if(graph.type=="knn"){
    graph <- suppressWarnings(buildKNNGraph( qr, d=10, k=k,
                                             BPPARAM=BPPARAM))
  }else{
    graph <- suppressWarnings(buildSNNGraph(qr, BPPARAM=BPPARAM, 
                              type=ifelse( trans %in% c("norm","none"),
                                           "number","rank") ))
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
  th <- doubletThresholding( d$ratio, d$type, clusters=clusters, dbr=dbr, 
                             dbr.sd=dbr.sd, do.plot=verbose )
  if(verbose) message("Threshold found:", round(th,3))

  d$classification <- ifelse(d$ratio >= th, "doublet", "singlet")
  if(fullTable) return(d)
  
  d <- d[seq_len(ncol(sce2)),]
  d$type <- NULL
  row.names(d) <- colnames(sce2)
  d <- d[colnames(sce),]
  sce$scDblFinder.neighbors <- d$nb.neighbors
  sce$scDblFinder.ratio <- d$ratio
  sce$scDblFinder.class <- d$classification
  sce
}


.getClusters <- function(sce, maxClusSize, minClusSize, ngenes=3000, verbose=TRUE){
    # we first simplify the dataset and identify rough clusters
    o <- order(Matrix::rowMeans(counts(sce)), decreasing=TRUE)
    sce <- sce[o[seq_len(min(nrow(sce),ngenes))],]
    if(is.null(maxClusSize) || !is.na(maxClusSize)){
        if(verbose) message("Overclustering...")
        clusters <- overcluster( counts(sce), 
                                 min.size=minClusSize, max.size=maxClusSize)
    }else{
        # for reasons of speed, with many cells we use the fast greedy algorithm
        clust.method <- ifelse(ncol(sce)>=5000,"igraph","hclust")
        if(verbose) message("Quick ", clust.method, " clustering:")
        clusters <- scran::quickCluster( sce, min.size=minClusSize, 
                                         method=clust.method, use.ranks=TRUE, 
                                         BPPARAM=BPPARAM)
    }
    if(verbose) print(table(clusters))
    clusters
}