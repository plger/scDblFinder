#' scDblFinder
#' 
#' Identification of doublets in single-cell RNAseq directly from counts using 
#' overclustering-based generation of artifical doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param artificialDoublets The approximate number of artificial doublets to 
#' create. If NULL, will be the maximum of the number of cells or 
#' `3*nbClusters^2`.
#' @param clusters The optional cluster assignments (if omitted, will run 
#' `quickCluster`). This is used to make doublets more efficiently.
#' @param samples A vector of the same length as cells (or the name of a column 
#' of `colData(sce)`), indicating to which sample each cell belongs. Here, a 
#' sample is understood as being processed independently. If omitted, doublets 
#' will be searched for with all cells together. If given, doublets will be 
#' searched for independently for each sample, which is preferable if they 
#' represent different captures.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if `clusters` is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if 
#' `clusters` is given. If NA, clustering will be performed using `quickCluster`,
#' otherwise via `overcluster`. If missing, the default value will be estimated 
#' by `overcluster`.
#' @param d The number of dimensions used to build the KNN network (default 10)
#' @param dbr The expected doublet rate. By default this is assumed to be 1\% 
#' per thousand cells captured (so 4\% among 4000 thousand cells), which is 
#' appropriate for 10x datasets.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.015.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param graph.type Either 'snn' or 'knn' (default).
#' @param fullTable Logical; whether to return the full table including artificial 
#' doublets (default FALSE), rather than the table for real cells only.
#' @param trans The transformation to use before computing the KNN network. The 
#' default, `scran::scaledColRanks`, gave the best result in our hands.
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param hybridScore Logical; whether to combine similarity to artificial 
#' doublets with information about library size and detection rate (default)
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when 
#' `samples!=NULL`); otherwise passed to scran (which doesn't work so well)...
#'
#' @return The `sce` object with the following additional colData columns: 
#' `scDblFinder.neighbors` (number of neighbors considered), `scDblFinder.ratio`
#' (ratio of aritifical doublets among neighbors), `scDblFinder.score` (an 
#' integrated doublet score, if `hybridScore=TRUE`), and `scDblFinder.class` 
#' (whether the cell is called as 'doublet' or 'singlet'). Alternatively, if 
#' `fullTable=TRUE`, a data.frame will be returned with information about real 
#' and artificial cells.
#' 
#' @import SingleCellExperiment scran Matrix BiocParallel
#' @importFrom dplyr bind_rows
#' @importFrom randomForest randomForest
#' @importFrom methods is
#' 
#' @examples
#' library(SingleCellExperiment)
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(50,x) ) )
#' sce <- SingleCellExperiment( list(counts=m) )
#' sce <- scDblFinder(sce, verbose=FALSE)
#' table(sce$scDblFinder.class)
#' 
#' @export
scDblFinder <- function( sce, artificialDoublets=NULL, clusters=NULL, 
                         samples=NULL, minClusSize=min(50,ncol(sce)/5), 
                         maxClusSize=NULL, d=10, dbr=NULL, dbr.sd=0.015, k=5, 
                         graph.type=c("knn","snn"), fullTable=FALSE, 
                         trans=c("rankTrans", "scran", "none", "lognorm"), 
                         verbose=is.null(samples), hybridScore=TRUE,
                         BPPARAM=SerialParam()
                        ){
  graph.type <- match.arg(graph.type)
  trans <- match.arg(trans)
  stopifnot(is(sce, "SingleCellExperiment"))
  if( !("counts" %in% assayNames(sce)) ) 
      stop("`sce` should have an assay named 'counts'")
  if(!is.null(samples)){
    # splitting by samples
    if(length(samples)==1 && samples %in% colnames(colData(sce)))
        samples <- colData(sce)[[samples]]
    if(ncol(sce)!=length(samples))
        stop("Samples vector matches number of cells")
    if(fullTable) warning("`fullTable` param ignored when splitting by samples")
    if(verbose) message("`verbose` param ignored when splitting by samples.")
    cs <- split(seq_along(samples), samples, drop=TRUE)
    CD <- bind_rows(bplapply(cs, BPPARAM=BPPARAM, FUN=function(x){ 
        if(!is.null(clusters) && length(clusters)>1) clusters <- clusters[x]
        x <- colData( scDblFinder(sce[,x], artificialDoublets=artificialDoublets, 
                                  clusters=clusters, minClusSize=minClusSize, 
                                  maxClusSize=maxClusSize, d=d, dbr=dbr, 
                                  dbr.sd=dbr.sd, k=k, graph.type=graph.type, 
                                  trans=trans, verbose=FALSE, 
				                  hybridScore=hybridScore) )
	fields <- paste0("scDblFinder.",c("neighbors","ratio","class"))
	if(hybridScore) fields <- c(fields, "scDblFinder.score")
        as.data.frame(x[,fields])
    }))
    for(f in colnames(CD)) colData(sce)[[f]] <- CD[unlist(cs),f]
    return(sce)
  }
  if(ncol(sce)<100) warning("scDblFinder might not work well with very low 
numbers of cells.")

  if(is.null(dbr)){
      ## dbr estimated as for chromium data, 1% per 1000 cells captured:
      dbr <- (0.01*ncol(sce)/1000)
  }
  if(is.null(clusters)) clusters <- .getClusters( sce, maxClusSize, minClusSize, 
                                                  BPPARAM=BPPARAM)
  if(length(unique(clusters)) == 1) stop("Only one cluster generated")
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
  
  if(is.null(colnames(sce)))
      colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  if(is.null(row.names(sce)))
      row.names(sce) <- paste0("f",seq_len(nrow(sce)))

  # get the artificial doublets
  if(is.null(artificialDoublets)){
    artificialDoublets <- max( ncol(sce), 
                               min(40000, 5*length(unique(clusters))^2))
  }
  if(verbose){
    message("Creating ~", artificialDoublets, " artifical doublets...")
    ad <- getArtificialDoublets( counts(sce), n=artificialDoublets, 
                                 clusters=clusters )
  }else{
    ad <- suppressWarnings( getArtificialDoublets(counts(sce), 
                                                  n=artificialDoublets, 
                                                  clusters=clusters ) )
  }
  
  # evaluate by library size and non-zero features
  lsizes <- c(colSums(counts(sce)),colSums(ad))
  nfeatures <- c(colSums(counts(sce)>0), colSums(ad>0))
  
  sce2 <- sce
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
  ad <- ad[row.names(sce2),]
  e <- cbind(counts(sce2), ad)

  # build graph and evaluate neigbhorhood
  if(verbose) message("Building ", toupper(graph.type), " graph...")
  qr <- switch(trans,
               scran=scran::scaledColRanks(e),
               rankTrans=rankTrans(e),
               none=e,
               norm=logcounts(scater::normalize(
                   SingleCellExperiment(list(counts=e)) ))
               )
  rm(e)
  
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
  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$type <- ctype

  if(hybridScore){
      rf <- randomForest(x=d[,-ncol(d)],y=as.factor(d$type), ntree=300, mtry=3)
      d$score <- rf$votes[,"artificial"]
      fscores <- d$score
  }else{
      fscores <- d$ratio
  }

  if(verbose) message("Finding threshold...")
  th <- doubletThresholding( fscores, d$type, clusters=clusters, dbr=dbr, 
                             dbr.sd=dbr.sd, do.plot=verbose )
  if(verbose) message("Threshold found:", round(th,3))

  d$classification <- ifelse(fscores >= th, "doublet", "singlet")
  if(fullTable) return(d)
  
  d <- d[seq_len(ncol(sce2)),]
  d$type <- NULL
  row.names(d) <- colnames(sce2)
  d <- d[colnames(sce),]
  sce$scDblFinder.neighbors <- d$nb.neighbors
  sce$scDblFinder.ratio <- d$ratio
  if(hybridScore) sce$scDblFinder.score <- d$score
  sce$scDblFinder.class <- d$classification
  sce
}


#' @import SingleCellExperiment scran Matrix BiocParallel
.getClusters <- function(sce, maxClusSize, minClusSize, ngenes=3000, 
                         verbose=TRUE, BPPARAM=SerialParam()){
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
