#' scDblFinder
#' 
#' Identification of doublets in single-cell RNAseq directly from counts using 
#' overclustering-based generation of artifical doublets.
#'
#' @param sce A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param artificialDoublets The approximate number of artificial doublets to 
#' create. If NULL, will be the maximum of the number of cells or 
#' `5*nbClusters^2`.
#' @param clusters The optional cluster assignments (if omitted, will run 
#' clustering). This is used to make doublets more efficiently. `clusters` 
#' should either be a vector of labels for each cell, or the name of a colData 
#' column of `sce`.
#' @param clust.method The clustering method if `clusters` is not given.
#' @param samples A vector of the same length as cells (or the name of a column 
#' of `colData(sce)`), indicating to which sample each cell belongs. Here, a 
#' sample is understood as being processed independently. If omitted, doublets 
#' will be searched for with all cells together. If given, doublets will be 
#' searched for independently for each sample, which is preferable if they 
#' represent different captures.
#' @param minClusSize The minimum cluster size for `quickCluster`/`overcluster` 
#' (default 50); ignored if `clusters` is given.
#' @param maxClusSize The maximum cluster size for `overcluster`. Ignored if 
#' `clusters` is given. If NA, clustering will be performed using 
#' `quickCluster`, otherwise via `overcluster`. If missing, the default value 
#' will be estimated by `overcluster`.
#' @param nfeatures The number of top features to use (default 1000)
#' @param dims The number of dimensions used to build the network (default 20)
#' @param dbr The expected doublet rate. By default this is assumed to be 1\% 
#' per thousand cells captured (so 4\% among 4000 thousand cells), which is 
#' appropriate for 10x datasets.
#' @param dbr.sd The standard deviation of the doublet rate, defaults to 0.015.
#' @param k Number of nearest neighbors (for KNN graph).
#' @param clust.graph.type Either 'snn' or 'knn'.
#' @param fullTable Logical; whether to return the full table including 
#' artificial doublets, rather than the table for real cells only (default).
#' @param trans The transformation to use before computing the KNN network. The 
#' default, `scran::scaledColRanks`, gave the best result in our hands.
#' @param verbose Logical; whether to print messages and the thresholding plot.
#' @param score Score to use for final classification; either 'weighted' 
#' (default), 'ratio' or 'hybrid' (includes information about library size and 
#' detection rate.
#' @param BPPARAM Used for multithreading when splitting by samples (i.e. when 
#' `samples!=NULL`); otherwise passed to eventual PCA and K/SNN calculations.
#'
#' @return The `sce` object with the following additional colData columns: 
#' `scDblFinder.ratio` (ratio of aritifical doublets among neighbors), 
#' `scDblFinder.weighted` (the ratio of artificial doublets among neighbors 
#' weigted by their distance), `scDblFinder.score` (the final score used,
#' by default the same as `scDblFinder.weighted`), and  `scDblFinder.class` 
#' (whether the cell is called as 'doublet' or 'singlet'). Alternatively, if 
#' `fullTable=TRUE`, a data.frame will be returned with information about real 
#' and artificial cells.
#' 
#' @import SingleCellExperiment Matrix BiocParallel
#' @importFrom SummarizedExperiment colData<- assayNames
#' @importFrom scater normalizeCounts calculatePCA
#' @importFrom dplyr bind_rows
#' @importFrom randomForest randomForest
#' @importFrom methods is
#' @importFrom DelayedArray as.matrix
#' @importFrom BiocNeighbors findKNN
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
                         clust.method=c("louvain","overcluster","fast_greedy"),
                         samples=NULL, minClusSize=min(50,ncol(sce)/5), 
                         maxClusSize=NULL, nfeatures=1000, dims=20, dbr=NULL, 
                         dbr.sd=0.015, k=20, clust.graph.type=c("snn","knn"), 
                         fullTable=FALSE, verbose=is.null(samples), 
                         score=c("weighted","ratio","hybrid"),
                         BPPARAM=SerialParam()
                        ){
  clust.graph.type <- match.arg(clust.graph.type)
  clust.method <- match.arg(clust.method)
  if(!is.null(clusters) && is.character(clusters) && length(clusters)==1){
      if(!(clusters %in% colnames(colData(sce)))) 
          stop("Could not find `clusters` column in colData!")
      clusters <- colData(sce)[[clusters]]
  }
  score <- match.arg(score)
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
                                  clusters=clusters[x], minClusSize=minClusSize, 
                                  maxClusSize=maxClusSize, dims=dims, dbr=dbr, 
                                  dbr.sd=dbr.sd, k=k, clust.method=clust.method,
                                  clust.graph.type=clust.graph.type, 
                                  score=score, nfeatures=nfeatures, 
                                  verbose=FALSE ) )
        fields <- paste0("scDblFinder.",c("weighted","ratio","class","score"))
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
  
  orig <- sce
  if(is.null(clusters)){
      if(verbose) message("Clustering cells...")
      sce <- fastClust( sce, method=clust.method, min.size=minClusSize, 
                        max.size=maxClusSize, nfeatures = nfeatures,
                        graph.type=clust.graph.type, BPPARAM=BPPARAM )
      clusters <- sce$scDblFinder.clusters
  }
  if(nrow(sce)>nfeatures){
      if(verbose) message("Identifying top genes per cluster...")
      # get mean expression across clusters
      cli <- split(seq_len(ncol(sce)), clusters)
      cl.means <- vapply(cli, FUN.VALUE=double(nrow(sce)), FUN=function(x){
          Matrix::rowMeans(counts(sce)[,x,drop=FALSE])
      })
      # grab the top genes in each cluster
      g <- unique(as.numeric(apply(cl.means, 2, FUN=function(x){
          order(x, decreasing=TRUE)[seq_len(nfeatures)]
      })))
      sce <- sce[g,]
  }
  
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
  cli <- split(seq_len(ncol(sce)), clusters)
  
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
    ad <- getArtificialDoublets( as.matrix(counts(sce)), n=artificialDoublets, 
                                 clusters=clusters )
  }else{
    ad <- suppressWarnings( getArtificialDoublets(as.matrix(counts(sce)), 
                                                  n=artificialDoublets, 
                                                  clusters=clusters ) )
  }

  e <- cbind(as.matrix(counts(sce)), ad[row.names(sce),])
  e <- normalizeCounts(e)
  pca <- calculatePCA(e, dims, subset_row=seq_len(nrow(e)))

  # evaluate by library size and non-zero features
  lsizes <- c(colSums(counts(sce)),colSums(ad))
  nfeatures <- c(colSums(counts(sce)>0), colSums(ad>0))
  
  if(verbose) message("Finding KNN...")
  ctype <- factor(rep(c("real","artificial"), c(ncol(sce),ncol(ad))))
  knn <- suppressWarnings(findKNN(pca, k, BPPARAM=BPPARAM))
  
  if(verbose) message("Evaluating cell neighborhoods...")
  knn$type <- matrix(as.numeric(ctype)[knn[[1]]]==1, nrow=nrow(knn[[1]]))
  if(any(knn$distance==0)) knn$distance <- knn$distance + 0.01
  w <- 1/knn$distance

  d <- data.frame( weighted=rowSums(knn$type*w/rowSums(w)),
                   distanceToNearest=knn$distance[,1],
                   ratio=rowSums(knn$type)/k )
  rm(knn)
  rm(pca)
  
  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$type <- ctype

  if(score=="hybrid"){
      rf <- randomForest(x=d[,-ncol(d)],y=d$type, ntree=150, mtry=3)
      fscores <- rf$votes[,"artificial"]
  }else{
      if(score=="ratio"){
          fscores <- d$ratio
      }else{
          fscores <- d$weighted
      }
  }
  d$score <- fscores

  if(verbose) message("Finding threshold...")
  th <- doubletThresholding( fscores, d$type, clusters=clusters, dbr=dbr, 
                             dbr.sd=dbr.sd, do.plot=verbose )
  if(verbose) message("Threshold found:", round(th,3))

  d$classification <- ifelse(fscores >= th, "doublet", "singlet")
  if(fullTable) return(d)
  
  d <- d[seq_len(ncol(orig)),]
  d$type <- NULL
  row.names(d) <- colnames(orig)
  d <- d[colnames(orig),]
  orig$scDblFinder.weighted <- d$weighted
  orig$scDblFinder.ratio <- d$ratio
  orig$scDblFinder.score <- d$score
  orig$scDblFinder.class <- d$classification
  orig
}
