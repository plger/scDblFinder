# get random cross-cluster pairs of cells from a cluster assignment vector
.getCellPairs <- function(clusters, n=1000){
  cli <- split(1:length(clusters), clusters)
  ca <- expand.grid(unique(clusters), unique(clusters))
  ca <- ca[which(ca[,1]!=ca[,2]),]
  ca <- do.call(rbind, lapply(1:nrow(ca), n=ceiling(n/nrow(ca)), FUN=function(i,n){ 
    cbind( sample(cli[[ca[i,1]]],size=n,replace=TRUE),
           sample(cli[[ca[i,2]]],size=n,replace=TRUE) )
  }))
  ca[!duplicated(ca),]
}  

# creates within-cluster meta-cells from a count matrix
.getMetaCells <- function(x, clusters, n.meta.cells=20, meta.cell.size=20){
  cli <- split(1:length(clusters), clusters)
  meta <- unlist(lapply(cli, FUN=function(x){
    lapply(1:n.meta.cells, FUN=function(y){
      sample(x,min(ceiling(0.6*length(x)),meta.cell.size),replace=FALSE)
    })
  }), recursive=FALSE)
  meta <- sapply(meta, FUN=function(y){ Matrix::rowMeans(x[,y]) })
  colnames(meta) <- paste0("metacell.",1:ncol(meta))
  meta
}
