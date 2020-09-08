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
    origins[w,1] <- as.character(known.origins[w])
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

#' getExpectedDoublets
#'
#' @param x A vector of cluster labels for each cell
#' @param dbr The expected doublet rate.
#' @param only.heterotypic Logical; whether to return expectations only for 
#' heterotypic doublets
#'
#' @return The expected number of doublets of each combination of clusters
#'
#' @examples
#' # random cluster labels
#' cl <- sample(LETTERS[1:4], size=2000, prob=c(.4,.2,.2,.2), replace=TRUE)
#' getExpectedDoublets(cl)
#' @export
getExpectedDoublets <- function(x, dbr=NULL, only.heterotypic=TRUE){
  if(is(x,"SingleCellExperiment")){
    clusters <- x$scDblFinder.clusters
  }else{
    clusters <- x
  }
  clusters <- factor(clusters)
  ncells <- length(clusters)
  if(is.null(dbr)) dbr <- (0.01*ncells/1000)
  if(only.heterotypic) dbr <- .adjust.dbr.homotypy(clusters, dbr)
  cs <- table(clusters)
  eg <- expand.grid(seq_along(cs), seq_along(cs))
  if(only.heterotypic){
    eg <- eg[eg[,1]<eg[,2],]
  }else{
    eg <- eg[eg[,1]<=eg[,2],]
  }
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
  g <- unique(as.numeric(t(apply(cl.means, 2, FUN=function(x){
    order(x, decreasing=TRUE)[seq_len(nfeatures)]
  }))))[seq_len(nfeatures)]
  g
}

#' mockDoubletSCE
#' 
#' Creates a mock random single-cell experiment object with doublets
#'
#' @param ncells A positive integer vector indicating the number of cells per
#' cluster (min 2 clusters)
#' @param ngenes The number of genes to simulate. Ignored if `mus` is given.
#' @param mus A list of cluster averages.
#' @param dbl.rate The doublet rate
#'
#' @return A SingleCellExperiment object, with the colData columns `type` 
#' indicating whether the cell is a singlet or doublet, and `cluster` 
#' indicating from which cluster (or cluster combination) it was simulated.
#' 
#' @export
#' @import SingleCellExperiment stats
#' @examples
#' sce <- mockDoubletSCE()
mockDoubletSCE <- function(ncells=c(200,300), ngenes=200, mus=NULL, 
                           dbl.rate=0.1){
  if(length(ncells)<2)
    stop("ncells should be a positive integer vector of length >=2")
  if(is.null(names(ncells))) names(ncells)<-paste0("cluster",seq_along(ncells))
  if(is.null(mus)){
    mus <- lapply(ncells, FUN=function(x) 2^rnorm(ngenes))
  }else{
    if(!is.list(mus) || length(mus)!=length(ncells) || 
       length(unique(lengths(mus)))!=1)
      stop("If provided, `mus` should be a list of length equal to that of ",
           "`ncells`, with each slot containing a numeric vector of averages ",
           "for each gene.")
    names(mus) <- names(ncells)
  }

  # non-doublets
  counts <- do.call(cbind, lapply(seq_along(mus), FUN=function(i){
    matrix(rpois(ncells[[i]]*ngenes, mus[[i]]), nrow=ngenes)
  }))
  sce <- SingleCellExperiment(list(counts=counts), 
                              colData=data.frame(type="singlet", 
                                           origin=rep(names(ncells), ncells)))
  
  # doublets
  expected <- getExpectedDoublets(sce$origin, dbl.rate, FALSE)
  simdbl <- rpois(length(expected), expected)
  mus <- lapply(strsplit(names(expected),"+",fixed=TRUE), FUN=function(x){
    mus[[x[[1]]]]+mus[[x[[2]]]]
  })
  doublets <- do.call(cbind, mapply(mu=mus, n=simdbl, FUN=function(mu, n){
    matrix(rpois(n*ngenes, mu), nrow=ngenes)
  }))
  sce2 <- SingleCellExperiment( list(counts=doublets), 
                                colData=data.frame(type="doublet", 
                                        origin=rep(names(expected), simdbl)))
  sce <- cbind(sce, sce2)
  
  # set original cluster for homotypic doublets
  sce$cluster <- factor(sce$origin, c(names(ncells),names(expected)))
  names(n) <- n <- levels(sce$cluster)
  n[paste(names(ncells),names(ncells),sep="+")] <- names(ncells)
  levels(sce$cluster) <- as.character(n)
  
  colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  row.names(sce) <- paste0("gene",seq_len(ngenes))
  sce$type <- factor(sce$type, c("singlet","doublet"))
  sce
}


.checkColArg <- function(sce, x, acceptNull=TRUE){
  arg <- deparse(substitute(x))
  if(is.null(x)){
    if(!acceptNull) stop("Missing argument `",arg,"`!")
    return(NULL)
  }
  if(is.character(x) && length(x)==1){
    if(!(x %in% colnames(colData(sce)))) 
      stop("Could not find `", arg, "` column in colData!")
    x <- colData(sce)[[x]]
  }else if(length(x)!=ncol(sce)){
    stop("`",arg,"` should have a length equal to the number of columns in ",
         "`sce`.")
  }
  x
}
